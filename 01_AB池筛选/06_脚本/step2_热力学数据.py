# -*- coding: utf-8 -*-
"""
规则1(续): 母池热力学数据补齐

对母池 1,860 项统一补齐(1610 项已有结果自动复用, 新增 ~446 项现场计算):
  (a) 电化学窗口 V_red / V_ox
      方法: 巨势相图(GrandPotentialPhaseDiagram), 开放元素 Li, 电位窗口 [-0.5, 7.0] V,
      0.25 V 粗扫 + 边界二分(7 次, 精度 ~0.005 V); 判据 e_above_hull <= 自身Ehull + 50 meV。
      归一化: 相图能量 per atom(MP S3 thermo manifest, 已含 MP2020 校正, 不二次校正)。
  (b) 对 Li 反应能 Er(Li,X)
      方法: 封闭组成相图两相混合反应(InterfacialReactivity, norm=True,
      use_hull_energy=False), 归一化为 eV/atom(混合原子)。化学空间 = X 的元素 + Li。
      说明: 封闭相图 vs Li 与 0 V 巨势界面反应(GPIR)等价(μ_Li = μ0), 此处用封闭法,
      对母池统一记录计算方式为 "closed_PD_vs_Li(=GPIR@0V)"。
  (c) 对目标 SSE 反应能 Er(X,SSE), SSE = Li2ZrCl6_beta 与 Li3YCl6(VASP 自算+Cl修正)
      方法: 封闭组成相图两相混合反应(InterfacialReactivity, norm=True), eV/atom。
      化学空间 > 6 元跳过(相图维度爆炸, 记 chemspace>6 → 待核)。

checkpoint 逐条落盘, 中断重跑自动续。
输出: ../02_热力学数据/窗口_checkpoint.jsonl, 反应能_checkpoint.jsonl
"""
import pickle, os, csv, json, time, warnings
warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.abspath(__file__))
WS = os.path.dirname(ROOT)
TOP = os.path.dirname(WS)
DATA = os.path.join(WS, "02_热力学数据")
os.makedirs(DATA, exist_ok=True)

from pymatgen.core import Element, Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram, GrandPotentialPhaseDiagram, GrandPotPDEntry
from pymatgen.analysis.interface_reactions import InterfacialReactivity

CKPT_WIN = os.path.join(DATA, "窗口_checkpoint.jsonl")
CKPT_RXN = os.path.join(DATA, "反应能_checkpoint.jsonl")


def main():
    raw = pickle.load(open(os.path.join(TOP, "multi_screen_Jeon2025", "02_数据获取", "pd_cache", "mp_gga_raw_entries.pkl"), "rb"))
    hal = pickle.load(open(os.path.join(TOP, "halide_vasp", "halide_entries_corrected.pkl"), "rb"))
    SSE_REAL = {"Li2ZrCl6": hal["Li2ZrCl6_beta"], "Li3YCl6": hal["Li3YCl6"]}
    entries = raw + list(SSE_REAL.values())
    by_formula = {}
    for e in entries:
        by_formula.setdefault(e.composition.reduced_formula, []).append(e)
    mu0 = min((e for e in entries if e.composition.reduced_formula == "Li"),
              key=lambda e: e.energy_per_atom).energy_per_atom

    pool = list(csv.DictReader(open(os.path.join(WS, "01_母池", "母池_无带隙门槛.csv"), encoding="utf-8-sig")))

    # ---------- 窗口 ----------
    def stable_at(rel, target, mu, tol):
        gpd = GrandPotentialPhaseDiagram(rel, {Element("Li"): mu})
        ge = GrandPotPDEntry(target, {Element("Li"): mu})
        _, ea = gpd.get_decomp_and_e_above_hull(ge, allow_negative=True)
        return ea <= tol

    def window_fast(rel, target, own, tol_extra=0.0501, dv=0.25):
        tol = own + tol_extra
        Vs = [round(-0.5 + i * dv, 3) for i in range(int(7.5 / dv) + 1)]
        flags = [stable_at(rel, target, mu0 - V, tol) for V in Vs]
        if not any(flags):
            return None
        first = flags.index(True)
        last = len(flags) - 1 - flags[::-1].index(True)

        def refine(lo, hi):
            for _ in range(7):
                mid = (lo + hi) / 2
                if stable_at(rel, target, mu0 - mid, tol):
                    lo = mid
                else:
                    hi = mid
            return (lo + hi) / 2

        vred = Vs[first] if first == 0 else refine(Vs[first - 1], Vs[first])
        vox = Vs[last] if last == len(Vs) - 1 else refine(Vs[last], Vs[last + 1])
        return round(vred, 2), round(vox, 2)

    done_w = set()
    if os.path.exists(CKPT_WIN):
        for line in open(CKPT_WIN, encoding="utf-8"):
            try:
                done_w.add(json.loads(line)["mp_id"])
            except Exception:
                pass
    # 从旧 1610 项 checkpoint 一次性播种(同一算法与能量源, 格式完全一致)
    OLD_WIN = os.path.join(TOP, "multi_screen_Jeon2025", "02_数据获取", "pd_cache", "window_checkpoint_v3.jsonl")
    n_seed = 0
    if os.path.exists(OLD_WIN):
        pool_ids = {r["mp_id"] for r in pool}
        with open(CKPT_WIN, "a", encoding="utf-8") as fseed:
            for line in open(OLD_WIN, encoding="utf-8"):
                try:
                    d = json.loads(line)
                except Exception:
                    continue
                if d["mp_id"] in pool_ids and d["mp_id"] not in done_w:
                    d["method"] = "GPD_openLi_mu=mu0-V_scan[-0.5,7.0]_tol=Ehull+50meV(reused_v3)"
                    fseed.write(json.dumps(d, ensure_ascii=False) + "\n")
                    done_w.add(d["mp_id"]); n_seed += 1
    print(f"[窗口] 旧checkpoint播种 {n_seed} 项, 待新算 {len(pool)-len(done_w)} 项")
    fw = open(CKPT_WIN, "a", encoding="utf-8")
    t0 = time.time(); n = 0
    for r in pool:
        if r["mp_id"] in done_w:
            continue
        n += 1
        try:
            cands = by_formula.get(r["formula"], [])
            if not cands:
                res = {"mp_id": r["mp_id"], "status": "no_entry_in_manifest"}
            else:
                target = min(cands, key=lambda e: e.energy_per_atom)
                els = {el.symbol for el in target.composition.elements}
                rel = [e for e in raw if all(el.symbol in els | {"Li"} for el in e.composition.elements)]
                pd = PhaseDiagram(rel)
                _, own = pd.get_decomp_and_e_above_hull(target)
                w = window_fast(rel, target, own)
                if w:
                    res = {"mp_id": r["mp_id"], "V_red": w[0], "V_ox": w[1],
                           "window": round(w[1] - w[0], 2), "Ehull_own_meV": round(own * 1000, 1),
                           "method": "GPD_openLi_mu=mu0-V_scan[-0.5,7.0]_tol=Ehull+50meV",
                           "status": "ok"}
                else:
                    res = {"mp_id": r["mp_id"], "status": "never_stable",
                           "method": "GPD_openLi_mu=mu0-V_scan[-0.5,7.0]"}
        except Exception as ex:
            res = {"mp_id": r["mp_id"], "status": f"error:{type(ex).__name__}:{str(ex)[:100]}"}
        fw.write(json.dumps(res, ensure_ascii=False) + "\n"); fw.flush()
        if n % 50 == 0:
            el = time.time() - t0
            print(f"[窗口] {n} done {el:.0f}s ({el/n:.2f}s/项)", flush=True)
    fw.close()
    print("[窗口] DONE", time.time() - t0)

    # ---------- 反应能(vs Li + vs 两个真实 SSE) ----------
    # vs-Li 从旧 1610 项 checkpoint 复用(同一封闭 PD 方法与能量源), 只新算缺失项
    OLD_RXN = os.path.join(TOP, "multi_screen_Jeon2025", "02_数据获取", "pd_cache", "reaction_checkpoint_v2.jsonl")
    old_rxn = {}
    if os.path.exists(OLD_RXN):
        for line in open(OLD_RXN, encoding="utf-8"):
            try:
                d = json.loads(line)
                if d.get("status") == "ok" and d.get("Er_Li") is not None:
                    old_rxn[d["mp_id"]] = (d["Er_Li"], d.get("rxn_Li", ""))
            except Exception:
                pass
    n_reuse_li = sum(1 for r in pool if r["mp_id"] in old_rxn)
    print(f"[反应能] vs-Li 旧checkpoint可复用 {n_reuse_li} 项")

    done_r = set()
    if os.path.exists(CKPT_RXN):
        for line in open(CKPT_RXN, encoding="utf-8"):
            try:
                done_r.add(json.loads(line)["mp_id"])
            except Exception:
                pass
    fr = open(CKPT_RXN, "a", encoding="utf-8")
    t0 = time.time(); n = 0
    for r in pool:
        if r["mp_id"] in done_r:
            continue
        n += 1
        res = {"mp_id": r["mp_id"], "method_vs_Li": "closed_PD_mix(norm=eV/atom) == GPIR@0V(mu_Li=mu0)",
               "method_vs_SSE": "closed_PD_mix(norm=eV/atom)"}
        try:
            cands = by_formula.get(r["formula"], [])
            if not cands:
                res["status"] = "no_entry_in_manifest"
            else:
                target = min(cands, key=lambda e: e.energy_per_atom)
                els = {el.symbol for el in target.composition.elements}
                # vs Li (封闭, 等价 0V 巨势): 优先复用旧结果
                if r["mp_id"] in old_rxn:
                    res["Er_Li"] = old_rxn[r["mp_id"]][0]
                    res["rxn_Li"] = old_rxn[r["mp_id"]][1]
                    res["src_Li"] = "reused_reaction_checkpoint_v2"
                else:
                    els_li = els | {"Li"}
                    rel = [e for e in raw if all(el.symbol in els_li for el in e.composition.elements)]
                    pd = PhaseDiagram(rel)
                    ir = InterfacialReactivity(target.composition, Composition("Li"), pd,
                                               norm=True, use_hull_energy=False)
                    kmin = min(ir.get_kinks(), key=lambda k: k[2])
                    res["Er_Li"] = round(float(kmin[2]), 4)
                    res["rxn_Li"] = str(kmin[3])
                    res["src_Li"] = "computed_now"
                # vs 真实 SSE
                for sname, sentry in SSE_REAL.items():
                    els2 = els | {el.symbol for el in sentry.composition.elements}
                    if len(els2) > 6:
                        res[f"Er_{sname}"] = None
                        res[f"note_{sname}"] = "chemspace>6"
                        continue
                    rel2 = [e for e in entries if all(el.symbol in els2 for el in e.composition.elements)]
                    pd2 = PhaseDiagram(rel2)
                    ir2 = InterfacialReactivity(target.composition, sentry.composition, pd2,
                                                norm=True, use_hull_energy=False)
                    kmin2 = min(ir2.get_kinks(), key=lambda k: k[2])
                    res[f"Er_{sname}"] = round(float(kmin2[2]), 4)
                    res[f"rxn_{sname}"] = str(kmin2[3])
                res["status"] = "ok"
        except Exception as ex:
            res["status"] = f"error:{type(ex).__name__}:{str(ex)[:80]}"
        fr.write(json.dumps(res, ensure_ascii=False) + "\n"); fr.flush()
        if n % 50 == 0:
            el = time.time() - t0
            print(f"[反应能] {n} done {el:.0f}s ({el/n:.2f}s/项) ETA {(len(pool)-len(done_r)-n)*el/n/60:.0f}min", flush=True)
    fr.close()
    print("[反应能] DONE", time.time() - t0)


if __name__ == "__main__":
    main()

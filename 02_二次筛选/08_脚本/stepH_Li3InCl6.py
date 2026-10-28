# -*- coding: utf-8 -*-
"""
Li3InCl6 路线全套筛选 (与 LZC/LYC 同方法)

步骤:
1. SSE 自身窗口: Li3InCl6 巨势窗口(已有 checkpoint, 直接读)
2. B 池: Er(B, Li3InCl6) 对母池 1860 项, 封闭两相混合, eV/atom
3. A|B 矩阵: A池340 × B池(入池) 组合反应能
4. 组合筛选: 同 step5 七规则, 输出通过/待核/动力学线索/未通过
5. 阈值敏感性 + 严格档 + Pareto

注意: Li3InCl6 的 Er_Li=-0.53, 它不能当 A 层/单层, 只能作被保护 SSE。
"""
import os, csv, json, pickle, time, warnings
warnings.filterwarnings("ignore")
import multiprocessing as mp

ROOT = os.path.dirname(os.path.abspath(__file__))
NEW = os.path.dirname(ROOT)
TOP = os.path.dirname(NEW)
WS = os.path.join(TOP, "AB池双层筛选")
OUT = os.path.join(NEW, "11_Li3InCl6路线")
os.makedirs(OUT, exist_ok=True)

NPROC = int(os.environ.get("NPROC", "12"))
SSE_NAME = "Li3InCl6"
SSE_MP_ID = "mp-676109"

from pymatgen.core import Composition, Element
from pymatgen.analysis.phase_diagram import PhaseDiagram, GrandPotentialPhaseDiagram, GrandPotPDEntry
from pymatgen.analysis.interface_reactions import InterfacialReactivity


def f(x):
    try:
        if x in ("", None, "inf"):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


# ---------- worker: Er(X, SSE) 与 Er(A,B) ----------
_W = {}

def _init(raw_pkl, sse_formula):
    raw = pickle.load(open(raw_pkl, "rb"))
    cands = [e for e in raw if e.composition.reduced_formula == sse_formula]
    sse = min(cands, key=lambda e: e.energy_per_atom)
    by_formula = {}
    for e in raw:
        by_formula.setdefault(e.composition.reduced_formula, []).append(e)
    li = min((e for e in raw if e.composition.reduced_formula == "Li"),
             key=lambda e: e.energy_per_atom)
    _W.update({"raw": raw, "sse": sse, "by_formula": by_formula, "li": li, "pd_cache": {}})


def _rxn_energy(f1, f2):
    """封闭两相混合反应能, eV/atom。f 为 reduced_formula。"""
    by_formula = _W["by_formula"]
    c1 = by_formula.get(f1)
    c2 = by_formula.get(f2)
    if not c1 or not c2:
        return None, "no_entry"
    e1 = min(c1, key=lambda x: x.energy_per_atom)
    e2 = min(c2, key=lambda x: x.energy_per_atom)
    els = {el.symbol for el in e1.composition.elements} | {el.symbol for el in e2.composition.elements}
    if len(els) > 6:
        return None, "chemspace>6"
    key = tuple(sorted(els))
    if key not in _W["pd_cache"]:
        rel = [e for e in _W["raw"] if all(el.symbol in els for el in e.composition.elements)]
        _W["pd_cache"][key] = PhaseDiagram(rel)
    pd = _W["pd_cache"][key]
    ir = InterfacialReactivity(e1.composition, e2.composition, pd, norm=True, use_hull_energy=False)
    kmin = min(ir.get_kinks(), key=lambda k: k[2])
    return round(float(kmin[2]), 4), "ok"


def _worker_b(job):
    mid, formula = job
    try:
        er, status = _rxn_energy(formula, SSE_NAME)
        return {"mp_id": mid, "formula": formula, f"Er_{SSE_NAME}": er, "status": status}
    except Exception as ex:
        return {"mp_id": mid, "formula": formula, f"Er_{SSE_NAME}": None,
                "status": f"error:{type(ex).__name__}"}


def _worker_ab(job):
    aid, a_formula, bid, b_formula = job
    try:
        er, status = _rxn_energy(a_formula, b_formula)
        return {"A": aid, "B": bid, "Er_AB": er, "method": "closed_PD_mix", "note": status}
    except Exception as ex:
        return {"A": aid, "B": bid, "Er_AB": None, "method": "closed_PD_mix",
                "note": f"error:{type(ex).__name__}"}


def main():
    raw_pkl = os.path.join(TOP, "multi_screen_Jeon2025", "02_数据获取", "pd_cache", "mp_gga_raw_entries.pkl")

    # ---- 0. 基础数据 ----
    merged = {r["mp_id"]: r for r in csv.DictReader(
        open(os.path.join(WS, "03_AB池", "母池_合并热力学数据.csv"), encoding="utf-8-sig"))}
    a_pool = [r["mp_id"] for r in csv.DictReader(
        open(os.path.join(WS, "03_AB池", "A池_靠Li.csv"), encoding="utf-8-sig")) if r["A池标记"] == "入池"]

    # SSE 窗口(读已有 checkpoint)
    sse_vred, sse_vox = None, None
    for line in open(os.path.join(WS, "02_热力学数据", "窗口_checkpoint.jsonl"), encoding="utf-8"):
        d = json.loads(line)
        if d["mp_id"] == SSE_MP_ID:
            sse_vred, sse_vox = d["V_red"], d["V_ox"]
    print(f"[{SSE_NAME}] 窗口 V_red={sse_vred} V_ox={sse_vox}", flush=True)

    # ---- 1. B 池: Er(B, Li3InCl6) ----
    ckpt_b = os.path.join(OUT, f"B池反应能_{SSE_NAME}.jsonl")
    done_b = set()
    if os.path.exists(ckpt_b):
        for line in open(ckpt_b, encoding="utf-8"):
            try:
                done_b.add(json.loads(line)["mp_id"])
            except Exception:
                pass
    todo_b = [(mid, merged[mid]["formula"]) for mid in merged if mid not in done_b]
    print(f"[{SSE_NAME}] B池反应能: 已完成 {len(done_b)}, 待算 {len(todo_b)}", flush=True)
    t0 = time.time()
    with open(ckpt_b, "a", encoding="utf-8") as fo:
        with mp.Pool(NPROC, initializer=_init, initargs=(raw_pkl, SSE_NAME)) as pool:
            for i, rec in enumerate(pool.imap_unordered(_worker_b, todo_b, chunksize=20)):
                fo.write(json.dumps(rec, ensure_ascii=False) + "\n")
                if (i + 1) % 200 == 0:
                    fo.flush()
                    el = time.time() - t0
                    print(f"  [B池] {i+1}/{len(todo_b)} {el:.0f}s ETA {(len(todo_b)-i-1)*el/(i+1)/60:.0f}min", flush=True)

    er_sse = {}
    for line in open(ckpt_b, encoding="utf-8"):
        d = json.loads(line)
        er_sse[d["mp_id"]] = d.get(f"Er_{SSE_NAME}")

    # B 池入池: Er >= -0.10
    b_pool = [mid for mid in merged if er_sse.get(mid) is not None and er_sse[mid] >= -0.10]
    print(f"[{SSE_NAME}] B池入池 {len(b_pool)} / {len(merged)}", flush=True)

    # ---- 2. A|B 矩阵(只对可能进通过/待核/动力学线索的组合) ----
    # 预筛: 排除 导体/导体 和 E3D>1.0 的组合, 这些必进未通过组, 不需要 Er_AB
    bv = {}
    for r in csv.DictReader(open(os.path.join(TOP, "softBV_Jeon2025", "06_结果汇总",
                                              "Jeon2025_1610项_Li迁移势垒汇总.csv"), encoding="utf-8-sig")):
        bv[r["mp_id"]] = r

    def _pre_pass(aid, bid):
        A, B = merged[aid], merged[bid]
        eg_a, eg_b = f(A["band_gap_eV"]), f(B["band_gap_eV"])
        if (eg_a is None or eg_a < EG_CUT) and (eg_b is None or eg_b < EG_CUT):
            return False
        for mid, inter_key in [(aid, "is_intermetallic"), (bid, "is_intermetallic")]:
            bvr = bv.get(mid)
            e3d = f(bvr["E_3D"]) if bvr else None
            inter = merged[mid].get(inter_key, "0") in ("1", "True", 1)
            if not inter and e3d is not None and e3d > E3D_MAX:
                return False
        return True

    EG_CUT, ER_TOL, E3D_MAX = 1.0, -0.10, 1.0
    DV_MARGIN, DV_EDGE = 0.10, 0.05

    ckpt_ab = os.path.join(OUT, f"AB反应能_{SSE_NAME}.jsonl")
    done_ab = set()
    if os.path.exists(ckpt_ab):
        for line in open(ckpt_ab, encoding="utf-8"):
            try:
                d = json.loads(line)
                done_ab.add((d["A"], d["B"]))
            except Exception:
                pass
    todo_ab = []
    for aid in a_pool:
        for bid in b_pool:
            if aid != bid and (aid, bid) not in done_ab and _pre_pass(aid, bid):
                todo_ab.append((aid, merged[aid]["formula"], bid, merged[bid]["formula"]))
    print(f"[{SSE_NAME}] AB矩阵(预筛后): 已完成 {len(done_ab)}, 待算 {len(todo_ab)}", flush=True)
    t0 = time.time()
    with open(ckpt_ab, "a", encoding="utf-8") as fo:
        with mp.Pool(NPROC, initializer=_init, initargs=(raw_pkl, SSE_NAME)) as pool:
            for i, rec in enumerate(pool.imap_unordered(_worker_ab, todo_ab, chunksize=20)):
                fo.write(json.dumps(rec, ensure_ascii=False) + "\n")
                if (i + 1) % 500 == 0:
                    fo.flush()
                    el = time.time() - t0
                    print(f"  [AB] {i+1}/{len(todo_ab)} {el:.0f}s ETA {(len(todo_ab)-i-1)*el/(i+1)/60:.0f}min", flush=True)

    er_ab = {}
    for line in open(ckpt_ab, encoding="utf-8"):
        d = json.loads(line)
        er_ab[(d["A"], d["B"])] = d.get("Er_AB")

    # ---- 3. 组合筛选(同 step5 规则) ----
    out = []
    for aid in a_pool:
        A = merged[aid]
        eg_a, vred_a, vox_a, er_li_a = f(A["band_gap_eV"]), f(A["V_red"]), f(A["V_ox"]), f(A["Er_Li"])
        bva = bv.get(aid)
        e3d_a = f(bva["E_3D"]) if bva else None
        a_inter = A.get("is_intermetallic", "0") in ("1", "True", 1)
        for bid in b_pool:
            if aid == bid:
                continue
            B = merged[bid]
            eg_b, vred_b, vox_b = f(B["band_gap_eV"]), f(B["V_red"]), f(B["V_ox"])
            er_b_sse = er_sse.get(bid)
            bvb = bv.get(bid)
            e3d_b = f(bvb["E_3D"]) if bvb else None
            b_inter = B.get("is_intermetallic", "0") in ("1", "True", 1)
            er_ab_v = er_ab.get((aid, bid))

            reasons_fail, reasons_pend = [], []
            block = (eg_a is not None and eg_a >= EG_CUT) or (eg_b is not None and eg_b >= EG_CUT)
            if not block:
                reasons_fail.append("规则3:导体/导体")
            if er_ab_v is None:
                reasons_pend.append("规则4:Er(A,B)缺失")
            elif min(er_li_a if er_li_a is not None else -9, er_ab_v,
                     er_b_sse if er_b_sse is not None else -9) < ER_TOL:
                reasons_fail.append("规则4:三界面Er<-0.10")
            if er_li_a is None or er_b_sse is None:
                reasons_pend.append("规则4:Li|A或B|SSE反应能缺失")

            kinetic_lead = False
            if None in (vred_a, vox_a, vred_b, vox_b):
                reasons_pend.append("规则5:窗口缺失")
            else:
                dv_ab = min(vox_a, vox_b) - max(vred_a, vred_b)
                dv_bs = min(vox_b, sse_vox) - max(vred_b, sse_vred)
                if dv_ab <= 0 or dv_bs <= 0:
                    kinetic_lead = True
                else:
                    if dv_ab <= DV_EDGE or dv_bs <= DV_EDGE:
                        reasons_pend.append("规则5:边界相接")
                    elif dv_ab < DV_MARGIN or dv_bs < DV_MARGIN:
                        reasons_pend.append("规则5:交叠<0.1V")
                    v1_lo = max(vred_a, vred_b, 0.0)
                    v2_hi = min(vox_b, sse_vox)
                    if v2_hi < v1_lo:
                        reasons_pend.append("规则5:层序电位连接不满足")

            for tag, inter, e3d in [("A", a_inter, e3d_a), ("B", b_inter, e3d_b)]:
                if inter:
                    reasons_pend.append(f"规则6:{tag}合金")
                elif e3d is None:
                    reasons_pend.append(f"规则6:{tag}层E3D缺失")
                elif e3d > E3D_MAX:
                    reasons_fail.append(f"规则6:{tag}层E3D>{E3D_MAX}")

            if reasons_fail:
                grp = "未通过组"
            elif kinetic_lead:
                grp = "动力学线索"
            elif reasons_pend:
                grp = "待核组"
            else:
                grp = "通过组"
            e3d_max = max([x for x in (e3d_a, e3d_b) if x is not None], default=None)
            dv_ab_s = round(min(vox_a, vox_b) - max(vred_a, vred_b), 2) if None not in (vred_a, vox_a, vred_b, vox_b) else ""
            dv_bs_s = round(min(vox_b, sse_vox) - max(vred_b, sse_vred), 2) if None not in (vred_b, vox_b) else ""
            out.append({"SSE": SSE_NAME, "A": aid, "A_formula": A["formula"], "B": bid, "B_formula": B["formula"],
                        "Eg_A": eg_a, "Eg_B": eg_b, "Er_Li_A": er_li_a, "Er_A_B": er_ab_v, "Er_B_SSE": er_b_sse,
                        "dV_AB": dv_ab_s, "dV_B_SSE": dv_bs_s, "E3D_A": e3d_a, "E3D_B": e3d_b,
                        "E3D_max": e3d_max, "group": grp,
                        "未通过条件": ";".join(reasons_fail), "待核原因": ";".join(reasons_pend)})

    from collections import Counter
    cnt = Counter(r["group"] for r in out)
    print(f"[{SSE_NAME}] 组合筛选: 总数 {len(out)} 分组 {dict(cnt)}", flush=True)
    with open(os.path.join(OUT, f"组合筛选_{SSE_NAME}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(out[0].keys()))
        w.writeheader(); w.writerows(out)

    # ---- 4. 严格档 ----
    strict = [r for r in out if r["group"] == "通过组" and
              min(r["Er_Li_A"], r["Er_A_B"], r["Er_B_SSE"]) >= -0.03 and
              r["E3D_max"] is not None and r["E3D_max"] <= 0.6]
    strict.sort(key=lambda r: (r["E3D_max"], -min(r["Er_Li_A"], r["Er_A_B"], r["Er_B_SSE"])))
    with open(os.path.join(OUT, f"strict_candidates_{SSE_NAME}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(out[0].keys()))
        w.writeheader(); w.writerows(strict)
    print(f"[{SSE_NAME}] 严格档 {len(strict)}", flush=True)


if __name__ == "__main__":
    main()

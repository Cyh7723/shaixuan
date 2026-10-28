# -*- coding: utf-8 -*-
"""
工作包C+D: 三界面反应产物 / 电子导通风险 / 反应体积变化

范围: 通过组 + 待核组中出现的唯一界面(按化学对去重), 控制规模。
  Li|A: 全部 A 池(340)
  B|SSE: 全部 B 池(LZC 793 / LYC 1070)
  A|B: 仅 通过组+待核组+动力学线索 中出现的唯一 (A,B) 对

对每个界面用 InterfacialReactivity 取全部 kinks(精确反应, 非 11 点采样),
保存: 配平反应、反应能、产物及系数、混合比例。
产物性质关联母池/MP 数据(带隙/金属性/体积); 未匹配到的标 unknown。

体积变化: δV = (V_products − V_reactants)/V_reactants, 用每化学式单元体积。
分别输出: 最负反应能 kink 的 δV; 全部 kink 中 max|δV|。
开放 Li 储库(仅 Li|A 界面)同时给 "计入金属Li体积" 与 "不计入(开放储库)" 两列。

输出:
  03_反应产物/interface_reactions.jsonl
  03_反应产物/reaction_product_properties.csv
  04_体积变化/reaction_volume_changes.csv
"""
import os, csv, json, pickle, time, warnings
warnings.filterwarnings("ignore")

# 多进程: 界面间相互独立, Windows 下用 spawn; 进程数由 NPROC 控制
import multiprocessing as mp
NPROC = int(os.environ.get("NPROC", "12"))

ROOT = os.path.dirname(os.path.abspath(__file__))
NEW = os.path.dirname(ROOT)
TOP = os.path.dirname(NEW)
WS = os.path.join(TOP, "AB池双层筛选")
OUT_R = os.path.join(NEW, "03_反应产物")
OUT_V = os.path.join(NEW, "04_体积变化")
os.makedirs(OUT_R, exist_ok=True); os.makedirs(OUT_V, exist_ok=True)

from pymatgen.core import Composition, Element
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.interface_reactions import InterfacialReactivity

RAW = None
SSE_REAL = None
BY_FORMULA = None
POOL_INFO = None  # mp_id -> dict(band_gap, volume...)
POOL_BY_FORMULA = None
_VOL_CACHE = {}


def f(x):
    try:
        if x in ("", None, "inf"):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def init():
    global RAW, SSE_REAL, BY_FORMULA, POOL_INFO, POOL_BY_FORMULA
    RAW = pickle.load(open(os.path.join(TOP, "multi_screen_Jeon2025", "02_数据获取",
                                        "pd_cache", "mp_gga_raw_entries.pkl"), "rb"))
    hal = pickle.load(open(os.path.join(TOP, "halide_vasp", "halide_entries_corrected.pkl"), "rb"))
    SSE_REAL = {"Li2ZrCl6": hal["Li2ZrCl6_beta"], "Li3YCl6": hal["Li3YCl6"]}
    BY_FORMULA = {}
    for e in RAW + list(SSE_REAL.values()):
        BY_FORMULA.setdefault(e.composition.reduced_formula, []).append(e)
    POOL_INFO = {}
    POOL_BY_FORMULA = {}
    for r in csv.DictReader(open(os.path.join(WS, "03_AB池", "母池_合并热力学数据.csv"),
                                 encoding="utf-8-sig")):
        POOL_INFO[r["mp_id"]] = r
        POOL_BY_FORMULA.setdefault(r["formula"], r)
    global _POOL_FORMULA
    _POOL_FORMULA = {mid: r["formula"] for mid, r in POOL_INFO.items()}


def props_by_formula(formula):
    """按化学式在母池清单查带隙/体积; 无匹配返回 None。"""
    if not hasattr(props_by_formula, "_mid_cache"):
        props_by_formula._mid_cache = {}
        for mid, rr in POOL_INFO.items():
            props_by_formula._mid_cache.setdefault(rr["formula"], mid)
    r = POOL_BY_FORMULA.get(formula)
    if r:
        return {"mp_id": props_by_formula._mid_cache.get(formula),
                "band_gap": f(r["band_gap_eV"]), "volume": f(r["volume"]),
                "nsites": f(r["nsites"]), "src": "母池"}
    return None


def vol_per_fu(formula):
    """每化学式单元体积: 母池有则用之; 否则用相图 entries 中该化学式最低能项的结构体积。"""
    if formula in _VOL_CACHE:
        return _VOL_CACHE[formula]
    out = None
    p = props_by_formula(formula)
    if p and p["volume"] and p["nsites"]:
        comp = Composition(formula)
        n_fu_cell = p["nsites"] / comp.num_atoms
        out = p["volume"] / n_fu_cell
    else:
        cands = BY_FORMULA.get(formula)
        if cands:
            e = min(cands, key=lambda x: x.energy_per_atom)
            try:
                s = e.structure
                n_fu_cell = s.num_sites / e.composition.num_atoms
                out = s.volume / n_fu_cell
            except Exception:
                out = None
    _VOL_CACHE[formula] = out
    return out


def gap_of(formula):
    p = props_by_formula(formula)
    if p:
        return p["band_gap"], p["mp_id"]
    return None, None


def analyze_interface(c1, c2, pd_entries, label, open_li=False):
    """InterfacialReactivity 全 kink 分析。
    pymatgen>=2024: get_kinks() 返回 (idx, x, Er, rxn, ...), rxn.products/reactants
    为 list[Composition] (无系数); 系数需从 rxn 字符串解析或用 rxn.get_coeff。"""
    ir = InterfacialReactivity(c1, c2, pd_entries, norm=True, use_hull_energy=False)
    out_kinks = []
    for k in ir.get_kinks():
        x, erxn, rxn = k[1], k[2], k[3]
        # 系数: pymatgen Reaction.get_coeff(comp) 对产物返回负值(反应约定),
        # 这里统一取绝对值; 旧版 products 为 dict。
        def _pairs(x_obj):
            if hasattr(x_obj, "items"):
                return list(x_obj.items())
            return [(c, abs(rxn.get_coeff(c)) if hasattr(rxn, "get_coeff") else 1.0) for c in x_obj]
        prods = {}
        for comp, coef in _pairs(rxn.products):
            prods[comp.reduced_formula] = float(coef)
        reacts = {}
        for comp, coef in _pairs(rxn.reactants):
            reacts[comp.reduced_formula] = float(coef)
        v_p = 0.0
        v_ok = True
        for comp, coef in _pairs(rxn.products):
            v = vol_per_fu(comp.reduced_formula)
            if v is None:
                v_ok = False
                continue
            v_p += float(coef) * v
        v_re = 0.0
        for comp, coef in _pairs(rxn.reactants):
            v = vol_per_fu(comp.reduced_formula)
            if v is None:
                v_ok = False
                continue
            v_re += float(coef) * v
        dV = (v_p - v_re) / v_re if (v_ok and v_re > 0) else None
        dV_noLi = None
        if open_li and v_ok:
            v_p2 = sum(float(coef) * (vol_per_fu(comp.reduced_formula) or 0.0)
                       for comp, coef in _pairs(rxn.products)
                       if comp.reduced_formula != "Li")
            dV_noLi = (v_p2 - v_re) / v_re if v_re > 0 else None
        out_kinks.append({
            "x": round(float(x), 4), "Er_eV_per_atom": round(float(erxn), 4),
            "reaction": str(rxn), "products": prods, "reactants": reacts,
            "dV": round(dV, 4) if dV is not None else None,
            "dV_openLi_excl_Li": round(dV_noLi, 4) if dV_noLi is not None else None,
            "vol_complete": v_ok,
        })
    return out_kinks


# ---- 多进程 worker ----
_W = {}


def _worker_init(raw_pkl, hal_pkl, pool_formula):
    import pickle
    from pymatgen.analysis.phase_diagram import PhaseDiagram
    from pymatgen.core import Composition
    global _POOL_FORMULA, POOL_BY_FORMULA, POOL_INFO, BY_FORMULA
    _POOL_FORMULA = pool_formula
    raw = pickle.load(open(raw_pkl, "rb"))
    hal = pickle.load(open(hal_pkl, "rb"))
    sse = {"Li2ZrCl6": hal["Li2ZrCl6_beta"], "Li3YCl6": hal["Li3YCl6"]}
    all_e = raw + list(sse.values())
    by_formula = {}
    for e in all_e:
        by_formula.setdefault(e.composition.reduced_formula, []).append(e)
    BY_FORMULA = by_formula
    li = min((e for e in raw if e.composition.reduced_formula == "Li"),
             key=lambda e: e.energy_per_atom)
    _W.update({"raw": raw, "sse": sse, "all_e": all_e, "by_formula": by_formula,
               "li": li, "pd_cache": {}})
    # analyze_interface -> vol_per_fu 依赖 POOL_BY_FORMULA; worker 里自构(只含 formula->volume/nsites)
    POOL_INFO = {}
    POOL_BY_FORMULA = {}
    for mid, fm in pool_formula.items():
        POOL_BY_FORMULA.setdefault(fm, {"formula": fm, "band_gap_eV": "", "volume": "", "nsites": ""})


def _worker_one(job):
    lab, s1, s2 = job
    try:
        raw, sse, all_e = _W["raw"], _W["sse"], _W["all_e"]
        by_formula, li = _W["by_formula"], _W["li"]
        if s1 == "Li_metal":
            e1 = li
        else:
            f1 = _POOL_FORMULA[s1]
            cands = by_formula.get(f1, [])
            if not cands:
                raise ValueError("no_entry")
            e1 = min(cands, key=lambda x: x.energy_per_atom)
        if s2 in sse:
            e2 = sse[s2]
        elif s2 == "Li_metal":
            e2 = li
        else:
            f2 = _POOL_FORMULA[s2]
            cands2 = by_formula.get(f2, [])
            if not cands2:
                raise ValueError("no_entry")
            e2 = min(cands2, key=lambda x: x.energy_per_atom)
        c1, c2 = e1.composition, e2.composition
        els = {el.symbol for el in c1.elements} | {el.symbol for el in c2.elements}
        if len(els) > 6:
            raise ValueError("chemspace>6")
        key = tuple(sorted(els))
        if key not in _W["pd_cache"]:
            rel = [e for e in all_e if all(el.symbol in els for el in e.composition.elements)]
            _W["pd_cache"][key] = PhaseDiagram(rel)
        pd = _W["pd_cache"][key]
        open_li = (s1 == "Li_metal")
        kinks = analyze_interface(c1, c2, pd, lab, open_li=open_li)
        worst_er = min(kinks, key=lambda k: k["Er_eV_per_atom"])
        return {
            "interface": lab, "side1": s1, "side2": s2,
            "side1_formula": c1.reduced_formula, "side2_formula": c2.reduced_formula,
            "n_kinks": len(kinks),
            "Er_min": worst_er["Er_eV_per_atom"],
            "Er_min_reaction": worst_er["reaction"],
            "Er_min_products": worst_er["products"],
            "Er_min_dV": worst_er["dV"],
            "Er_min_dV_openLi_excl_Li": worst_er.get("dV_openLi_excl_Li"),
            "max_abs_dV": round(max(abs(k["dV"]) for k in kinks if k["dV"] is not None), 4)
                            if any(k["dV"] is not None for k in kinks) else None,
            "kinks": kinks,
        }
    except Exception as ex:
        return {"interface": lab, "side1": s1, "side2": s2,
                "error": f"{type(ex).__name__}:{str(ex)[:100]}"}


_POOL_FORMULA = {}


def main():
    init()
    # ---- 收集需计算的界面 ----
    a_ids = [r["mp_id"] for r in csv.DictReader(open(os.path.join(WS, "03_AB池", "A池_靠Li.csv"),
                                                     encoding="utf-8-sig")) if r["A池标记"] == "入池"]
    interfaces = []
    for aid in a_ids:
        interfaces.append(("Li|A", "Li_metal", aid))
    for sse in ("Li2ZrCl6", "Li3YCl6"):
        b_ids = [r["mp_id"] for r in csv.DictReader(
            open(os.path.join(WS, "03_AB池", f"B池_靠SSE_{sse}.csv"), encoding="utf-8-sig"))
            if r["B池标记"] == "入池"]
        for bid in b_ids:
            interfaces.append(("B|SSE", bid, sse))
    for sse in ("Li2ZrCl6", "Li3YCl6"):
        seen = set()
        with open(os.path.join(WS, "05_结果", f"组合筛选_{sse}.csv"), encoding="utf-8-sig") as fo:
            for r in csv.DictReader(fo):
                if r["group"] in ("通过组", "待核组", "动力学线索"):
                    key = (r["A"], r["B"])
                    if key not in seen:
                        seen.add(key)
                        interfaces.append(("A|B", r["A"], r["B"]))
    uniq = {}
    for lab, s1, s2 in interfaces:
        uniq[(lab, s1, s2)] = True
    interfaces = sorted(uniq.keys())
    print(f"[C+D] 唯一界面数 {len(interfaces)} (Li|A {sum(1 for i in interfaces if i[0]=='Li|A')}, "
          f"A|B {sum(1 for i in interfaces if i[0]=='A|B')}, B|SSE {sum(1 for i in interfaces if i[0]=='B|SSE')})",
          flush=True)

    # checkpoint 续跑: 已完成的界面跳过
    CKPT = os.path.join(OUT_R, "interface_reactions.jsonl")
    done = set()
    if os.path.exists(CKPT):
        for line in open(CKPT, encoding="utf-8"):
            try:
                d = json.loads(line)
                done.add((d.get("interface"), d.get("side1"), d.get("side2")))
            except Exception:
                pass
    todo = [i for i in interfaces if i not in done]
    print(f"[C+D] 已完成 {len(done)}, 待算 {len(todo)}", flush=True)

    # 多进程: 每个 worker 预加载 entries 并自建相图缓存
    li_entry = min((e for e in RAW if e.composition.reduced_formula == "Li"),
                   key=lambda e: e.energy_per_atom)
    all_entries = RAW + list(SSE_REAL.values())

    fout = open(CKPT, "a", encoding="utf-8")
    t0 = time.time()
    n_done = n_err = 0
    prod_props = {}

    with mp.Pool(NPROC, initializer=_worker_init,
                 initargs=(os.path.join(TOP, "multi_screen_Jeon2025", "02_数据获取",
                                        "pd_cache", "mp_gga_raw_entries.pkl"),
                           os.path.join(TOP, "halide_vasp", "halide_entries_corrected.pkl"),
                           _POOL_FORMULA)) as pool:
        for rec in pool.imap_unordered(_worker_one, todo, chunksize=20):
            n_done += 1
            if "error" in rec:
                n_err += 1
            fout.write(json.dumps(rec, ensure_ascii=False) + "\n")
            if n_done % 200 == 0:
                fout.flush()
                el = time.time() - t0
                print(f"[C+D] {n_done}/{len(todo)} err={n_err} {el:.0f}s "
                      f"ETA {(len(todo)-n_done)*el/max(1,n_done)/60:.0f}min", flush=True)
    fout.close()
    print(f"[C+D] 界面计算完成 {n_done} err {n_err}", flush=True)

    # ---- 汇总: 从 jsonl 重读, 生成 vol csv 与产物性质表 ----
    vol_rows = []
    for line in open(CKPT, encoding="utf-8"):
        try:
            d = json.loads(line)
        except Exception:
            continue
        if "error" in d:
            continue
        all_prods = set()
        for k in d.get("kinks", []):
            all_prods.update(k.get("products", {}).keys())
        for pf in all_prods:
            if pf not in prod_props:
                g, mid = gap_of(pf)
                prod_props[pf] = {"band_gap": g, "mp_id": mid,
                                  "vol_per_fu": vol_per_fu(pf),
                                  "is_metal": (g is not None and g < 0.05),
                                  "src": "母池" if mid else ("entries" if vol_per_fu(pf) else "unknown")}
        vol_rows.append({
            "interface": d["interface"], "side1": d["side1"], "side2": d["side2"],
            "side1_formula": d.get("side1_formula"), "side2_formula": d.get("side2_formula"),
            "Er_min": d.get("Er_min"),
            "dV_at_Er_min": d.get("Er_min_dV"),
            "dV_at_Er_min_openLi_excl_Li": d.get("Er_min_dV_openLi_excl_Li"),
            "max_abs_dV_all_kinks": d.get("max_abs_dV"),
            "vol_complete_at_Er_min": d.get("kinks", [{}])[0].get("vol_complete") if d.get("kinks") else None,
            "has_gas_product": False,
            "conducting_products": ";".join(d.get("conducting_products", [])),
            "min_product_gap": d.get("min_product_gap"),
            "product_gap_coverage": d.get("product_gap_coverage"),
        })

    with open(os.path.join(OUT_V, "reaction_volume_changes.csv"), "w", newline="",
              encoding="utf-8-sig") as fo:
        fields = ["interface", "side1", "side2", "side1_formula", "side2_formula", "Er_min",
                  "dV_at_Er_min", "dV_at_Er_min_openLi_excl_Li", "max_abs_dV_all_kinks",
                  "vol_complete_at_Er_min", "has_gas_product", "conducting_products",
                  "min_product_gap", "product_gap_coverage"]
        w = csv.DictWriter(fo, fieldnames=fields)
        w.writeheader(); w.writerows(vol_rows)

    with open(os.path.join(OUT_R, "reaction_product_properties.csv"), "w", newline="",
              encoding="utf-8-sig") as fo:
        fields = ["product_formula", "mp_id", "band_gap", "vol_per_fu", "is_metal", "src"]
        w = csv.DictWriter(fo, fieldnames=fields)
        w.writeheader()
        for pf, m in sorted(prod_props.items()):
            w.writerow({"product_formula": pf, "mp_id": m["mp_id"], "band_gap": m["band_gap"],
                        "vol_per_fu": None if m["vol_per_fu"] is None else round(m["vol_per_fu"], 3),
                        "is_metal": int(m["is_metal"]), "src": m["src"]})
    print(f"[C+D] DONE 产物种类 {len(prod_props)}")


if __name__ == "__main__":
    main()

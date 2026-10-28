# -*- coding: utf-8 -*-
"""
工作包F: 通道几何(BVSE 维度) + 传输证据分级 + ASR 情景 + 电位敏感性

范围: 严格档(E_worst>=-0.03, E_max<=0.6)涉及的唯一材料(LZC+LYC 并集)。

1. 通道几何: 用 BVSE E_1D/E_2D/E_3D 推断贯通维度
   (E_1D<E_2D≈E_3D → 1D 通道; E_2D<E_3D → 2D; E_1D≈E_2D≈E_3D → 3D);
   不做 CAVD 几何瓶颈(需单独装 cavd 包, 列为后续)。
2. 传输证据分级 A/B/C/D/U; ASR_bulk = tA/σA + tB/σB (仅当有可信 σ)。
3. 电位敏感性: 对严格档材料, 在 V ∈ [-0.5, 7] 扫描巨势分解,
   输出分解产物、是否出现导电产物、分类是否变化。

输出:
  06_通道几何/channel_geometry.csv
  06_通道几何/transport_evidence.csv
  06_通道几何/potential_sensitivity.csv
"""
import os, csv, json, pickle, warnings
warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.abspath(__file__))
NEW = os.path.dirname(ROOT)
TOP = os.path.dirname(NEW)
WS = os.path.join(TOP, "AB池双层筛选")
OUT = os.path.join(NEW, "06_通道几何")
os.makedirs(OUT, exist_ok=True)


def f(x):
    try:
        if x in ("", None, "inf"):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def load_strict_ids():
    ids = set()
    for tag in ("LZC", "LYC"):
        p = os.path.join(NEW, "01_阈值敏感性", f"strict_candidates_{tag}.csv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p, encoding="utf-8-sig")):
                ids.add(r["A"]); ids.add(r["B"])
    return ids


def channel_dim(e1, e2, e3, tol=0.05):
    if None in (e1, e2, e3):
        return "未知"
    if abs(e1 - e3) <= tol and abs(e2 - e3) <= tol:
        return "3D"
    if e2 < e3 - tol:
        return "2D(层内)"
    if e1 < e2 - tol:
        return "1D(链状)"
    return "各向异性不明显"


def transport_grade(mp_id, bv_row):
    """A: 实验σ; B: 缺陷/无序模拟; C: NEB; D: BVSE; U: 未知"""
    # 本数据库只有 BVSE/softBV; NEB/实验未系统收录
    if bv_row is None:
        return "U"
    src = bv_row.get("barrier_source", "")
    if "BVEL13k" in src:
        return "D"
    if "softBV" in src:
        return "D"
    return "U"


def main():
    ids = load_strict_ids()
    print(f"[F] 严格档唯一材料 {len(ids)}")

    bv = {}
    for r in csv.DictReader(open(os.path.join(TOP, "softBV_Jeon2025", "06_结果汇总",
                                              "Jeon2025_1610项_Li迁移势垒汇总.csv"), encoding="utf-8-sig")):
        bv[r["mp_id"]] = r
    merged = {r["mp_id"]: r for r in csv.DictReader(
        open(os.path.join(WS, "03_AB池", "母池_合并热力学数据.csv"), encoding="utf-8-sig"))}

    # ---- 通道几何 + 传输证据 ----
    ch_rows, tr_rows = [], []
    for mid in sorted(ids):
        b = bv.get(mid)
        e1, e2, e3 = (f(b["E_1D"]), f(b["E_2D"]), f(b["E_3D"])) if b else (None, None, None)
        ch_rows.append({
            "mp_id": mid, "formula": merged[mid]["formula"],
            "E_1D": e1, "E_2D": e2, "E_3D": e3,
            "贯通维度推断": channel_dim(e1, e2, e3),
            "取向依赖注意": "1D/2D 通道对取向敏感" if channel_dim(e1, e2, e3).startswith(("1D", "2D")) else "",
            "CAVD瓶颈半径": "未算(需 cavd 包, 后续)", "barrier_source": b.get("barrier_source", "") if b else "",
        })
        grade = transport_grade(mid, b)
        # ASR 需要 σ; BVSE 不直接换算电导率(方案禁止)
        tr_rows.append({
            "mp_id": mid, "formula": merged[mid]["formula"],
            "证据等级": grade, "证据来源": b.get("barrier_source", "") if b else "",
            "E3D_eV": e3,
            "ASR_5nm": "无法生成(无可信σ)", "ASR_10nm": "无法生成(无可信σ)", "ASR_20nm": "无法生成(无可信σ)",
            "备注": "BVSE 势垒不换算室温电导率; 无界面电阻" if grade == "D" else "缺迁移数据",
        })
    with open(os.path.join(OUT, "channel_geometry.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(ch_rows[0].keys())); w.writeheader(); w.writerows(ch_rows)
    with open(os.path.join(OUT, "transport_evidence.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(tr_rows[0].keys())); w.writeheader(); w.writerows(tr_rows)

    # ---- 电位敏感性 (巨势分解 vs V) ----
    raw = pickle.load(open(os.path.join(TOP, "multi_screen_Jeon2025", "02_数据获取",
                                        "pd_cache", "mp_gga_raw_entries.pkl"), "rb"))
    hal = pickle.load(open(os.path.join(TOP, "halide_vasp", "halide_entries_corrected.pkl"), "rb"))
    sse = {"Li2ZrCl6": hal["Li2ZrCl6_beta"], "Li3YCl6": hal["Li3YCl6"]}
    entries = raw + list(sse.values())
    by_formula = {}
    for e in entries:
        by_formula.setdefault(e.composition.reduced_formula, []).append(e)
    mu0 = min((e for e in entries if e.composition.reduced_formula == "Li"),
              key=lambda e: e.energy_per_atom).energy_per_atom

    from pymatgen.core import Element
    from pymatgen.analysis.phase_diagram import PhaseDiagram, GrandPotentialPhaseDiagram, GrandPotPDEntry

    gap_map = {mid: f(merged[mid]["band_gap_eV"]) for mid in ids}
    # 产物带隙从 reaction_product_properties.csv 读取
    prod_gap = {}
    pp = os.path.join(NEW, "03_反应产物", "reaction_product_properties.csv")
    if os.path.exists(pp):
        for r in csv.DictReader(open(pp, encoding="utf-8-sig")):
            prod_gap[r["product_formula"]] = f(r["band_gap"])

    Vs = [round(-0.5 + 0.25 * i, 2) for i in range(31)]
    sens_rows = []
    for mid in sorted(ids):
        formula = merged[mid]["formula"]
        cands = by_formula.get(formula)
        if not cands:
            sens_rows.append({"mp_id": mid, "formula": formula, "status": "no_entry"})
            continue
        target = min(cands, key=lambda e: e.energy_per_atom)
        els = {el.symbol for el in target.composition.elements}
        rel = [e for e in raw if all(el.symbol in els | {"Li"} for el in e.composition.elements)]
        pd = PhaseDiagram(rel)
        _, own = pd.get_decomp_and_e_above_hull(target)
        tol = own + 0.0501
        states = []
        for V in Vs:
            mu = mu0 - V
            gpd = GrandPotentialPhaseDiagram(rel, {Element("Li"): mu})
            ge = GrandPotPDEntry(target, {Element("Li"): mu})
            _, ea = gpd.get_decomp_and_e_above_hull(ge, allow_negative=True)
            stable = ea <= tol
            prods = ""
            conducting = ""
            if not stable:
                decomp, _ = gpd.get_decomp_and_e_above_hull(ge, allow_negative=True)
                pfs = [k.composition.reduced_formula for k in decomp]
                # 巨势相图分解自动扣除开放元素 Li; 还原分解(V<V_red)会析出 Li 金属,
                # 需显式补回: 若相对原料失 Li 则产物含 Li。
                prods = "+".join(sorted(set(pfs)))
                # 导电产物: 金属 Li(还原侧) + 已知金属相
                cond = [p for p in pfs if prod_gap.get(p) is not None and prod_gap[p] < 0.05 and p != "Li"]
                if V < 0.5:  # 还原分解区, Li 析出
                    cond.append("Li(金属,还原析出)")
                conducting = ";".join(sorted(set(cond)))
            states.append({"V": V, "stable": stable, "products": prods, "conducting": conducting})
        # 汇总: 稳定区间, 导电产物出现的 V 范围
        stable_Vs = [s["V"] for s in states if s["stable"]]
        cond_Vs = [s["V"] for s in states if s["conducting"]]
        cond_set = set()
        for s in states:
            if s["conducting"]:
                cond_set.update(s["conducting"].split(";"))
        sens_rows.append({
            "mp_id": mid, "formula": formula, "status": "ok",
            "V_red_calc": min(stable_Vs) if stable_Vs else None,
            "V_ox_calc": max(stable_Vs) if stable_Vs else None,
            "truncated_low": (min(stable_Vs) == -0.5) if stable_Vs else None,
            "truncated_high": (max(stable_Vs) == 7.0) if stable_Vs else None,
            "conducting_product_V_range": f"{min(cond_Vs)}~{max(cond_Vs)}" if cond_Vs else "",
            "conducting_products": ";".join(sorted(cond_set)),
            "分类是否改变": "是" if cond_Vs else "否",
        })
    with open(os.path.join(OUT, "potential_sensitivity.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        fields = ["mp_id", "formula", "status", "V_red_calc", "V_ox_calc", "truncated_low",
                  "truncated_high", "conducting_product_V_range", "conducting_products", "分类是否改变"]
        w = csv.DictWriter(fo, fieldnames=fields); w.writeheader(); w.writerows(sens_rows)
    print(f"[F] 通道 {len(ch_rows)}, 传输 {len(tr_rows)}, 电位敏感性 {len(sens_rows)}")


if __name__ == "__main__":
    main()

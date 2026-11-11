# -*- coding: utf-8 -*-
"""
工作包A: 阈值敏感性 + 严格档输出 (不改原管线, 全部新算)

对 LZC/LYC 两路线统计:
  E_worst = min[Er(Li,A), Er(A,B), Er(B,SSE)] 下限 ∈ {-0.10, -0.05, -0.03} eV/atom
  E_max   = max(E3D_A, E3D_B)                  上限 ∈ {0.5, 0.6, 0.7, 1.0} eV
  Eg 阻隔门槛 ∈ {1.0, 2.0} eV (至少一层)
  窗口交叠下限敏感性 ΔV ∈ {0.1, 0.3, 0.5} V

只统计"数据齐全"的组合(Er 三值齐、E3D 两值齐、Eg 两值齐、窗口四值齐);
缺失单独计数, 不把缺失当通过。

输出(01_阈值敏感性/):
  threshold_sensitivity.csv      全阈值交叉统计
  attrition_by_stage.csv         逐级淘汰(基线门槛顺序)
  strict_candidates_LZC.csv      严格优先档(E_worst>=-0.03, E_max<=0.6, 原Eg/窗口规则)
  strict_candidates_LYC.csv
"""
import os, csv, json
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
NEW = os.path.dirname(ROOT)
TOP = os.path.dirname(NEW)
WS = os.path.join(TOP, "AB池双层筛选")
OUT = os.path.join(TOP, "99_探索_支线", "01_阈值敏感性")
os.makedirs(OUT, exist_ok=True)

ER_CUTS = [-0.10, -0.05, -0.03]
E3D_CUTS = [0.5, 0.6, 0.7, 1.0]
EG_CUTS = [1.0, 2.0]
DV_CUTS = [0.1, 0.3, 0.5]


def f(x):
    try:
        if x in ("", None, "inf"):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def load(sse):
    swin = json.load(open(os.path.join(WS, "02_热力学数据", "SSE窗口.json"), encoding="utf-8"))
    sse_vred, sse_vox = swin[sse]["V_red"], swin[sse]["V_ox"]
    merged = {r["mp_id"]: r for r in csv.DictReader(
        open(os.path.join(WS, "03_AB池", "母池_合并热力学数据.csv"), encoding="utf-8-sig"))}
    a_pool = [r["mp_id"] for r in csv.DictReader(
        open(os.path.join(WS, "03_AB池", "A池_靠Li.csv"), encoding="utf-8-sig")) if r["A池标记"] == "入池"]
    b_pool = [r["mp_id"] for r in csv.DictReader(
        open(os.path.join(WS, "03_AB池", f"B池_靠SSE_{sse}.csv"), encoding="utf-8-sig")) if r["B池标记"] == "入池"]
    er_ab = {}
    for line in open(os.path.join(WS, "04_组合筛选", f"AB反应能_checkpoint_{sse}.jsonl"), encoding="utf-8"):
        d = json.loads(line)
        er_ab[(d["A"], d["B"])] = d.get("Er_AB")
    bv = {}
    for r in csv.DictReader(open(os.path.join(TOP, "softBV_Jeon2025", "06_结果汇总",
                                              "Jeon2025_1610项_Li迁移势垒汇总.csv"), encoding="utf-8-sig")):
        bv[r["mp_id"]] = r
    return sse_vred, sse_vox, merged, a_pool, b_pool, er_ab, bv


def family(formula, merged_row):
    """粗材料家族标签"""
    ch = merged_row.get("chemsys", "")
    f_ = formula
    if merged_row.get("is_intermetallic") in ("1", 1, "True"):
        return "金属间化合物"
    for el, fam in [("F", "氟化物"), ("Cl", "氯化物"), ("Br", "溴化物"), ("I", "碘化物"),
                    ("O", "氧化物/含氧盐"), ("S", "硫化物"), ("Se", "硒化物"), ("Te", "碲化物"),
                    ("N", "氮化物"), ("P", "磷化物/磷酸盐"), ("B", "硼化物/硼酸盐"),
                    ("C", "含碳"), ("H", "含氢")]:
        if el in ch.split("-"):
            return fam
    return "其他"


def run(sse, tag):
    sse_vred, sse_vox, merged, a_pool, b_pool, er_ab, bv = load(sse)

    rows_combo = []   # 每个数据齐全组合的最小信息行
    n_missing = Counter()
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
            er_b_sse = f(B[f"Er_{sse}"])
            bvb = bv.get(bid)
            e3d_b = f(bvb["E_3D"]) if bvb else None
            b_inter = B.get("is_intermetallic", "0") in ("1", "True", 1)
            er_ab_v = er_ab.get((aid, bid))

            # 缺失统计(任一关键项缺失 → 不进阈值统计, 记缺失类型)
            miss = []
            if None in (er_li_a, er_b_sse) or er_ab_v is None:
                miss.append("Er")
            if e3d_a is None or e3d_b is None:
                miss.append("E3D" + ("+合金" if (a_inter or b_inter) else ""))
            if eg_a is None or eg_b is None:
                miss.append("Eg")
            if None in (vred_a, vox_a, vred_b, vox_b):
                miss.append("窗口")
            if a_inter or b_inter:
                miss.append("含合金层(BVSE不适用)")
            if miss:
                n_missing["+".join(sorted(set(miss)))] += 1
                continue

            e_worst = min(er_li_a, er_ab_v, er_b_sse)
            e_max = max(e3d_a, e3d_b)
            dv_ab = min(vox_a, vox_b) - max(vred_a, vred_b)
            dv_bs = min(vox_b, sse_vox) - max(vred_b, sse_vred)
            dv_min = min(dv_ab, dv_bs)
            rows_combo.append({
                "A": aid, "A_formula": A["formula"], "B": bid, "B_formula": B["formula"],
                "Eg_A": eg_a, "Eg_B": eg_b, "E_worst": e_worst, "E_max": e_max,
                "dV_AB": round(dv_ab, 3), "dV_B_SSE": round(dv_bs, 3), "dV_min": round(dv_min, 3),
                "Er_Li_A": er_li_a, "Er_A_B": er_ab_v, "Er_B_SSE": er_b_sse,
                "E3D_A": e3d_a, "E3D_B": e3d_b,
                "架构": ("双绝缘" if (eg_a >= 1.0 and eg_b >= 1.0) else
                        ("绝缘A/导电B" if eg_a >= 1.0 else "导电A/绝缘B")),
                "家族A": family(A["formula"], A), "家族B": family(B["formula"], B),
            })

    # ---- 阈值交叉统计 ----
    sens = []
    for erc in ER_CUTS:
        for e3c in E3D_CUTS:
            for egc in EG_CUTS:
                for dvc in DV_CUTS:
                    sel = [r for r in rows_combo
                           if r["E_worst"] >= erc and r["E_max"] <= e3c
                           and (r["Eg_A"] >= egc or r["Eg_B"] >= egc)
                           and r["dV_min"] >= dvc]
                    sens.append({
                        "SSE": sse, "Er下限": erc, "E3D上限": e3c, "Eg门槛": egc, "dV下限": dvc,
                        "配对数": len(sel),
                        "独立A数": len({r["A"] for r in sel}), "独立B数": len({r["B"] for r in sel}),
                        "独立化学式对数": len({(r["A_formula"], r["B_formula"]) for r in sel}),
                    })

    with open(os.path.join(OUT, f"threshold_sensitivity_{tag}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(sens[0].keys()))
        w.writeheader(); w.writerows(sens)

    # ---- 逐级淘汰(基线门槛顺序, 固定 Eg>=1.0, dV>=0.1) ----
    stages = []
    def count(sel):
        return (len(sel), len({r["A"] for r in sel}), len({r["B"] for r in sel}),
                len({(r["A_formula"], r["B_formula"]) for r in sel}))
    n0 = len(a_pool) * len(b_pool) - len(set(a_pool) & set(b_pool))
    stages.append({"SSE": sse, "阶段": "0_全部组合(排除A==B)", "配对数": n0, "独立A数": len(a_pool),
                   "独立B数": len(b_pool), "独立化学式对数": "", "相对上一步保留%": 100.0, "备注": ""})
    stages.append({"SSE": sse, "阶段": "1_关键数据齐全", "配对数": count(rows_combo)[0],
                   "独立A数": count(rows_combo)[1], "独立B数": count(rows_combo)[2],
                   "独立化学式对数": count(rows_combo)[3],
                   "相对上一步保留%": round(100 * len(rows_combo) / n0, 2),
                   "备注": f"缺失明细: {dict(n_missing)}"})
    s2 = [r for r in rows_combo if r["E_worst"] >= -0.10]
    stages.append({"SSE": sse, "阶段": "2_三界面Er>=-0.10", "配对数": count(s2)[0],
                   "独立A数": count(s2)[1], "独立B数": count(s2)[2], "独立化学式对数": count(s2)[3],
                   "相对上一步保留%": round(100 * len(s2) / max(1, len(rows_combo)), 2), "备注": ""})
    s3 = [r for r in s2 if r["Eg_A"] >= 1.0 or r["Eg_B"] >= 1.0]
    stages.append({"SSE": sse, "阶段": "3_电子阻隔(Eg>=1至少一层)", "配对数": count(s3)[0],
                   "独立A数": count(s3)[1], "独立B数": count(s3)[2], "独立化学式对数": count(s3)[3],
                   "相对上一步保留%": round(100 * len(s3) / max(1, len(s2)), 2),
                   "备注": "架构分布: " + str(dict(Counter(r["架构"] for r in s3)))})
    s4 = [r for r in s3 if r["dV_min"] >= 0.1]
    stages.append({"SSE": sse, "阶段": "4_窗口交叠>=0.1V", "配对数": count(s4)[0],
                   "独立A数": count(s4)[1], "独立B数": count(s4)[2], "独立化学式对数": count(s4)[3],
                   "相对上一步保留%": round(100 * len(s4) / max(1, len(s3)), 2), "备注": ""})
    s5 = [r for r in s4 if r["E_max"] <= 1.0]
    stages.append({"SSE": sse, "阶段": "5_BVSE E3D<=1.0", "配对数": count(s5)[0],
                   "独立A数": count(s5)[1], "独立B数": count(s5)[2], "独立化学式对数": count(s5)[3],
                   "相对上一步保留%": round(100 * len(s5) / max(1, len(s4)), 2),
                   "备注": "对应原'通过组'(条件顺序不同, 数量应≈报告通过组)"})
    s6 = [r for r in s5 if r["E_worst"] >= -0.03 and r["E_max"] <= 0.6]
    stages.append({"SSE": sse, "阶段": "6_严格档(Er>=-0.03 & E3D<=0.6)", "配对数": count(s6)[0],
                   "独立A数": count(s6)[1], "独立B数": count(s6)[2], "独立化学式对数": count(s6)[3],
                   "相对上一步保留%": round(100 * len(s6) / max(1, len(s5)), 2), "备注": "工程分档, 非普适阈值"})

    with open(os.path.join(OUT, f"attrition_by_stage_{tag}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(stages[0].keys()))
        w.writeheader(); w.writerows(stages)

    # ---- 严格档候选清单 ----
    strict = sorted(s6, key=lambda r: (r["E_max"], -r["E_worst"]))
    fields = list(rows_combo[0].keys()) if rows_combo else []
    with open(os.path.join(OUT, f"strict_candidates_{tag}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=fields)
        w.writeheader(); w.writerows(strict)

    print(f"[{tag}] 数据齐全 {len(rows_combo)} / {n0}; 基线条件通过 {len(s5)}; 严格档 {len(s6)}; 缺失 {sum(n_missing.values())}")
    return rows_combo, s6


if __name__ == "__main__":
    run("Li2ZrCl6", "LZC")
    run("Li3YCl6", "LYC")

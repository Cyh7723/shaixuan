# -*- coding: utf-8 -*-
"""
Li3InCl6 路线后处理: 阈值敏感性 + MatTen/产物关联 + Pareto 排名 + 三 SSE 对比
"""
import os, csv, json
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
NEW = os.path.dirname(ROOT)
TOP = os.path.dirname(NEW)
OUT = os.path.join(NEW, "11_Li3InCl6路线")
SSE = "Li3InCl6"

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


def dominates(a, b, keys_dirs):
    better = False
    for k, d in keys_dirs:
        va, vb = a.get(k), b.get(k)
        if va is None or vb is None:
            return False
        if d > 0:
            if va < vb: return False
            if va > vb: better = True
        else:
            if va > vb: return False
            if va < vb: better = True
    return better


def main():
    rows = list(csv.DictReader(open(os.path.join(OUT, f"组合筛选_{SSE}.csv"), encoding="utf-8-sig")))
    ela = {r["mp_id"]: r for r in csv.DictReader(
        open(os.path.join(NEW, "02_MatTen力学", "material_elasticity.csv"), encoding="utf-8-sig"))}

    # 并入弹性 + 派生
    combos = []
    for r in rows:
        if r["group"] != "通过组":
            continue
        ga, gb = f(ela.get(r["A"], {}).get("GVRH")), f(ela.get(r["B"], {}).get("GVRH"))
        e_worst = min(f(r["Er_Li_A"]), f(r["Er_A_B"]), f(r["Er_B_SSE"]))
        e_max = f(r["E3D_max"])
        dv_min = min(f(r["dV_AB"]), f(r["dV_B_SSE"]))
        if None in (e_worst, e_max, dv_min):
            continue
        combos.append({
            "A": r["A"], "A_formula": r["A_formula"], "B": r["B"], "B_formula": r["B_formula"],
            "Eg_A": f(r["Eg_A"]), "Eg_B": f(r["Eg_B"]),
            "E_worst": e_worst, "E_max": e_max, "dV_min": dv_min,
            "GVRH_A": ga, "GVRH_B": gb,
            "GVRH_min": min([v for v in (ga, gb) if v is not None], default=None),
            "架构": ("双绝缘" if (f(r["Eg_A"]) or 0) >= 1 and (f(r["Eg_B"]) or 0) >= 1 else
                    ("绝缘A/导电B" if (f(r["Eg_A"]) or 0) >= 1 else "导电A/绝缘B")),
        })
    print(f"[{SSE}] 通过组带完整数据 {len(combos)}")

    # ---- 阈值敏感性 ----
    sens = []
    for erc in ER_CUTS:
        for e3c in E3D_CUTS:
            for egc in EG_CUTS:
                for dvc in DV_CUTS:
                    sel = [r for r in combos
                           if r["E_worst"] >= erc and r["E_max"] <= e3c
                           and (r["Eg_A"] >= egc or r["Eg_B"] >= egc)
                           and r["dV_min"] >= dvc]
                    sens.append({"SSE": SSE, "Er下限": erc, "E3D上限": e3c, "Eg门槛": egc, "dV下限": dvc,
                                 "配对数": len(sel),
                                 "独立A数": len({r["A"] for r in sel}),
                                 "独立B数": len({r["B"] for r in sel}),
                                 "独立化学式对数": len({(r["A_formula"], r["B_formula"]) for r in sel})})
    with open(os.path.join(OUT, f"threshold_sensitivity_{SSE}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(sens[0].keys())); w.writeheader(); w.writerows(sens)

    # ---- 严格档 ----
    strict = [r for r in combos if r["E_worst"] >= -0.03 and r["E_max"] <= 0.6]
    with open(os.path.join(OUT, f"strict_full_{SSE}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(combos[0].keys())); w.writeheader(); w.writerows(strict)
    print(f"[{SSE}] 严格档 {len(strict)}")

    # ---- Pareto ----
    keys_dirs = [("E_max", -1), ("E_worst", 1), ("dV_min", 1), ("GVRH_min", 1)]
    full = [r for r in strict if r["GVRH_min"] is not None]
    front = []
    for a in full:
        if not any(dominates(b, a, keys_dirs) for b in full if b is not a):
            front.append(a)
    front.sort(key=lambda r: (r["E_max"], -r["E_worst"]))
    with open(os.path.join(OUT, f"pareto_front_{SSE}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(combos[0].keys())); w.writeheader(); w.writerows(front)
    print(f"[{SSE}] Pareto 前沿 {len(front)}")

    # ---- 三 SSE 对比 ----
    lzc = {(r["A_formula"], r["B_formula"]) for r in csv.DictReader(
        open(os.path.join(NEW, "01_阈值敏感性", "strict_candidates_LZC.csv"), encoding="utf-8-sig"))}
    lyc = {(r["A_formula"], r["B_formula"]) for r in csv.DictReader(
        open(os.path.join(NEW, "01_阈值敏感性", "strict_candidates_LYC.csv"), encoding="utf-8-sig"))}
    lic = {(r["A_formula"], r["B_formula"]) for r in strict}
    triple = lzc & lyc & lic
    print(f"[三SSE] 严格档交集化学式对: LZC {len(lzc)}, LYC {len(lyc)}, Li3InCl6 {len(lic)}, 三者交集 {len(triple)}")

    # Pareto 交集
    plzc = {(r["A_formula"], r["B_formula"]) for r in csv.DictReader(
        open(os.path.join(NEW, "07_排名候选", "pareto_front_LZC.csv"), encoding="utf-8-sig"))}
    plyc = {(r["A_formula"], r["B_formula"]) for r in csv.DictReader(
        open(os.path.join(NEW, "07_排名候选", "pareto_front_LYC.csv"), encoding="utf-8-sig"))}
    plic = {(r["A_formula"], r["B_formula"]) for r in front}
    ptriple = plzc & plyc & plic
    print(f"[三SSE] Pareto 前沿交集: {len(ptriple)} -> {sorted(ptriple)}")

    with open(os.path.join(OUT, "三SSE对比.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.writer(fo)
        w.writerow(["类别", "LZC", "LYC", "Li3InCl6", "三者交集"])
        w.writerow(["严格档化学式对", len(lzc), len(lyc), len(lic), len(triple)])
        w.writerow(["Pareto前沿化学式对", len(plzc), len(plyc), len(plic), len(ptriple)])


if __name__ == "__main__":
    main()

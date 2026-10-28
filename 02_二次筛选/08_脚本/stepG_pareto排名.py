# -*- coding: utf-8 -*-
"""
阶段3: Pareto 排序 + 代表候选 + 新旧排名对比

流程: 必要条件 → 数据质量 → 功能路线标签 → Pareto → 代表候选

必要条件(严格档):
  E_worst >= -0.03, E_max <= 0.6, 至少一层 Eg>=1.0, dV_min >= 0.1 (原窗口规则)
数据质量分:
  BVSE 来源(查表>自算)、窗口完整、反应能完整、弹性正常(无质检标记)
功能路线(互斥标签, 优先级):
  R1 双绝缘+低反应+高模量 (GVRH_A>=35 且 GVRH_B>=35, 待核区 25-45 标"边界")
  R2 双绝缘+低反应+柔顺 (GVRH 较小)
  R3 合金辅助传输+明确电子阻隔 (一层金属间化合物, 另一层 Eg>=2)
  R4 受控反应成膜 (E_worst 在 [-0.10,-0.03), 产物无导电相)
Pareto 字段: E_max(↓), |E_worst|(↓), dV_min(↑), GVRH_min(↑)
  GVRH_min = min(GVRH_A, GVRH_B) — 双层有效模量不等于 max; 用 min 作保守指标。
缺失数据多的组合单独输出, 不进前沿。

输出:
  07_排名候选/pareto_front_LZC.csv, pareto_front_LYC.csv
  07_排名候选/route_candidates.csv     (分路线代表)
  07_排名候选/old_vs_new_ranking.csv   (新旧排名对比)
  07_排名候选/ranking_notes.md
"""
import os, csv, json
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
NEW = os.path.dirname(ROOT)
TOP = os.path.dirname(NEW)
WS = os.path.join(TOP, "AB池双层筛选")
OUT = os.path.join(NEW, "07_排名候选")
os.makedirs(OUT, exist_ok=True)

ER_STRICT = -0.03
E3D_STRICT = 0.6
G_HIGH = 35.0
G_GRAY_LO, G_GRAY_HI = 25.0, 45.0


def f(x):
    try:
        if x in ("", None, "inf"):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def dominates(a, b, keys_dirs):
    """a 支配 b: 所有目标不差且至少一个更好。keys_dirs: [(key, +1越大越好/-1越小越好)]"""
    better = False
    for k, d in keys_dirs:
        va, vb = a.get(k), b.get(k)
        if va is None or vb is None:
            return False
        if d > 0:
            if va < vb:
                return False
            if va > vb:
                better = True
        else:
            if va > vb:
                return False
            if va < vb:
                better = True
    return better


def run(sse, tag):
    strict_path = os.path.join(NEW, "01_阈值敏感性", f"strict_candidates_{tag}.csv")
    rows = list(csv.DictReader(open(strict_path, encoding="utf-8-sig")))
    # 并入弹性
    ela = {r["mp_id"]: r for r in csv.DictReader(
        open(os.path.join(NEW, "02_MatTen力学", "material_elasticity.csv"), encoding="utf-8-sig"))}
    # 并入反应产物(导电产物风险)
    rxn_risk = {}
    jpath = os.path.join(NEW, "03_反应产物", "interface_reactions.jsonl")
    if os.path.exists(jpath):
        for line in open(jpath, encoding="utf-8"):
            d = json.loads(line)
            if "error" in d:
                continue
            rxn_risk[(d["interface"], d["side1"], d["side2"])] = d
    # 并入单层对照判定
    bilayer = {}
    bp = os.path.join(NEW, "05_单双层对照", f"bilayer_vs_monolayer_{tag}.csv")
    if os.path.exists(bp):
        for r in csv.DictReader(open(bp, encoding="utf-8-sig")):
            bilayer[(r["A"], r["B"])] = r

    out = []
    for r in rows:
        aid, bid = r["A"], r["B"]
        ea, eb = ela.get(aid, {}), ela.get(bid, {})
        ga, gb = f(ea.get("GVRH")), f(eb.get("GVRH"))
        e_worst = f(r["E_worst"]); e_max = f(r["E_max"]); dv_min = f(r["dV_min"])
        if None in (e_worst, e_max, dv_min):
            continue
        # 导电产物风险(三界面)
        cond = []
        for lab, s1, s2 in [("Li|A", "Li_metal", aid), ("A|B", aid, bid), ("B|SSE", bid, sse)]:
            rec = rxn_risk.get((lab, s1, s2))
            if rec:
                cond.extend(rec.get("conducting_products", []))
        cond = sorted(set(c for c in cond if c))
        # 缺失数据计数
        miss = sum(1 for v in (ga, gb) if v is None)
        g_min = min([v for v in (ga, gb) if v is not None], default=None)
        # 路线
        arch = r["架构"]
        if arch == "双绝缘" and g_min is not None and g_min >= G_HIGH:
            route = "R1_双绝缘低反应高模量"
        elif arch == "双绝缘":
            route = "R2_双绝缘低反应柔顺"
        elif "导电" in arch:
            route = "R3_合金辅助+电子阻隔"
        else:
            route = "R4_其他"
        if cond:
            route += "|导电产物风险"
        gray = (ga is not None and G_GRAY_LO <= ga <= G_GRAY_HI) or \
               (gb is not None and G_GRAY_LO <= gb <= G_GRAY_HI)
        bl = bilayer.get((aid, bid), {})
        out.append({
            "SSE": sse, "A": aid, "A_formula": r["A_formula"], "B": bid, "B_formula": r["B_formula"],
            "架构": arch, "路线": route,
            "E_worst": e_worst, "E_max": e_max, "dV_min": dv_min,
            "GVRH_A": ga, "GVRH_B": gb, "GVRH_min": g_min,
            "GVRH_灰区": int(gray),
            "导电产物": ";".join(cond),
            "双层判定": bl.get("判定", ""), "双层判定说明": bl.get("判定说明", ""),
            "缺失计数": miss,
            "A弹性质检": ea.get("质量标记", ""), "B弹性质检": eb.get("质量标记", ""),
        })

    # Pareto 前沿(只在缺失计数==0 的组合上)
    keys_dirs = [("E_max", -1), ("E_worst", 1), ("dV_min", 1), ("GVRH_min", 1)]
    full = [r for r in out if r["缺失计数"] == 0 and r["GVRH_min"] is not None]
    front = []
    for i, a in enumerate(full):
        if not any(dominates(b, a, keys_dirs) for b in full if b is not a):
            front.append(a)
    front.sort(key=lambda r: (r["E_max"], -r["E_worst"]))
    fields = list(out[0].keys())
    with open(os.path.join(OUT, f"pareto_front_{tag}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=fields); w.writeheader(); w.writerows(front)

    # 分路线代表: 每路线取 Pareto 前沿内 E_max 最小的前 5, 同 A 化学式限 2
    reps = []
    by_route = defaultdict(list)
    for r in front:
        by_route[r["路线"].split("|")[0]].append(r)
    for route, lst in sorted(by_route.items()):
        seen_a = Counter()
        n = 0
        for r in lst:
            if seen_a[r["A_formula"]] >= 2:
                continue
            seen_a[r["A_formula"]] += 1
            reps.append({**r, "代表备注": f"{route} 第{n+1}名"})
            n += 1
            if n >= 5:
                break
    with open(os.path.join(OUT, f"route_candidates_{tag}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(reps[0].keys()) if reps else fields)
        w.writeheader(); w.writerows(reps)

    # 新旧排名对比
    old_rank = {}
    op = os.path.join(WS, "05_结果", f"排名评分_{sse}.csv")
    if os.path.exists(op):
        for i, r in enumerate(csv.DictReader(open(op, encoding="utf-8-sig"))):
            old_rank[(r["A"], r["B"])] = i + 1
    new_rank = { (r["A"], r["B"]): i+1 for i, r in enumerate(front) }
    cmp_rows = []
    for r in out:
        key = (r["A"], r["B"])
        cmp_rows.append({
            "A": r["A_formula"], "B": r["B_formula"], "SSE": sse,
            "旧名次": old_rank.get(key, ""), "新Pareto名次": new_rank.get(key, ""),
            "变化原因": "" if key in old_rank else "旧排名外(未进通过组前30或新筛出)",
        })
    with open(os.path.join(OUT, f"old_vs_new_{tag}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(cmp_rows[0].keys())); w.writeheader(); w.writerows(cmp_rows)

    print(f"[{tag}] 严格档 {len(out)} 完整数据 {len(full)} Pareto前沿 {len(front)} 代表 {len(reps)}")
    return out, front, reps


if __name__ == "__main__":
    run("Li2ZrCl6", "LZC")
    run("Li3YCl6", "LYC")

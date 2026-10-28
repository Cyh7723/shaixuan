# -*- coding: utf-8 -*-
"""
规则7(扩展): 通过组的优先级分档 + 多指标排名 + 短名单

不改变宽筛通过组, 只在其上加一个更严格的优先档(试运行阈值, 尚未认定为最佳):

  指标            宽筛通过条件          优先档
  E_3D,max        <= 1.0 eV             <= 0.6 eV      传输潜力优先
  E_worst         >= -0.10 eV/atom      >= -0.03       减少对反应钝化的依赖
  ΔV_min          > 0 且避开边界敏感    >= 0.2 V       窗口匹配余量
  电子阻隔        至少一层 Eg>=1.0      保持原要求
  Ehull           <= 5 meV/atom         保持原要求

  E_worst = min(Er_Li,A, Er_A,B, Er_B,SSE)
  ΔV_min  = min(ΔV_AB, ΔV_B,SSE)

排名(优先档内, 综合分 0-100):
  s_transport = max(0, 1 - E_3D,max/0.6)        权重 0.45
  s_rxn       = clip((E_worst+0.03)/0.03, 0,1)  权重 0.35
  s_window    = min(ΔV_min, 2.0)/2.0            权重 0.20
  权重为试运行值, 报告中注明可调整。

去重与限额:
  - 同化学式不同晶型: 保留晶型(mp_id)信息, 选综合分最高的代表, 其余列入"晶型备选";
  - A|B 与 B|A 层序不同, 分别保留(A/B 角色由池定义);
  - 短名单每个 A 化学式限额(默认 5), 避免被 LiCl 组合占满。

输出:
  05_结果/优先档_Li2ZrCl6.csv, 优先档_Li3YCl6.csv   (优先档全量, 按综合分降序)
  05_结果/短名单_双SSE_每A限额5.csv                 (两 SSE 均入优先档, 每 A 限额)
"""
import os, csv

ROOT = os.path.dirname(os.path.abspath(__file__))
WS = os.path.dirname(ROOT)

E3D_PRIO = 0.6
EWORST_PRIO = -0.03
DV_PRIO = 0.2
W_T, W_R, W_W = 0.45, 0.35, 0.20
QUOTA_PER_A = 5


def f(x):
    try:
        if x in ("", None, "inf"):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def score(e3d_max, e_worst, dv_min):
    s_t = max(0.0, 1.0 - e3d_max / E3D_PRIO)
    s_r = min(1.0, max(0.0, (e_worst - EWORST_PRIO) / (-EWORST_PRIO)))
    s_w = min(dv_min, 2.0) / 2.0
    return round(100 * (W_T * s_t + W_R * s_r + W_W * s_w), 1)


def load_pass(sse):
    rows = []
    for r in csv.DictReader(open(os.path.join(WS, "05_结果", f"组合筛选_{sse}.csv"), encoding="utf-8-sig")):
        if r["group"] != "通过组":
            continue
        e3d = f(r["E3D_max"]); er_ab = f(r["Er_A_B"]); er_li = f(r["Er_Li_A"]); er_bs = f(r["Er_B_SSE"])
        dv_ab = f(r["dV_AB"]); dv_bs = f(r["dV_B_SSE"])
        if None in (e3d, er_ab, er_li, er_bs, dv_ab, dv_bs):
            continue
        r["_E3D_max"] = e3d
        r["_E_worst"] = min(er_li, er_ab, er_bs)
        r["_dV_min"] = min(dv_ab, dv_bs)
        rows.append(r)
    return rows


def prio_tier(rows):
    return [r for r in rows if r["_E3D_max"] <= E3D_PRIO
            and r["_E_worst"] >= EWORST_PRIO and r["_dV_min"] >= DV_PRIO]


def dedup_polymorph(rows):
    """同 (A化学式, B化学式) 选综合分最高的晶型代表, 其余标为晶型备选。"""
    best = {}
    for r in rows:
        k = (r["A_formula"], r["B_formula"])
        if k not in best or r["_score"] > best[k]["_score"]:
            best[k] = r
    reps = set(id(v) for v in best.values())
    for r in rows:
        r["晶型标记"] = "代表" if id(r) in reps else "晶型备选"
    return rows


def main():
    per_sse = {}
    for sse in ["Li2ZrCl6", "Li3YCl6"]:
        rows = load_pass(sse)
        for r in rows:
            r["_score"] = score(r["_E3D_max"], r["_E_worst"], r["_dV_min"])
        pr = prio_tier(rows)
        pr = dedup_polymorph(pr)
        pr.sort(key=lambda r: -r["_score"])
        reps = [r for r in pr if r["晶型标记"] == "代表"]
        per_sse[sse] = (rows, pr, reps)
        print(f"[{sse}] 通过组 {len(rows)} → 优先档 {len(pr)} (去重化学式对后 {len(reps)})")

        fields = ["排名", "A_formula", "B_formula", "A", "B", "综合分",
                  "E3D_max", "E3D_A", "E3D_B", "Er_Li_A", "Er_A_B", "Er_B_SSE", "E_worst",
                  "dV_AB", "dV_B_SSE", "dV_min", "Eg_A", "Eg_B", "层序标记", "晶型标记"]
        with open(os.path.join(WS, "05_结果", f"优先档_{sse}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
            w = csv.DictWriter(fo, fieldnames=fields)
            w.writeheader()
            for i, r in enumerate(reps, 1):
                w.writerow({
                    "排名": i, "A_formula": r["A_formula"], "B_formula": r["B_formula"],
                    "A": r["A"], "B": r["B"], "综合分": r["_score"],
                    "E3D_max": r["_E3D_max"], "E3D_A": r["E3D_A"], "E3D_B": r["E3D_B"],
                    "Er_Li_A": r["Er_Li_A"], "Er_A_B": r["Er_A_B"], "Er_B_SSE": r["Er_B_SSE"],
                    "E_worst": round(r["_E_worst"], 4),
                    "dV_AB": r["dV_AB"], "dV_B_SSE": r["dV_B_SSE"], "dV_min": round(r["_dV_min"], 2),
                    "Eg_A": r["Eg_A"], "Eg_B": r["Eg_B"],
                    "层序标记": r["层序标记"], "晶型标记": r["晶型标记"]})

    # 双 SSE 交集: 两路线都入优先档(按代表化学式对)
    reps_lzc = {(r["A_formula"], r["B_formula"]): r for r in per_sse["Li2ZrCl6"][2]}
    reps_lyc = {(r["A_formula"], r["B_formula"]): r for r in per_sse["Li3YCl6"][2]}
    both = sorted(set(reps_lzc) & set(reps_lyc),
                  key=lambda k: -(reps_lzc[k]["_score"] + reps_lyc[k]["_score"]))
    print(f"双SSE均入优先档(去重化学式对): {len(both)}")

    # 短名单: 每个 A 限额
    count_a = {}
    shortlist = []
    for k in both:
        fa = k[0]
        if count_a.get(fa, 0) >= QUOTA_PER_A:
            continue
        count_a[fa] = count_a.get(fa, 0) + 1
        r1, r2 = reps_lzc[k], reps_lyc[k]
        shortlist.append({
            "A": k[0], "B": k[1],
            "综合分_LZC": r1["_score"], "综合分_LYC": r2["_score"],
            "E3D_max": r1["_E3D_max"],
            "E_worst_LZC": round(r1["_E_worst"], 4), "E_worst_LYC": round(r2["_E_worst"], 4),
            "dV_min_LZC": round(r1["_dV_min"], 2), "dV_min_LYC": round(r2["_dV_min"], 2),
            "Er_A_B": r1["Er_A_B"], "Eg_A": r1["Eg_A"], "Eg_B": r1["Eg_B"],
            "层序标记": r1["层序标记"],
        })
    with open(os.path.join(WS, "05_结果", f"短名单_双SSE_每A限额{QUOTA_PER_A}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(shortlist[0].keys()))
        w.writeheader(); w.writerows(shortlist)
    print(f"短名单(每A限额{QUOTA_PER_A}): {len(shortlist)} 对, 覆盖 A 材料 {len(count_a)} 种")
    print("--- 短名单前 15 ---")
    for r in shortlist[:15]:
        print("  {:12s}|{:16s} 分LZC={} 分LYC={} E3D={} Eworst={}/{}".format(
            r["A"], r["B"], r["综合分_LZC"], r["综合分_LYC"], r["E3D_max"],
            r["E_worst_LZC"], r["E_worst_LYC"]))


if __name__ == "__main__":
    main()

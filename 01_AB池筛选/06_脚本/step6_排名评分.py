# -*- coding: utf-8 -*-
"""
通过组多性质排名评分(规则7 排序的量化实现)

方法: Aykol et al., Nat. Commun. 2016 (电池保护涂层高通量筛选)的
Weighted sum 一路——各性质 min-max 归一化到 0-1, 统一"越大越好",
等权合成总分 S = (f_E + f_R + f_V)/3, 总分越高越优。
(原文另有 Rank aggregation 一路, 本管线不采用。)

指标(当前数据最完整的三项, 均来自通过组已有列):
  E_max   = max(E3D_A, E3D_B)             越低越优 (规则6 传输)
  E_worst = min(Er_Li_A, Er_A_B, Er_B_SSE) 越接近 0 越优 (规则4 三界面)
  dV_min  = min(dV_AB, dV_B_SSE)          越大越优 (规则5 窗口余量)

权重敏感性: 在 (E_max, E_worst, dV_min) 权重单形上以 0.25 步长采样
15 组正权重, 记录各组合名次的最好/最差区间, 检验前列名单稳定性。

输出(每个 SSE 一份):
  05_结果/排名评分_<sse>.csv   含综合分与名次
  05_结果/双SSE交集排名.csv    两 SSE 均通过组合的排名对照
"""
import os, csv
from collections import defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
WS = os.path.dirname(ROOT)

METRICS = ["E_max", "E_worst", "dV_min"]
# 方向: +1 越大越优, -1 越小越优
DIRECTION = {"E_max": -1, "E_worst": +1, "dV_min": +1}


def fnum(x):
    try:
        if x in ("", None):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def load_pass(sse):
    rows = []
    for r in csv.DictReader(open(os.path.join(WS, "05_结果", f"组合筛选_{sse}.csv"), encoding="utf-8-sig")):
        if r["group"] != "通过组":
            continue
        e3d_a, e3d_b = fnum(r["E3D_A"]), fnum(r["E3D_B"])
        er_li, er_ab, er_bs = fnum(r["Er_Li_A"]), fnum(r["Er_A_B"]), fnum(r["Er_B_SSE"])
        dv_ab, dv_bs = fnum(r["dV_AB"]), fnum(r["dV_B_SSE"])
        if None in (e3d_a, e3d_b, er_li, er_ab, er_bs, dv_ab, dv_bs):
            continue  # 通过组原则上数据完整, 防御性跳过
        rows.append({
            "A": r["A"], "B": r["B"], "A_formula": r["A_formula"], "B_formula": r["B_formula"],
            "E_max": max(e3d_a, e3d_b),
            "E_worst": min(er_li, er_ab, er_bs),
            "dV_min": min(dv_ab, dv_bs),
            "Eg_A": r["Eg_A"], "Eg_B": r["Eg_B"],
            "Er_Li_A": er_li, "Er_A_B": er_ab, "Er_B_SSE": er_bs,
            "dV_AB": dv_ab, "dV_B_SSE": dv_bs,
            "层序标记": r["层序标记"],
        })
    return rows


def ranks_desc(values):
    """平均秩, 数值大者名次靠前(第 1 名最优)。"""
    order = sorted(range(len(values)), key=lambda i: values[i], reverse=True)
    rk = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            rk[order[k]] = avg
        i = j + 1
    return rk


def score(rows, weights):
    """weights: dict metric->w。返回总分列表(0-1, 越大越优)。"""
    n = len(rows)
    norm = {}
    for m in METRICS:
        vals = [r[m] for r in rows]
        lo, hi = min(vals), max(vals)
        if hi == lo:
            norm[m] = [0.5] * n   # 无区分度, 置中性
        else:
            if DIRECTION[m] > 0:
                norm[m] = [(v - lo) / (hi - lo) for v in vals]
            else:
                norm[m] = [(hi - v) / (hi - lo) for v in vals]
    wsum = sum(weights.values())
    return [sum(weights[m] * norm[m][i] for m in METRICS) / wsum for i in range(n)]


def main(sse):
    rows = load_pass(sse)
    print(f"[{sse}] 通过组可评分组合 {len(rows)}")

    sa = score(rows, {m: 1.0 for m in METRICS})
    ra = ranks_desc(sa)
    assert abs(sa[ra.index(1.0)] - max(sa)) < 1e-9, "名次方向错误"
    for i, r in enumerate(rows):
        r["score"] = round(sa[i], 4)
        r["rank"] = ra[i]

    # 权重敏感性: 单形上网格采样(步长 0.25), 记录名次区间
    grid = [w / 4 for w in range(5)]
    weight_sets = [(we, wr, wv) for we in grid for wr in grid for wv in grid
                   if abs(we + wr + wv - 1.0) < 1e-9 and we > 0 and wr > 0 and wv > 0]
    rank_span = defaultdict(lambda: [10**9, -1])
    for we, wr, wv in weight_sets:
        ra_w = ranks_desc(score(rows, {"E_max": we, "E_worst": wr, "dV_min": wv}))
        for i in range(len(rows)):
            k = (rows[i]["A"], rows[i]["B"])
            rank_span[k][0] = min(rank_span[k][0], ra_w[i])
            rank_span[k][1] = max(rank_span[k][1], ra_w[i])
    for r in rows:
        lo, hi = rank_span[(r["A"], r["B"])]
        r["rank_min_w"], r["rank_max_w"] = lo, hi

    rows.sort(key=lambda r: r["rank"])
    fields = ["rank", "score", "A_formula", "B_formula", "A", "B",
              "E_max", "E_worst", "dV_min",
              "Eg_A", "Eg_B", "Er_Li_A", "Er_A_B", "Er_B_SSE", "dV_AB", "dV_B_SSE",
              "rank_min_w", "rank_max_w", "层序标记"]
    with open(os.path.join(WS, "05_结果", f"排名评分_{sse}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        wcsv = csv.DictWriter(fo, fieldnames=fields)
        wcsv.writeheader(); wcsv.writerows(rows)
    stable_top30 = sum(1 for r in rows[:30] if r["rank_max_w"] <= 30)
    print(f"[{sse}] 前30名中权重扰动下始终在前30: {stable_top30}/30")
    return rows


if __name__ == "__main__":
    all_rows = {}
    for sse in ["Li2ZrCl6", "Li3YCl6"]:
        all_rows[sse] = load_pass(sse)

    # 以化学式对为键取最优代表(E_max 最小)
    for sse in all_rows:
        best = {}
        for r in all_rows[sse]:
            k = (r["A_formula"], r["B_formula"])
            if k not in best or r["E_max"] < best[k]["E_max"]:
                best[k] = r
        all_rows[sse] = list(best.values())

    scored = {sse: main(sse) for sse in all_rows}

    # 双 SSE 交集对照表(按两路线名次最大值升序)
    lzc = {(r["A_formula"], r["B_formula"]): r for r in scored["Li2ZrCl6"]}
    lyc = {(r["A_formula"], r["B_formula"]): r for r in scored["Li3YCl6"]}
    common = sorted(set(lzc) & set(lyc),
                    key=lambda k: max(lzc[k]["rank"], lyc[k]["rank"]))
    out = []
    for k in common:
        r1, r2 = lzc[k], lyc[k]
        out.append({"排名": len(out) + 1, "A": k[0], "B": k[1],
                    "score_LZC": r1["score"], "score_LYC": r2["score"],
                    "rank_LZC": r1["rank"], "rank_LYC": r2["rank"],
                    "E_max": r1["E_max"], "E_worst_min": min(r1["E_worst"], r2["E_worst"]),
                    "dV_min_min": min(r1["dV_min"], r2["dV_min"]),
                    "rank_span_LZC": f"{r1['rank_min_w']}-{r1['rank_max_w']}",
                    "rank_span_LYC": f"{r2['rank_min_w']}-{r2['rank_max_w']}"})
    with open(os.path.join(WS, "05_结果", "双SSE交集排名.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        wcsv = csv.DictWriter(fo, fieldnames=list(out[0].keys()))
        wcsv.writeheader(); wcsv.writerows(out)
    print(f"双SSE交集(去重化学式对) {len(common)} 对, 已写 双SSE交集排名.csv")

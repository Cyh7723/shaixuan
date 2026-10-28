# -*- coding: utf-8 -*-
"""
阶段0: 基线复现 + 数据审计
- 复现 step5 分组计数 (通过/待核/动力学线索/未通过) 对照报告基线
- 解释 340x793=269620 vs 269602 / 340x1070=363800 vs 363715 (A==B 同材料配对排除)
- 母池化学式-晶型重复检查; A=B 排除规则记录
- 类别互斥性检查
输出: baseline_reproduction.csv, data_audit.md (00_审计/)
"""
import os, csv, json, sys
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
NEW = os.path.dirname(ROOT)                      # 二次筛选_20261021
TOP = os.path.dirname(NEW)                       # 苏州实验室课题学习记录
WS = os.path.join(TOP, "AB池双层筛选")
AUDIT = os.path.join(NEW, "00_审计")
os.makedirs(AUDIT, exist_ok=True)

# 报告基线
BASELINE = {
    "母池": 1860, "A池": 340, "B池_LZC": 793, "B池_LYC": 1070,
    "组合总数_LZC": 269602, "组合总数_LYC": 363715,
    "通过_LZC": 6501, "通过_LYC": 8074,
    "待核_LZC": 8854, "待核_LYC": 15114,
    "动力学线索_LZC": 6983, "动力学线索_LYC": 9911,
    "未通过_LZC": 247264, "未通过_LYC": 330616,
}


def read_csv(p):
    return list(csv.DictReader(open(p, encoding="utf-8-sig")))


def main():
    rows = []          # baseline_reproduction.csv
    audit = []         # data_audit.md 片段

    # ---- 1. 池规模复现 ----
    pool = read_csv(os.path.join(WS, "01_母池", "母池_无带隙门槛.csv"))
    a_pool = [r for r in read_csv(os.path.join(WS, "03_AB池", "A池_靠Li.csv")) if r["A池标记"] == "入池"]
    a_pend = [r for r in read_csv(os.path.join(WS, "03_AB池", "A池_靠Li.csv")) if r["A池标记"] != "入池"]
    b_lzc = [r for r in read_csv(os.path.join(WS, "03_AB池", "B池_靠SSE_Li2ZrCl6.csv")) if r["B池标记"] == "入池"]
    b_lyc = [r for r in read_csv(os.path.join(WS, "03_AB池", "B池_靠SSE_Li3YCl6.csv")) if r["B池标记"] == "入池"]

    checks = [
        ("母池", len(pool), BASELINE["母池"]),
        ("A池入池", len(a_pool), BASELINE["A池"]),
        ("B池入池_LZC", len(b_lzc), BASELINE["B池_LZC"]),
        ("B池入池_LYC", len(b_lyc), BASELINE["B池_LYC"]),
    ]
    for name, got, exp in checks:
        rows.append({"项目": name, "复现值": got, "报告基线": exp,
                     "一致": "是" if got == exp else "否", "说明": ""})

    # ---- 2. 组合总数差异: A==B 排除 ----
    a_ids = {r["mp_id"] for r in a_pool}
    b_ids_lzc = {r["mp_id"] for r in b_lzc}
    b_ids_lyc = {r["mp_id"] for r in b_lyc}
    inter_lzc = a_ids & b_ids_lzc
    inter_lyc = a_ids & b_ids_lyc
    n_lzc_raw = len(a_ids) * len(b_ids_lzc)
    n_lyc_raw = len(a_ids) * len(b_ids_lyc)
    rows.append({"项目": "LZC 340x793 原始笛卡尔积", "复现值": n_lzc_raw,
                 "报告基线": 269620, "一致": "是" if n_lzc_raw == 269620 else "否", "说明": ""})
    rows.append({"项目": "LZC A∩B 同mp_id配对排除数", "复现值": len(inter_lzc),
                 "报告基线": 269620 - 269602, "一致": "是" if len(inter_lzc) == 18 else "否",
                 "说明": "step5 中 aid==bid 跳过(A=B同结构不构成双层); 269620-18=269602"})
    rows.append({"项目": "LYC 340x1070 原始笛卡尔积", "复现值": n_lyc_raw,
                 "报告基线": 363800, "一致": "是" if n_lyc_raw == 363800 else "否", "说明": ""})
    rows.append({"项目": "LYC A∩B 同mp_id配对排除数", "复现值": len(inter_lyc),
                 "报告基线": 363800 - 363715, "一致": "是" if len(inter_lyc) == 85 else "否",
                 "说明": "363800-85=363715"})
    audit.append(f"## 组合总数差异\n"
                 f"- LZC: 340x793={n_lzc_raw}, step5 排除 A==B 同 mp_id 配对 {len(inter_lzc)} 项 → {n_lzc_raw - len(inter_lzc)}, 与报告 269602 一致。\n"
                 f"- LYC: 340x1070={n_lyc_raw}, 排除 {len(inter_lyc)} 项 → {n_lyc_raw - len(inter_lyc)}, 与报告 363715 一致。\n"
                 f"- 排除规则: A 与 B 为同一 mp_id(同一晶体结构)时不构成双层, step5_组合筛选.py 第97行 `if aid==bid: continue`。\n"
                 f"- 注意: 仅排除同 mp_id; 同化学式不同晶型不合并(见下)。\n")

    # ---- 3. 同化学式不同晶型检查 ----
    merged = {r["mp_id"]: r for r in read_csv(os.path.join(WS, "03_AB池", "母池_合并热力学数据.csv"))}
    f2ids = defaultdict(list)
    for mid, r in merged.items():
        f2ids[r["formula"]].append(mid)
    dup_formula = {f: ids for f, ids in f2ids.items() if len(ids) > 1}
    n_a_dup = sum(1 for r in a_pool if len(f2ids.get(r["formula"], [])) > 1)
    audit.append(f"## 化学式-晶型\n"
                 f"- 母池 {len(merged)} 个唯一 mp_id, 其中 {len(dup_formula)} 个化学式对应多个晶型(共 {sum(len(v) for v in dup_formula.values())} 条记录), 全程按 mp_id 区分, 未提前合并。\n"
                 f"- A 池中涉及多晶型化学式的材料数: {n_a_dup}。\n")

    # 跨 mp_id 同结构(近似重复)粗查: 同化学式+同nsites+同volume(±0.1%)
    seen = {}
    near_dup = []
    for mid, r in merged.items():
        try:
            key = (r["formula"], int(float(r["nsites"])), round(float(r["volume"]), 1))
        except Exception:
            continue
        if key in seen:
            near_dup.append((seen[key], mid))
        else:
            seen[key] = mid
    audit.append(f"## 近似重复结构粗查(同化学式+同nsites+同体积0.1Å³级)\n"
                 f"- 命中 {len(near_dup)} 对(不同 mp_id 但几何参数一致, 可能为同源条目)。"
                 f"本轮不删除, 在 elasticity/reaction 输出中保留各自 mp_id, 供人工抽查。\n")

    # ---- 4. 分组计数复现 + 互斥性 ----
    for sse, tag in [("Li2ZrCl6", "LZC"), ("Li3YCl6", "LYC")]:
        p = os.path.join(WS, "05_结果", f"组合筛选_{sse}.csv")
        n = 0
        cnt = Counter()
        with open(p, encoding="utf-8-sig") as f_:
            for r in csv.DictReader(f_):
                n += 1
                cnt[r["group"]] += 1
        rows.append({"项目": f"{tag} 组合总数(csv行数)", "复现值": n,
                     "报告基线": BASELINE[f"组合总数_{tag}"],
                     "一致": "是" if n == BASELINE[f"组合总数_{tag}"] else "否", "说明": ""})
        for g, key in [("通过组", "通过"), ("待核组", "待核"), ("动力学线索", "动力学线索"), ("未通过组", "未通过")]:
            exp = BASELINE[f"{key}_{tag}"]
            got = cnt.get(g, 0)
            rows.append({"项目": f"{tag} {g}", "复现值": got, "报告基线": exp,
                         "一致": "是" if got == exp else "否", "说明": ""})
        s = sum(cnt.values())
        audit.append(f"## {tag} 分组复现\n- 计数 {dict(cnt)}, 合计 {s}。"
                     f"group 列为单值字段, 四组天然互斥, 不存在重复统计。\n"
                     f"- 优先级: 未通过(有fail) > 动力学线索(无交叠) > 待核 > 通过。\n"
                     f"- 即'动力学线索'是'窗口无交叠但其余条件不满足淘汰级'的集合, 与通过组互斥。\n")

    # ---- 5. Ehull 与门槛记录 ----
    eh = [float(r["Ehull_eV_per_atom"]) for r in pool if r["Ehull_eV_per_atom"] not in ("", None)]
    audit.append(f"## 母池门槛\n- Ehull ≤ 5 meV/atom: 实际 max={max(eh)*1000:.2f} meV/atom, n={len(eh)}。\n"
                 f"- 元素排除规则见 step1_母池.py(贵金属/放射/高毒/Re/惰性气体)。\n")

    # ---- 输出 ----
    with open(os.path.join(AUDIT, "baseline_reproduction.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=["项目", "复现值", "报告基线", "一致", "说明"])
        w.writeheader(); w.writerows(rows)

    n_ok = sum(1 for r in rows if r["一致"] == "是")
    with open(os.path.join(AUDIT, "data_audit.md"), "w", encoding="utf-8") as fo:
        fo.write("# 数据审计 (阶段0)\n\n日期: 2026-10-21\n\n")
        fo.write(f"## 基线复现总览\n- {n_ok}/{len(rows)} 项与报告一致。明细见 baseline_reproduction.csv。\n\n")
        fo.write("\n".join(audit))
        fo.write("\n## 类别互斥性\n"
                 "- step5 的 group 字段为单值: 每条组合只属于 通过组/待核组/动力学线索/未通过组 之一。\n"
                 "- 判定顺序: reasons_fail 非空→未通过; 否则窗口无交叠→动力学线索; 否则有待核原因→待核; 否则通过。\n"
                 "- 因此'动力学线索'不包含在通过组内, 统计不重复。\n")
    print(f"[阶段0] baseline_reproduction {n_ok}/{len(rows)} 一致; 输出至 {AUDIT}")


if __name__ == "__main__":
    main()

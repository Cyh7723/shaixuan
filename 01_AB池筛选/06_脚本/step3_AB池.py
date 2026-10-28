# -*- coding: utf-8 -*-
"""
规则2: 按接触对象建立 A 池(靠 Li)与 B 池(靠 SSE)

A 池(Li 侧):
  - E_rxn(Li,A) >= -0.10 eV/atom (必需)
  - V_red^A <= 0.3 V (暂时保留; 0.3 V 是对部分还原驱动力的经验容忍度,
    不得表述为"在 0 V 严格稳定")
  - 电子绝缘推迟到组合阶段判断, 此步不限带隙
B 池(SSE 侧):
  - E_rxn(B,SSE) >= -0.10 eV/atom (必需; SSE = Li2ZrCl6_beta 与 Li3YCl6 分别成池)
  - 电子绝缘同样推迟到组合阶段
数据缺失(status 非 ok / chemspace>6)者进"待核"标记, 不视为淘汰。

输入: 01_母池/母池_无带隙门槛.csv, 02_热力学数据/窗口_checkpoint.jsonl,
      02_热力学数据/反应能_checkpoint.jsonl
输出: 03_AB池/A池_靠Li.csv, B池_靠SSE_Li2ZrCl6.csv, B池_靠SSE_Li3YCl6.csv,
      03_AB池/母池_合并热力学数据.csv (供后续组合筛选)
"""
import os, csv, json

ROOT = os.path.dirname(os.path.abspath(__file__))
WS = os.path.dirname(ROOT)

ER_TOL = -0.10          # eV/atom, 相容性容忍阈值
VRED_A_MAX = 0.3        # A 层还原电位上限(经验容忍, 非 0V 严格稳定)


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    pool = list(csv.DictReader(open(os.path.join(WS, "01_母池", "母池_无带隙门槛.csv"), encoding="utf-8-sig")))
    win = {}
    for line in open(os.path.join(WS, "02_热力学数据", "窗口_checkpoint.jsonl"), encoding="utf-8"):
        d = json.loads(line)
        win[d["mp_id"]] = d
    rxn = {}
    for line in open(os.path.join(WS, "02_热力学数据", "反应能_checkpoint.jsonl"), encoding="utf-8"):
        d = json.loads(line)
        rxn[d["mp_id"]] = d

    rows = []
    for r in pool:
        m = dict(r)
        w = win.get(r["mp_id"], {})
        x = rxn.get(r["mp_id"], {})
        m["V_red"] = w.get("V_red", "")
        m["V_ox"] = w.get("V_ox", "")
        m["window_method"] = w.get("method", "")
        m["win_status"] = w.get("status", "missing")
        m["Er_Li"] = x.get("Er_Li", "")
        m["Er_Li2ZrCl6"] = x.get("Er_Li2ZrCl6", "")
        m["Er_Li3YCl6"] = x.get("Er_Li3YCl6", "")
        m["rxn_method_vs_Li"] = x.get("method_vs_Li", "")
        m["rxn_method_vs_SSE"] = x.get("method_vs_SSE", "")
        m["rxn_status"] = x.get("status", "missing")
        rows.append(m)

    fields = list(rows[0].keys())
    with open(os.path.join(WS, "03_AB池", "母池_合并热力学数据.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        wcsv = csv.DictWriter(fo, fieldnames=fields)
        wcsv.writeheader(); wcsv.writerows(rows)

    # ---- A 池 ----
    a_pass, a_pending = [], []
    for m in rows:
        er = f(m["Er_Li"]); vred = f(m["V_red"])
        if er is None or m["rxn_status"] != "ok":
            mm = dict(m); mm["A池标记"] = "待核:Er_Li缺失"
            a_pending.append(mm); continue
        if er < ER_TOL:
            continue  # 不满足 A 池必需条件, 不进池
        if vred is None or m["win_status"] != "ok":
            mm = dict(m); mm["A池标记"] = "待核:窗口缺失"
            a_pending.append(mm); continue
        if vred > VRED_A_MAX:
            continue
        mm = dict(m); mm["A池标记"] = "入池"
        a_pass.append(mm)

    # ---- B 池(每个 SSE 一列池) ----
    b_stats = {}
    for sse in ["Li2ZrCl6", "Li3YCl6"]:
        col = f"Er_{sse}"
        b_pass, b_pending = [], []
        for m in rows:
            er = f(m[col])
            if er is None:
                note = m.get(f"note_{sse}", "") or ("rxn_status=" + m["rxn_status"])
                mm = dict(m); mm["B池标记"] = f"待核:{note or col+'缺失'}"
                b_pending.append(mm); continue
            if er < ER_TOL:
                continue
            mm = dict(m); mm["B池标记"] = "入池"
            b_pass.append(mm)
        b_stats[sse] = (len(b_pass), len(b_pending))
        out = b_pass + b_pending
        with open(os.path.join(WS, "03_AB池", f"B池_靠SSE_{sse}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
            wcsv = csv.DictWriter(fo, fieldnames=fields + ["B池标记"])
            wcsv.writeheader(); wcsv.writerows(out)

    with open(os.path.join(WS, "03_AB池", "A池_靠Li.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        wcsv = csv.DictWriter(fo, fieldnames=fields + ["A池标记"])
        wcsv.writeheader(); wcsv.writerows(a_pass + a_pending)

    print(f"A池(靠Li): 入池 {len(a_pass)}, 待核 {len(a_pending)}")
    for sse, (p, pd_) in b_stats.items():
        print(f"B池(靠{sse}): 入池 {p}, 待核 {pd_}")


if __name__ == "__main__":
    main()

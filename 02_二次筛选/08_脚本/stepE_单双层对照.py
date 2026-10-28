# -*- coding: utf-8 -*-
"""
工作包E: 单层 vs 双层对照

对每对 (A,B,SSE) 比较三种结构:
  Li|A|SSE     (A 兼作 A 层与 B 层功能)
  Li|B|SSE     (B 直接接触 Li)
  Li|A|B|SSE   (双层)

比较量: 最差界面反应能、电子阻隔、窗口(电稳定性)、传输证据、缺失数据。
判定双层增益:
  1. 替换了原单层中不相容接触(单层最差Er<-0.10 而双层三界面均>=-0.10)
  2. 改善最差界面反应驱动力(双层 E_worst > 两个单层的最差值)
  3. 增加必要电子阻隔(单层无阻隔而双层有)
  4. 明确功能分工(A 贴 Li 稳定 + B 贴 SSE 稳定, 单层各自不稳定)

输出: 05_单双层对照/bilayer_vs_monolayer.csv
只对通过组+严格档组合输出(控制规模), 判定规则逐列给出。
"""
import os, csv, json
from collections import Counter

ROOT = os.path.dirname(os.path.abspath(__file__))
NEW = os.path.dirname(ROOT)
TOP = os.path.dirname(NEW)
WS = os.path.join(TOP, "AB池双层筛选")
OUT = os.path.join(NEW, "05_单双层对照")
os.makedirs(OUT, exist_ok=True)

ER_TOL = -0.10
EG_CUT = 1.0


def f(x):
    try:
        if x in ("", None, "inf"):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def mono_metrics(x_id, sse, merged, er_li_map):
    """单层 X 直接夹 Li 与 SSE 之间"""
    r = merged.get(x_id)
    if r is None:
        return None
    er_li = f(r.get("Er_Li"))
    er_sse = f(r.get(f"Er_{sse}"))
    eg = f(r.get("band_gap_eV"))
    vred, vox = f(r.get("V_red")), f(r.get("V_ox"))
    worst = min([v for v in (er_li, er_sse) if v is not None], default=None)
    return {"Er_Li": er_li, "Er_SSE": er_sse, "Eg": eg,
            "worst_Er": worst, "V_red": vred, "V_ox": vox}


def judge(mA, mB, e_worst_bi, block_bi, gains):
    """返回 (判定, 说明)。注意判据优先级:
    - 判据2(双层最差界面优于任一单层)是实质性热力学改善, 优先报告;
    - 判据1 中 'Li|B|SSE 单层失败' 是结构性的(B 池未筛 Li|B 相容),
      单独标注为 'B层本不能贴Li' 而非双层独有增益;
    - 只有 'Li|A|SSE 单层失败(即 A 与 SSE 直接接触不相容)' 才是双层解决的真问题。"""
    if mA is None or mB is None or e_worst_bi is None:
        return "无法比较", "单层或双层数据缺失"
    a_mono_fail = mA["worst_Er"] is not None and mA["worst_Er"] < ER_TOL   # Li|A|SSE 失败 = A|SSE 不相容
    b_li_fail = mB["Er_Li"] is not None and mB["Er_Li"] < ER_TOL           # B 不能贴 Li (结构性)
    b_sse_fail = mB["Er_SSE"] is not None and mB["Er_SSE"] < ER_TOL        # B 连 SSE 都不相容(异常)

    mono_worst_best = max([m["worst_Er"] for m in (mA, mB) if m["worst_Er"] is not None],
                          default=None)
    gain2 = mono_worst_best is not None and e_worst_bi > mono_worst_best + 1e-9
    block_mono = (mA["Eg"] is not None and mA["Eg"] >= EG_CUT) and \
                 (mB["Eg"] is not None and mB["Eg"] >= EG_CUT)
    gain3 = (not block_mono) and block_bi

    if gain2:
        return "双层增益明确", f"改善最差界面驱动力(单层最好 {mono_worst_best:.3f} → 双层 {e_worst_bi:.3f})"
    if a_mono_fail and e_worst_bi >= ER_TOL:
        return "双层增益明确", "替换不相容接触(Li|A|SSE中A|SSE不相容, B层隔离后三界面达标)"
    if gain3:
        return "双层增益明确", "增加电子阻隔(单层均<1eV, 双层至少一层>=1eV)"
    # 功能分工: A 贴 Li 稳定但贴 SSE 不稳; B 贴 SSE 稳定但贴 Li 不稳
    if mA["Er_Li"] is not None and mA["Er_SSE"] is not None and \
       mB["Er_Li"] is not None and mB["Er_SSE"] is not None:
        if mA["Er_Li"] >= ER_TOL and mA["Er_SSE"] < ER_TOL and \
           mB["Er_SSE"] >= ER_TOL and mB["Er_Li"] < ER_TOL and e_worst_bi >= ER_TOL:
            return "双层增益明确", "功能分工: A只稳Li侧, B只稳SSE侧, 双层互补"
    if b_li_fail and not a_mono_fail and not b_sse_fail and e_worst_bi >= ER_TOL:
        return "可能有增益但证据不足", "B层本不能直接贴Li(结构性), A层承担Li侧接触; 双层相对Li|A|SSE单层无额外热力学改善"
    if gains.get("证据不足"):
        return "可能有增益但证据不足", gains["证据不足"]
    return "暂无已证明增益", "单层已兼容或双层未改善最差指标"


def run(sse, tag):
    merged = {r["mp_id"]: r for r in csv.DictReader(
        open(os.path.join(WS, "03_AB池", "母池_合并热力学数据.csv"), encoding="utf-8-sig"))}
    rows = []
    with open(os.path.join(WS, "05_结果", f"组合筛选_{sse}.csv"), encoding="utf-8-sig") as fo:
        for r in csv.DictReader(fo):
            if r["group"] != "通过组":
                continue
            aid, bid = r["A"], r["B"]
            mA = mono_metrics(aid, sse, merged, None)
            mB = mono_metrics(bid, sse, merged, None)
            e_worst_bi = f(r["Er_Li_A"]) if r["Er_Li_A"] != "" else None
            vals = [f(r["Er_Li_A"]), f(r["Er_A_B"]), f(r["Er_B_SSE"])]
            e_worst_bi = min([v for v in vals if v is not None], default=None)
            eg_a, eg_b = f(r["Eg_A"]), f(r["Eg_B"])
            block_bi = (eg_a is not None and eg_a >= EG_CUT) or (eg_b is not None and eg_b >= EG_CUT)
            gains = {}
            if mA is None or mB is None:
                gains["证据不足"] = "单层热力学数据缺失"
            verdict, note = judge(mA, mB, e_worst_bi, block_bi, gains)
            rows.append({
                "SSE": sse, "A": aid, "A_formula": r["A_formula"], "B": bid, "B_formula": r["B_formula"],
                "双层E_worst": e_worst_bi,
                "单层Li|A|SSE最差Er": mA["worst_Er"] if mA else None,
                "单层Li|B|SSE最差Er": mB["worst_Er"] if mB else None,
                "A单层Eg": mA["Eg"] if mA else None, "B单层Eg": mB["Eg"] if mB else None,
                "双层电子阻隔": int(block_bi),
                "A单层可传输": "", "B单层可传输": "",
                "判定": verdict, "判定说明": note,
            })
    with open(os.path.join(OUT, f"bilayer_vs_monolayer_{tag}.csv"), "w", newline="",
              encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"[{tag}] 对照 {len(rows)} 组: {dict(Counter(r['判定'] for r in rows))}")
    return rows


if __name__ == "__main__":
    run("Li2ZrCl6", "LZC")
    run("Li3YCl6", "LYC")

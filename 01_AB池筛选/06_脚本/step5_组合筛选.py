# -*- coding: utf-8 -*-
"""
规则3+5+6+7: A×B 组合筛选与三组输出

对每个 SSE(Li2ZrCl6_beta / Li3YCl6)独立评估 A池×B池 全部组合:

规则3(电子阻隔, 或): Eg_A >= 1.0 eV 或 Eg_B >= 1.0 eV。
  - 绝缘/绝缘、导体/绝缘、绝缘/导体 均保留(绝缘/导体标记层序注意: 导体在 SSE 侧);
  - 导体/导体(Eg_A<1 且 Eg_B<1)在电子阻隔路线上淘汰。
  - A/B 对调算不同结构, 需重新检查(本管线 A/B 角色由池定义, 对调即另一组合)。
  - "通过带隙"= 该层具备阻电子潜力, 不代表全堆栈电子绝缘。

规则4(三界面, 且): min[Er(Li,A), Er(A,B), Er(B,SSE)] >= -0.10 eV/atom。
  输出统一称"通过相容性容忍阈值"。Er(A,B) 来自 step4 矩阵(封闭两相混合);
  Er(Li,A) 封闭 PD(=0V 巨势); Er(B,SSE) 封闭 PD。方法逐条记录。

规则5(窗口匹配, 相邻层): W_X=[V_red^X, V_ox^X];
  ΔV_XY = min(V_ox^X,V_ox^Y) - max(V_red^X,V_red^Y);
  要求 ΔV_AB > 0 且 ΔV_B,SSE > 0(静态热力学保护路线);
  层序连接: 存在 v1 ∈ W_A∩W_B, v2 ∈ W_B∩W_SE 且 0 <= v1 <= v2。
  不要求三重交叠, 不要求 A-SSE 直接交叠。
  窗口缺失 / 仅边界相接(ΔV<=0.05 V) / 交叠过小(<0.1 V) → 待核;
  无交叠组合不进本路线通过清单, 可保留为动力学保护线索(kinetic_lead 标记)。

规则6(锂传输, 两层皆需):
  - 离子晶体: BVSE E_3D <= 1.0 eV(双层候选阈值; softBV 已算, BVEL13k 查表);
  - 金属/合金(is_intermetallic=1): 不得直接套用离子晶体 BVSE;
    有适用扩散计算或可信文献证据才可入通过组, 否则待核;
  - 方法不适用/数据缺失 → 待核, 不得把缺失当零势垒, 也不得因此淘汰。

规则7(三组输出):
  - 通过组: 全部条件满足且数据完整;
  - 待核组: 关键数据缺失 / 窗口边界敏感 / 传输方法未验证(合金);
  - 未通过组: 记录具体未通过条件。
  通过组内两层同为 BVSE 方法时按 E_3D,max = max(E_3D^A, E_3D^B) 升序排列,
  再比较三界面反应能、窗口交叠余量; 含合金组合单独排序。

输出: 05_结果/组合筛选_<sse>.csv (group 列: 通过组/待核组/未通过组/动力学线索)
"""
import os, csv, json

ROOT = os.path.dirname(os.path.abspath(__file__))
WS = os.path.dirname(ROOT)
TOP = os.path.dirname(WS)

EG_CUT = 1.0
ER_TOL = -0.10
DV_MARGIN = 0.10      # 交叠余量 < 0.1 V → 待核(边界敏感)
DV_EDGE = 0.05        # 仅边界相接阈值
E3D_MAX = 1.0
SSE_WIN = {"Li2ZrCl6": None, "Li3YCl6": None}  # 由 real_sse 数据填入


def f(x):
    try:
        if x in ("", None, "inf"):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def main(sse="Li2ZrCl6"):
    # SSE 窗口: 读 step2b 结果(SSE 自身 VASP+Cl修正条目的巨势窗口)
    swin = json.load(open(os.path.join(WS, "02_热力学数据", "SSE窗口.json"), encoding="utf-8"))
    sse_vred, sse_vox = swin[sse]["V_red"], swin[sse]["V_ox"]
    merged = {r["mp_id"]: r for r in csv.DictReader(open(os.path.join(WS, "03_AB池", "母池_合并热力学数据.csv"), encoding="utf-8-sig"))}
    a_pool = [r["mp_id"] for r in csv.DictReader(open(os.path.join(WS, "03_AB池", "A池_靠Li.csv"), encoding="utf-8-sig")) if r["A池标记"] == "入池"]
    b_pool = [r["mp_id"] for r in csv.DictReader(open(os.path.join(WS, "03_AB池", f"B池_靠SSE_{sse}.csv"), encoding="utf-8-sig")) if r["B池标记"] == "入池"]

    er_ab = {}
    ckpt = os.path.join(WS, "04_组合筛选", f"AB反应能_checkpoint_{sse}.jsonl")
    for line in open(ckpt, encoding="utf-8"):
        d = json.loads(line)
        er_ab[(d["A"], d["B"])] = d

    # BVSE 势垒
    bv = {}
    bvpath = os.path.join(TOP, "softBV_Jeon2025", "06_结果汇总", "Jeon2025_1610项_Li迁移势垒汇总.csv")
    for r in csv.DictReader(open(bvpath, encoding="utf-8-sig")):
        bv[r["mp_id"]] = r

    out = []
    for aid in a_pool:
        A = merged[aid]
        eg_a = f(A["band_gap_eV"])
        vred_a, vox_a = f(A["V_red"]), f(A["V_ox"])
        er_li_a = f(A["Er_Li"])
        e3d_a = None
        bva = bv.get(aid)
        if bva and bva.get("E_3D") not in ("", "inf", None):
            e3d_a = float(bva["E_3D"])
        a_inter = A.get("is_intermetallic", "0") in ("1", "True", 1)

        for bid in b_pool:
            if aid == bid:
                continue
            B = merged[bid]
            eg_b = f(B["band_gap_eV"])
            vred_b, vox_b = f(B["V_red"]), f(B["V_ox"])
            er_b_sse = f(B[f"Er_{sse}"])
            e3d_b = None
            bvb = bv.get(bid)
            if bvb and bvb.get("E_3D") not in ("", "inf", None):
                e3d_b = float(bvb["E_3D"])
            b_inter = B.get("is_intermetallic", "0") in ("1", "True", 1)

            dab = er_ab.get((aid, bid), {})
            er_ab_v = dab.get("Er_AB")
            reasons_fail, reasons_pend = [], []

            # 规则3
            block = (eg_a is not None and eg_a >= EG_CUT) or (eg_b is not None and eg_b >= EG_CUT)
            if not block:
                reasons_fail.append("规则3:导体/导体组合,电子阻隔路线淘汰")
            order_flag = ""
            if block and eg_a is not None and eg_b is not None:
                if eg_a >= EG_CUT and eg_b < EG_CUT:
                    order_flag = "层序注意:B层(SSE侧)为低带隙,电子阻隔靠A层(Li侧)"

            # 规则4
            if er_ab_v is None:
                reasons_pend.append(f"规则4:Er(A,B)缺失({dab.get('note','未算')})")
            elif min(er_li_a if er_li_a is not None else -9, er_ab_v,
                     er_b_sse if er_b_sse is not None else -9) < ER_TOL:
                reasons_fail.append("规则4:三界面反应能低于-0.10 eV/atom")
            if er_li_a is None or er_b_sse is None:
                reasons_pend.append("规则4:Li|A或B|SSE反应能缺失")

            # 规则5
            kinetic_lead = False
            if None in (vred_a, vox_a, vred_b, vox_b):
                reasons_pend.append("规则5:窗口缺失")
            else:
                dv_ab = min(vox_a, vox_b) - max(vred_a, vred_b)
                dv_bs = min(vox_b, sse_vox) - max(vred_b, sse_vred)
                if dv_ab <= 0 or dv_bs <= 0:
                    kinetic_lead = True  # 无交叠: 不进静态路线, 保留动力学线索
                else:
                    if dv_ab <= DV_EDGE or dv_bs <= DV_EDGE:
                        reasons_pend.append("规则5:窗口仅边界相接")
                    elif dv_ab < DV_MARGIN or dv_bs < DV_MARGIN:
                        reasons_pend.append("规则5:窗口交叠余量<0.1V")
                    # 层序连接: 存在 v1 ∈ W_A∩W_B, v2 ∈ W_B∩W_SE, 0<=v1<=v2
                    # ΔV_AB>0 且 ΔV_B,SE>0 已保证两交叠区间非空;
                    # 可选 v1 = min(0 以上, W_A∩W_B 下界), 只需 v2 区间上界 >= v1
                    v1_lo = max(vred_a, vred_b, 0.0)
                    v2_hi = min(vox_b, sse_vox)
                    if v2_hi < v1_lo:
                        reasons_pend.append("规则5:层序电位连接不满足(0<=v1<=v2)")

            # 规则6
            for tag, inter, e3d in [("A", a_inter, e3d_a), ("B", b_inter, e3d_b)]:
                if inter:
                    reasons_pend.append(f"规则6:{tag}层为金属间化合物,BVSE不适用,需合金扩散证据")
                elif e3d is None:
                    reasons_pend.append(f"规则6:{tag}层E_3D缺失")
                elif e3d > E3D_MAX:
                    reasons_fail.append(f"规则6:{tag}层E_3D={e3d:.2f}>1.0eV")

            if reasons_fail:
                grp = "未通过组"
            elif kinetic_lead:
                grp = "动力学线索"
            elif reasons_pend:
                grp = "待核组"
            else:
                grp = "通过组"

            e3d_max = max([x for x in (e3d_a, e3d_b) if x is not None], default=None)
            dv_ab_s = round(min(vox_a, vox_b) - max(vred_a, vred_b), 2) if None not in (vred_a, vox_a, vred_b, vox_b) else ""
            dv_bs_s = round(min(vox_b, sse_vox) - max(vred_b, sse_vred), 2) if None not in (vred_b, vox_b) else ""
            out.append({
                "SSE": sse, "A": aid, "A_formula": A["formula"], "B": bid, "B_formula": B["formula"],
                "Eg_A": eg_a, "Eg_B": eg_b, "层序标记": order_flag,
                "Er_Li_A": er_li_a, "Er_A_B": er_ab_v, "Er_B_SSE": er_b_sse,
                "V_red_A": vred_a, "V_ox_A": vox_a, "V_red_B": vred_b, "V_ox_B": vox_b,
                "dV_AB": dv_ab_s, "dV_B_SSE": dv_bs_s,
                "E3D_A": e3d_a, "E3D_B": e3d_b, "E3D_max": e3d_max,
                "A_is_intermetallic": int(a_inter), "B_is_intermetallic": int(b_inter),
                "group": grp, "未通过条件": ";".join(reasons_fail), "待核原因": ";".join(reasons_pend),
                "Er_AB_method": dab.get("method", ""), "Er_AB_note": dab.get("note", ""),
            })

    # 排序: 通过组按 E3D_max 升序(无 E3D_max 排后), 含合金组合单独排
    def sort_key(r):
        alloy = r["A_is_intermetallic"] or r["B_is_intermetallic"]
        return (r["group"] != "通过组", alloy, r["E3D_max"] if r["E3D_max"] is not None else 9e9,
                -(r["Er_A_B"] if r["Er_A_B"] is not None else 0))
    out.sort(key=sort_key)

    os.makedirs(os.path.join(WS, "05_结果"), exist_ok=True)
    fields = list(out[0].keys()) if out else []
    with open(os.path.join(WS, "05_结果", f"组合筛选_{sse}.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        wcsv = csv.DictWriter(fo, fieldnames=fields)
        wcsv.writeheader(); wcsv.writerows(out)
    from collections import Counter
    print(f"[{sse}] 组合总数 {len(out)} 分组: {dict(Counter(r['group'] for r in out))}")


if __name__ == "__main__":
    import sys
    sse = sys.argv[1] if len(sys.argv) > 1 else "Li2ZrCl6"
    main(sse)

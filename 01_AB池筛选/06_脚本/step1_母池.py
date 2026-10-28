# -*- coding: utf-8 -*-
"""
规则1: 候选母池(取消统一带隙门槛)

新规则:
  - 元素限制: 沿用现有毒性/放射性/成本规则(与 multi_screen_Jeon2025 screen1 一致)
  - 材料自身稳定性: Ehull <= 5 meV/atom(本项目范围设定)
  - 必须含 Li(当前路线保留; 纯金属前驱体暂不纳入)
  - 带隙: 母池不筛除低带隙材料(绝缘含锂化合物与含锂金属间化合物同时保留)
  - 氧化还原窗口: 保存上下限, 入池不统一要求宽度 >= 1 V(移到组合阶段判断)

输入: 开源数据集/A_Honrao2021_全体含锂池_含provenance_v2.csv (21,712 含锂体系)
输出: ../01_母池/母池_无带隙门槛.csv + ../07_日志/step1.log
"""
import csv, os, time

ROOT = os.path.dirname(os.path.abspath(__file__))
WS = os.path.dirname(ROOT)                      # AB池双层筛选/
TOP = os.path.dirname(WS)                       # 苏州实验室课题学习记录/
SRC = os.path.join(TOP, "开源数据集")

EXCLUDE_NOBLE = {"Pt", "Au", "Pd", "Ir", "Rh", "Ru", "Os", "Ag"}
EXCLUDE_RADIO = {"Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Pm", "Tc", "Ra"}
EXCLUDE_TOXIC = {"Hg", "Pb", "Cd", "As", "Tl", "Be"}
EXCLUDE_RARE = {"Re"}
EXCLUDE_INERT = {"Xe", "Ar", "Kr", "Ne", "He", "Rn"}
# 阴离子(用于"金属间化合物"分类标记, 不用于筛除)
ANIONS = {"O", "S", "F", "Cl", "Br", "I", "N", "Se", "Te", "H", "C", "P", "B"}
TM = {"Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
      "Y","Zr","Nb","Mo","Ru","Rh","Pd","Ag","Cd",
      "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"}

EHULL_MAX = 0.005  # eV/atom, 本项目范围设定

def main():
    rows = list(csv.DictReader(open(os.path.join(SRC, "A_Honrao2021_全体含锂池_含provenance_v2.csv"), encoding="utf-8-sig")))
    out = []
    stats = {"输入含锂体系": len(rows), "Ehull>5meV剔除": 0, "元素规则剔除": 0, "入池": 0,
             "入池_低带隙Eg<0.5": 0, "入池_金属间化合物": 0}
    for r in rows:
        eh = float(r["Ehull_eV_per_atom"])
        if eh > EHULL_MAX:
            stats["Ehull>5meV剔除"] += 1
            continue
        els = set(r["chemsys"].split("-"))
        reason = []
        for tag, group in (("贵金属", EXCLUDE_NOBLE), ("放射性", EXCLUDE_RADIO),
                           ("高毒", EXCLUDE_TOXIC), ("稀有Re", EXCLUDE_RARE), ("惰性气", EXCLUDE_INERT)):
            hit = els & group
            if hit:
                reason.append(tag + ":" + ",".join(sorted(hit)))
        if reason:
            stats["元素规则剔除"] += 1
            continue
        eg = float(r["band_gap_eV"]) if r["band_gap_eV"] not in ("", None) else 0.0
        intermetallic = not (els & ANIONS)
        r["exclude_reason"] = ""
        r["is_intermetallic"] = int(intermetallic)
        r["has_TM"] = int(bool(els & (TM - {"Zn"})))
        r["Eg_cut_1.0"] = int(eg >= 1.0)
        out.append(r)
        stats["入池"] += 1
        stats["入池_低带隙Eg<0.5"] += int(eg < 0.5)
        stats["入池_金属间化合物"] += int(intermetallic)

    fields = list(out[0].keys())
    os.makedirs(os.path.join(WS, "01_母池"), exist_ok=True)
    with open(os.path.join(WS, "01_母池", "母池_无带隙门槛.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(out)

    os.makedirs(os.path.join(WS, "07_日志"), exist_ok=True)
    with open(os.path.join(WS, "07_日志", "step1.log"), "w", encoding="utf-8") as f:
        f.write("规则1 母池统计(取消统一带隙门槛)\n")
        for k, v in stats.items():
            f.write(f"{k}: {v}\n")
    print(stats)

if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"耗时 {time.time()-t0:.1f}s")

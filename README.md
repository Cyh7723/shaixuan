# 负极双层保护膜高通量筛选

对 Li | A | B | SSE 堆栈的双层负极保护膜进行多级高通量计算筛选，目标 SSE 为 Li₂ZrCl₆ (LZC)、Li₃YCl₆ (LYC) 两种卤化物固态电解质（Li₃InCl₆ 作为扩展探索）。

**主线叙事与全流程总结见 [`00_主线_61材料/主线工作思路_61材料双层筛选.md`](00_主线_61材料/主线工作思路_61材料双层筛选.md)**——从 21,712 含锂体系三级收缩到 196 对组合 / 61 个材料（NEB 计算与实验验证的输入清单）。

## 目录结构（主线 / 探索分区）

```
00_主线_61材料/          ★ 主线交付物与全流程总结
    ├── 主线工作思路_61材料双层筛选.md   （主叙事：三级筛选→196对/61材料→NEB/实验）
    ├── pairs_196.csv                   （196 对组合明细）
    └── materials_61.csv                （61 材料清单）

01_AB池筛选/            主线·第一轮：七规则双层初筛 (1,860 → 6,501/8,074 → 81对短名单)
02_二次筛选/            主线·第二轮：严格档复筛 (→302 → 196对/61材料) + 单层对照 + 力学/反应产物
    ├── 07_排名候选/双层严格档_196对全性质表.csv  （每对全性质，含GVRH剪切模量信息列）
    ├── 07_排名候选/双层严格档与单层对照_说明.md
    └── 07_排名候选/单层严格档_*.csv            （单层对照：双层必要性证据）
99_探索_支线/           支线探索：Li3InCl6路线、通道几何、阈值敏感性、Pareto替代排名
```

## 筛选流程（主线）

```
母池 1,860 项 (Honrao2021 含锂池, Ehull≤5 meV/atom, 无带隙门槛)
    │
第一轮 (01_AB池筛选/): 七规则双层初筛
    ├─ A池 (贴Li): 340 项,  Er(Li,A) ≥ -0.10 eV/atom
    ├─ B池 (贴SSE): LZC 793 / LYC 1070 项
    └─ 组合筛选 → 通过组 LZC 6,501 / LYC 8,074 → 双SSE短名单 81 对
    │
第二轮 (02_二次筛选/): 严格档复筛（主线核心步骤）
    ├─ 严格档 (E_worst≥-0.03, E_max≤0.6): 302 对
    └─ ★ 主线交付: E_worst≥-0.02 & E_max≤0.6 & Vred_A≤0.1V → 196 对 / 61 材料
    │
单层对照（同口径）: LZC 2 种 / LYC 0 种 / 双SSE交集 0 种 —— 双层必要性的量化证据
    │
下游: 61 材料可合成性分层 (P0 实验梯队 ~25 种) → VASP 体相优化 → CI-NEB 势垒验证 (进行中)
```

## 核心方法

| 性质 | 方法 | 来源 |
|---|---|---|
| 电化学窗口 | 巨势相图 (GrandPotentialPhaseDiagram) | pymatgen |
| 界面反应能 | 封闭两相混合 (InterfacialReactivity) | pymatgen |
| 离子迁移势垒 | BVSE (BVEL13k 查表 + softBV 自算) | 已有数据库 |
| 弹性张量 | MatTen 等变图神经网络 | wengroup/matten |
| 电子阻隔 | PBE 带隙 ≥1.0 eV (至少一层) | MP |
| SSE 能量 | VASP 自算 + MP2020 Cl 修正 (−0.614 eV/Cl) | halide_vasp |

## 主要结果

### 主线：196 对 / 61 材料（严格档，2026-11 定稿）

- **61 材料构成**：A 池 17 种（LiCl/LiBr/LiI/LiF + Li₂O/Li₂S/Li₂Te/Li₃N + LiREO₂×5 + 氮化物 3）贴 Li；B 池 46 种（LiCl/LiBr + 硫酸盐 + 氟化复盐 + Ga/In 卤化物等）贴 SSE。
- 排名前列：LiI|LiCl、LiBr|LiCl、Li₂O|LiCl（双 SSE 交集前三，E_max≤0.31 eV，B 层高度收敛于 LiCl）。
- 每对全性质（含 GVRH 剪切模量信息列）见 `02_二次筛选/07_排名候选/双层严格档_196对全性质表.csv`。

### 同口径单层对照（双层必要性的直接证据）

| 口径 | 单层通过 | 同口径双层 |
|---|---|---|
| LZC | **2 种**（LiCl、LiBr） | 196 对 |
| LYC | **0 种** | — |
| 双 SSE 交集 | **0 种** | — |

单层近乎归零而双层近两百对通过：把"对 Li 稳定"与"对 SSE 稳定"拆给两层是解空间存在的关键。从 1,860 全母池复跑结论一致（139 池严格包含全部严格档解；16 项 no-entry 缺口中仅 4 项 Eg≥1.0，不改变定性结论）。

## 复现

```bash
# 第一轮（七规则）
cd 01_AB池筛选/06_脚本
python step1_母池.py && python step2_热力学数据.py
python step3_AB池.py && python step4_AB界面矩阵.py Li2ZrCl6
python step5_组合筛选.py Li2ZrCl6 && python step6_排名评分.py

# 第二轮（严格档+对照）
cd 02_二次筛选/08_脚本
python step00_审计.py && python stepA_阈值敏感性.py
python stepB_matten.py && python stepCD_反应产物体积.py
python stepE_单双层对照.py && python stepG_pareto排名.py
```

依赖: `pymatgen>=2024`, `pandas`, `numpy`, `torch>=2.0`, `e3nn`, `matten` (GitHub wengroup/matten)

## 大文件说明

以下数据文件超过 GitHub 单文件限制，未包含在本仓库，可通过上述脚本复现：
- 完整组合筛选表（LZC 27 万行 / LYC 36 万行）
- 三界面反应产物明细 (interface_reactions.jsonl, 60 MB)
- MatTen 预训练模型权重 (41 MB, 可从 https://github.com/wengroup/matten 获取)
- 弹性张量矩阵 (elastic_tensors.json)

## 引用

- MatTen: Wen et al., Digital Discovery, 2024, DOI: 10.1039/D3DD00233K
- 双层设计先例: Wang, Aoyagi & Mueller, JMCA, 2021, DOI: 10.1039/D1TA04987A
- 枝晶力学: Monroe & Newman, JES, 2005, DOI: 10.1149/1.1850854
- 电子泄漏: Han et al., Nature Energy, 2019, DOI: 10.1038/s41560-018-0312-z

## 局限声明

- BVSE 势垒为 D 级证据，不换算室温电导率（主线 NEB 验证正是为此）
- MatTen 为机器学习预测非 DFT，非正定条目已标记未截断
- GVRH 剪切模量仅作信息列，不作普适淘汰条件（Monroe-Newman 模型针对聚合物/锂界面）
- 体积变化为体相加和近似，不等于真实界面应力
- −0.02 eV/atom 严格档为动力学容忍口径，非严格热力学稳定

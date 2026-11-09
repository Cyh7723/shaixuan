# 负极双层保护膜高通量筛选

对 Li | A | B | SSE 堆栈的双层负极保护膜进行两轮高通量计算筛选，目标 SSE 为 Li₂ZrCl₆ (LZC)、Li₃YCl₆ (LYC)、Li₃InCl₆ (LIC) 三种卤化物固态电解质。

## 筛选流程

```
母池 1,860 项 (Honrao2021 含锂池, Ehull≤5 meV/atom)
    │
第一轮 (AB池双层筛选/): 七规则初筛
    ├─ A池 (贴Li): 340 项,  Er(Li,A) ≥ -0.10 eV/atom
    ├─ B池 (贴SSE): LZC 793 / LYC 1070 / LIC 937 项
    └─ 组合筛选 → 通过组 LZC 6,501 / LYC 8,074 / LIC 6,775
    │
第二轮 (二次筛选/): 新性质 + 严格档 + Pareto
    ├─ 严格档 (E_worst≥-0.03, E_max≤0.6): 302 / 465 / 587
    ├─ Pareto 前沿: 27 / 25 / 32
    ├─ 全严格档 (+GVRH_min≥35 + 无导电产物 + |dV|≤20%): 22 / 85 / 186
    └─ 三SSE全严格档交集: 15 对化学式
```

## 目录结构

- `01_AB池筛选/` — 第一轮：七规则初筛管线
  - `06_脚本/` — step1~step6 完整可复现脚本
  - `README.md` — 第一轮方法详解与已知局限
- `02_二次筛选/` — 第二轮：新性质筛选
  - `08_脚本/` — step00 审计 ~ stepH Li₃InCl₆ + stepI 后处理
  - `execution_report.md` — 第二轮执行报告
  - `00_审计/` ~ `11_Li3InCl6路线/` — 各工作包结果

## 核心方法

| 性质 | 方法 | 来源 |
|---|---|---|
| 电化学窗口 | 巨势相图 (GrandPotentialPhaseDiagram) | pymatgen |
| 界面反应能 | 封闭两相混合 (InterfacialReactivity) | pymatgen |
| 离子迁移势垒 | BVSE (BVEL13k 查表 + softBV 自算) | 已有数据库 |
| 弹性张量 | MatTen 等变图神经网络 | wengroup/matten |
| 电子阻隔 | PBE 带隙 ≥1.0 eV (至少一层) | MP |
| 排序 | 多目标 Pareto 前沿 (E_max/E_worst/dV_min/GVRH_min) | 本项目 |

## 主要结果

### 三种 SSE 全严格档共同候选（15 对）

A 层只有 LiF 和 Li₂S 两种材料能满足全部条件：

| A | B (代表) |
|---|---|
| LiF | Li₃AlF₆, LiTi₂(PO₄)₃, LiZr₂(PO₄)₃, LiAl(SiO₃)₂, Li₂Mn₃CoO₈ 等 12 种 |
| Li₂S | Li₃AlF₆, LiAl(SiO₃)₂, LiAl(Si₂O₅)₂ |

### Pareto 前沿三者交集（5 对）

Li₂O|Li₂SO₄, Li₂S|Li₃AlF₆, LiCl|LiBF₄, LiErO₂|Li₂SO₄, LiHoO₂|Li₂SO₄

### 双层严格档 196 对全性质表（2026-11-09 补充）

81 对双 SSE 短名单经严格复筛（Vred_A≤0.1V、三界面 Er≥−0.02 eV/atom、Emax≤0.6 eV）得 196 对、覆盖 61 种材料，作为 Li 迁移 NEB 计算的输入清单。每对组合的全部性质（含 MatTen 剪切模量 GVRH，仅作信息列不筛选）见 `02_二次筛选/07_排名候选/双层严格档_196对全性质表.csv`，方法与结果解读见同目录 `双层严格档与单层对照_说明.md`。

### 同口径单层对照（双层必要性的直接证据）

对 139 种单层候选用同一套严格标准（ErLi≥−0.02 & ErSSE≥−0.02 & E3D≤0.6 & Eg≥1.0 & Vred≤0.1）筛选：

| 口径 | 单层通过 | 同口径双层 |
|---|---|---|
| LZC | **2 种**（LiCl、LiBr） | 196 对 |
| LYC | **0 种**（LiCl/LiBr/LiI 卡在 −0.022~−0.027） | — |
| 双 SSE 交集 | **0 种** | — |

单层近乎归零而双层近两百对通过，量化印证了把"对 Li 稳定"与"对 SSE 稳定"拆给两层的必要性。数据见 `02_二次筛选/07_排名候选/单层严格档_*.csv`。

**全母池口径复核**：139 池是 Jeon2025 单层漏斗（Er≥−0.1 入池）的历史产物；从 1,860 母池直接跑同一套严格漏斗结果完全一致（LZC 152→2、LYC 0），因严格档阈值（−0.02）远严于 139 池入池条件，139 池严格包含全部可行解。唯一不确定来源是 16 项 2025.03 后新入 MP、无 S3 缓存能量的条目（其中 Eg≥1.0 仅 4 项），已记录于说明文档，不改变定性结论。见 `单层严格档_全母池1860_通过者_LZC_LYC.csv`。

## 复现

```bash
# 第一轮
cd 01_AB池筛选/06_脚本
python step1_母池.py && python step2_热力学数据.py
python step3_AB池.py && python step4_AB界面矩阵.py Li2ZrCl6
python step5_组合筛选.py Li2ZrCl6 && python step6_排名评分.py

# 第二轮
cd 02_二次筛选/08_脚本
python step00_审计.py && python stepA_阈值敏感性.py
python stepB_matten.py && python stepCD_反应产物体积.py
python stepE_单双层对照.py && python stepF_通道传输电位.py
python stepG_pareto排名.py
NPROC=14 python stepH_Li3InCl6.py && python stepI_Li3InCl6后处理.py
```

依赖: `pymatgen>=2024`, `pandas`, `numpy`, `torch>=2.0`, `e3nn`, `matten` (GitHub wengroup/matten)

## 大文件说明

以下数据文件超过 GitHub 单文件限制，未包含在本仓库，可通过上述脚本复现：
- 完整组合筛选表（LZC 27 万行 / LYC 36 万行 / LIC 32 万行）
- 三界面反应产物明细 (interface_reactions.jsonl, 60 MB)
- MatTen 预训练模型权重 (41 MB, 可从 https://github.com/wengroup/matten 获取)
- 弹性张量矩阵 (elastic_tensors.json)

## 引用

- MatTen: Wen et al., Digital Discovery, 2024, DOI: 10.1039/D3DD00233K
- 双层设计先例: Wang, Aoyagi & Mueller, JMCA, 2021, DOI: 10.1039/D1TA04987A
- 枝晶力学: Monroe & Newman, JES, 2005, DOI: 10.1149/1.1850854
- 电子泄漏: Han et al., Nature Energy, 2019, DOI: 10.1038/s41560-018-0312-z

## 局限声明

- BVSE 势垒为 D 级证据，不换算室温电导率
- MatTen 为机器学习预测非 DFT，非正定条目已标记未截断
- GVRH≥35 GPa 仅作高模量标签，不作普适淘汰条件（Monroe-Newman 模型针对聚合物/锂界面）
- 体积变化为体相加和近似，不等于真实界面应力

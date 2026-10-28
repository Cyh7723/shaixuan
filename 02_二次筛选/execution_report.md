# 二次筛选执行报告

日期: 2026-10-21
体系: Li | A | B | Li2ZrCl6 / Li3YCl6
输入: `AB池双层筛选` 第一轮管线 (母池 1860, A池 340, B池 LZC 793 / LYC 1070)

---

## 1. 基线复现 (阶段0)

全部 18 项基线计数与报告一致 (`00_审计/baseline_reproduction.csv`):

| 项目 | LZC | LYC |
|---|---|---|
| 组合总数 | 269,602 | 363,715 |
| 通过组 | 6,501 | 8,074 |
| 待核组 | 8,854 | 15,114 |
| 动力学线索 | 6,983 | 9,911 |
| 未通过组 | 247,264 | 330,616 |

**数量差异解释** (`00_审计/data_audit.md`):
- LZC: 340×793=269,620, step5 排除 A==B 同 mp_id 配对 18 项 → 269,602。
- LYC: 340×1070=363,800, 排除 85 项 → 363,715。
- 排除规则: 同一 mp_id 不能同时作为 A 层和 B 层(同一晶体结构不构成双层);
  化学式相同但晶型不同不合并, 全程按 mp_id 区分。
- 类别互斥: group 为单值字段, 判定顺序 未通过>动力学线索>待核>通过, 无重复统计。

## 2. 工作包A: 阈值敏感性

`01_阈值敏感性/threshold_sensitivity_{LZC,LYC}.csv`: 3×4×2×3=72 组阈值交叉统计。

关键阈值组合(配对数):

| Er下限 | E3D上限 | LZC | LYC |
|---|---|---|---|
| -0.10 | 1.0 (基线) | 6,501 | 8,074 |
| -0.05 | 0.7 | 812 | 1,162 |
| -0.03 | 0.6 (严格档) | **302** | **465** |
| -0.03 | 0.5 | 152 | 238 |

- 严格档两 SSE 共同化学式对 247 对。
- 严格档中双绝缘(Eg_A>=1 且 Eg_B>=1)占 274/302; 全部组合至少一层 Eg>=2。
- 逐级淘汰见 `attrition_by_stage_{LZC,LYC}.csv`。

## 3. 工作包B: MatTen 弹性张量

- 母池 1,860 个唯一晶体结构全部完成 MatTen 预训练模型预测 (CPU, 约 6 分钟)。
- 张量存于 `02_MatTen力学/elastic_tensors.json` (6x6 Voigt, GPa, pymatgen IEEE 约定)。
- 派生指标 `material_elasticity.csv`: KVRH/GVRH/E/泊松比/AU。
- 质检: 正常 1,842; 带标记 18 (非正定/病态/泊松比超界/AU异常); 失败 0。
- GVRH 分布: 中位 41.2 GPa, min -8.6, max 410.4。
- 本地修补(不改数值): torch.load weights_only 兼容; Windows 长路径(>260字符)目录名缩短。
  详见 `matten_run_report.md`。

## 4. 工作包C+D: 反应产物与体积变化

- 唯一界面 35,884 个 (Li|A 340, A|B 33,681, B|SSE 1,863)。
- 成功 32,130; chemspace>6 跳过 3,754 (合法限制, 相图维度爆炸)。
- `03_反应产物/interface_reactions.jsonl`: 每界面全部 kink 的配平反应/反应能/产物系数/dV。
- `03_反应产物/reaction_product_properties.csv`: 5,387 种产物的带隙/体积/金属性。
- `04_体积变化/reaction_volume_changes.csv`: 每界面最负 Er 处 dV 与全 kink 最大 |dV|。
- 开放 Li 储库: Li|A 界面同时输出 "计入金属Li体积" 与 "不计入" 两列。

## 5. 工作包E: 单层 vs 双层

`05_单双层对照/bilayer_vs_monolayer_{LZC,LYC}.csv` (通过组 6,501/8,074 全部对照):

| 判定 | LZC | LYC |
|---|---|---|
| 双层增益明确 | 4,698 | 5,499 |
| 可能有增益但证据不足 | 1,721 | 2,407 |
| 暂无已证明增益 | 82 | 168 |

增益类型(LZC): 改善最差界面驱动力 4,226; 替换不相容接触 364; 增加电子阻隔 108。
判定时区分了 "B层本不能贴Li(结构性)" 与 "A|SSE 不相容(双层解决的真问题)",
避免把池定义导致的必然失败算作双层独有增益。

## 6. 工作包F: 通道/传输/电位

`06_通道几何/`:
- `channel_geometry.csv`: 严格档 104 个唯一材料的 BVSE E1D/E2D/E3D 与贯通维度推断。
- `transport_evidence.csv`: 证据等级全部为 D(仅 BVSE); 无 A/B/C 级 σ, ASR 不生成伪精确值。
- `potential_sensitivity.csv`: V∈[-0.5,7] 巨势分解扫描。
  - 96/104 材料在还原侧(V<0.5)分解析出金属 Li(导电), 这是卤化物/硫族化合物的普遍热力学行为;
  - 9 种材料窗口扫到边界(LiF/LiCl/LiBr/LiI/Li2O/Li2S/Li3N 等还原侧截断于 -0.5 V), 已标记 truncated;
  - CAVD 几何瓶颈半径未算(需单独安装 cavd 包), 列为后续工作。

## 7. 阶段3: Pareto 排序与代表候选

Pareto 字段(仅物理含义可比): E_max(↓), E_worst(↑), dV_min(↑), GVRH_min(↑)。
GVRH_min = min(G_A, G_B) 作保守双层模量指标, 不用 max。

| | LZC | LYC |
|---|---|---|
| 严格档 | 302 | 465 |
| 完整数据(GVRH 无缺失) | 302 | 465 |
| **Pareto 前沿** | **27** | **25** |
| 分路线代表(每路线≤5, 同A≤2) | 9 | 10 |

双 SSE 共同 Pareto 前沿化学式对 10 对:
Li2O|Li2SO4, Li2O|LiBr, Li2O|LiCl, Li2S|Li3AlF6, LiCl|LiBF4,
LiErO2|Li2SO4, LiF|Li3AlF6, LiF|LiAl(SiO3)2, LiHoO2|Li2SO4, LiI|LiCl

新旧排名对比 `old_vs_new_{LZC,LYC}.csv`: 旧等权评分前三名(LiI|LiCl, LiCl|LiBF4, Li2O|LiCl)
在新 Pareto 前沿中分别列第 2/1/5(LZC), 一致性良好; 新排名额外考虑了模量与导电产物风险。

## 8. 从 6501/8074 到两位数的路径

```
通过组(基线)              6501 / 8074
  └─ 严格档(Er>=-0.03, E3D<=0.6)    302 / 465
       └─ Pareto 前沿(E_max/E_worst/dV/GVRH_min)   27 / 25
            └─ 双 SSE 共同前沿化学式对              10
            └─ 分路线代表候选(去重后)                ~15
```

## 9. 主要局限与诚实声明

1. **反应能为封闭两相混合** (InterfacialReactivity, norm=True), 归一化 eV/atom;
   Li|A 与 0 V 巨势等价, 但 A|B 与 B|SSE 不是开放体系。
2. **BVSE 势垒系统性偏高** (~0.2-0.3 eV), E3D<=0.6 为宽档; 不换算室温电导率。
3. **MatTen 为机器学习预测**, 非 DFT; 18 个非正定条目仅标记, 未判定真实不稳定。
4. **体积变化为体相加和近似**, 不等于真实界面应力; 仅用于风险排序。
5. **电位敏感性还原侧导电产物** 标注了 Li 金属析出, 但未包含动力学钝化;
   96/104 材料在 V<0.5 分解是热力学预期, 不代表实际界面一定失效。
6. **CAVD 几何瓶颈** 未算(需 cavd 包), 通道维度仅由 BVSE E1D/E2D/E3D 推断。
7. **GVRH>=35 GPa** 仅作高模量标签, 不作淘汰条件; 25-45 GPa 灰区无置信含义。

## 10. 复现

```bash
cd 二次筛选_20261021/08_脚本
python step00_审计.py          # 基线复现(秒级)
python stepA_阈值敏感性.py     # 阈值统计(分钟级)
python stepB_matten.py         # MatTen 批量(~6 min CPU)
python stepCD_反应产物体积.py  # 反应产物+体积(多进程 ~17 min, NPROC=14)
python stepE_单双层对照.py     # 单层对照(分钟级)
python stepF_通道传输电位.py   # 通道/传输/电位(分钟级)
python stepG_pareto排名.py     # Pareto+代表候选(秒级)
```

所有脚本只读 `AB池双层筛选/` 与 `开源数据集/`, 不覆盖旧数据。

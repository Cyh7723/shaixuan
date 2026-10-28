# A/B 双层负极涂层初筛管线 —— 规则与运行说明

日期: 2026-10-16。框架: 用户指定的"组合能否完成保护功能"七规则初筛, 堆栈 Li|A|B|SSE,
目标 SSE = Li₂ZrCl₆(β, 在壳) 与 Li₃YCl₆(亚稳 22.7 meV/atom), 均为 VASP 自算 + Cl 校正。
本阶段**不引入**缺陷浓度、显式界面、膜厚优化; 原位反应前驱体/聚合物/多相复合另案处理。

## 核心改动(相对旧 1610 池管线)

取消母池统一绝缘(Eg≥0.5)与统一窗口宽度(≥1 V)两个硬门槛;
电子阻隔移到组合阶段(规则3, 或条件); 电化学窗口移到相邻界面匹配(规则5);
三界面相容性(规则4)与两层锂传输(规则6)保留。

## 七规则实现与阈值

### 规则1 母池(step1_母池.py)
- 输入: Honrao2021 全体含锂池 21,712 体系(S3 Open Data, 能量已含 MP2020 校正)。
- 门槛: Ehull ≤ 5 meV/atom(项目范围); 含 Li(纯金属前驱体暂排除);
  元素规则(贵金属 Pt/Au/Pd/Ir/Rh/Ru/Os/Ag; 放射/锕系 Ac/Th/Pa/U/Np/Pu/Am/Pm/Tc/Ra;
  高毒 Hg/Pb/Cd/As/Tl/Be; Re; 惰性气体)。
- **不设带隙门槛**: 同时保留绝缘含锂化合物与含锂金属间化合物。
- 记录 redox 窗口上下限(不强制宽度 ≥1 V)。
- 结果: 1,860 项(含低带隙 446, 金属间化合物 270)。
- 输出: `01_母池/母池_无带隙门槛.csv`

### 热力学数据补齐(step2_热力学数据.py, step2b_SSE窗口.py)
- 电化学窗口 V_red/V_ox: 巨势相图(GrandPotentialPhaseDiagram), 开放 Li,
  μ_Li = μ0 − V, 电位窗 [−0.5, 7.0] V, 0.25 V 粗扫 + 边界二分(7 次);
  判据 e_above_hull ≤ 自身Ehull + 50 meV(Honrao2021 亚稳余量)。
  归一化: 相图能量 per atom, MP S3 manifest 原值(**不二次校正**)。
- Er(Li,X): 封闭组成相图两相混合(InterfacialReactivity, norm=True,
  use_hull_energy=False), eV/atom(混合原子)。封闭 vs Li 等价 0 V 巨势界面反应,
  统一记 `closed_PD_vs_Li(=GPIR@0V)`。1,610 项旧结果复用, 余者新算。
- Er(X,SSE): 同法对 Li₂ZrCl₆_beta / Li₃YCl₆; 化学空间 >6 元跳过(记 chemspace>6 → 待核)。
- SSE 自身窗口(step2b): Li₂ZrCl₆ V_red=1.25 V, V_ox=4.46 V(在壳);
  Li₃YCl₆ V_red=0.75 V, V_ox=4.37 V(Ehull 22.7 meV)。
- 输出: `02_热力学数据/窗口_checkpoint.jsonl`(1,860) `反应能_checkpoint.jsonl`(1,860) `SSE窗口.json`

### 规则2 A/B 池(step3_AB池.py)
- A 池(靠 Li): Er(Li,A) ≥ −0.10 eV/atom **必需**; V_red^A ≤ 0.3 V 暂时保留
  (**0.3 V 是对部分还原驱动力的经验容忍度, 不得写成"0 V 严格稳定"**)。
  电子绝缘推迟到组合阶段。入池 340, 待核 21。
- B 池(靠 SSE): Er(B,SSE) ≥ −0.10 eV/atom 必需。电子绝缘同样推迟。
  入池 793(Li₂ZrCl₆) / 1,070(Li₃YCl₆), 待核 295/298(多为 chemspace>6 或窗口缺失)。
- 输出: `03_AB池/A池_靠Li.csv`, `B池_靠SSE_Li2ZrCl6.csv`, `B池_靠SSE_Li3YCl6.csv`,
  `母池_合并热力学数据.csv`

### 规则4 三界面相容性(step4_AB界面矩阵.py + step5)
min[Er(Li,A), Er(A,B), Er(B,SSE)] ≥ −0.10 eV/atom, 输出统一称"通过相容性容忍阈值"。
- Er(A,B): 封闭组成相图两相混合。实现: 按并集元素集合分桶建相图,
  沿 A-B 混合线在 11 个原子分数 x∈[0.05,0.95] 采样
  Er(x) = E_hull((1−x)A + xB) − [(1−x)E_A + xE_B], 取最小。
  与 InterfacialReactivity 最低 kink 对照 10 组最大偏差 0.002 eV/atom;
  Er=0 即 A、B 两相平衡(同一棱)。**不以"所有界面都在 0 V 计算"作为唯一依据**:
  A|B 为无锂库封闭混合, Li|A 另有封闭 vs Li(=0V 巨势), B|SSE 为封闭 vs SSE。
- 方法逐条记录于 checkpoint 的 `method` 字段。
- 输出: `04_组合筛选/AB反应能_checkpoint_Li2ZrCl6.jsonl` / `_Li3YCl6.jsonl`

### 规则3 电子阻隔(step5, 或条件)
E_g^A ≥ 1.0 eV **或** E_g^B ≥ 1.0 eV。
- 绝缘/绝缘、导体/绝缘、绝缘/导体 均保留(绝缘/导体标记层序注意: 导体在 SSE 侧);
- 导体/导体(Eg_A<1 且 Eg_B<1)在电子阻隔路线淘汰;
- A/B 对调算不同结构, 需重新检查(本管线 A/B 角色由池定义, 对调即另一组合);
- 通过带隙 = 该层具备阻电子潜力, 不代表全堆栈电子绝缘。

### 规则5 窗口匹配(step5, 相邻层)
W_X = [V_red^X, V_ox^X]; ΔV_XY = min(V_ox^X,V_ox^Y) − max(V_red^X,V_red^Y)。
- 静态热力学保护路线要求 ΔV_AB > 0 且 ΔV_B,SSE > 0;
- 层序连接: 存在 v1 ∈ W_A∩W_B, v2 ∈ W_B∩W_SE, 0 ≤ v1 ≤ v2
  (实现: v1 下界 = max(V_red^A, V_red^B, 0), 需 min(V_ox^B, V_ox^SSE) ≥ v1 下界);
- **不要求**三重交叠, **不要求** A-SSE 直接交叠;
- 窗口缺失 / 仅边界相接(ΔV ≤ 0.05 V) / 交叠 < 0.1 V → 待核;
- 无交叠组合不进本路线通过清单, 保留为动力学保护线索(`动力学线索` 组)。

### 规则6 锂传输(step5, 两层皆需)
- 离子晶体: BVSE/softBV E_3D ≤ 1.0 eV(双层候选阈值);
- 金属/合金(is_intermetallic): **不直接套用**离子晶体 BVSE,
  需适用扩散计算或可信文献证据, 否则待核;
- 方法不适用/数据缺失 → 待核, 不把缺失当零势垒也不因此淘汰。

### 规则7 三组输出(step5)
- 通过组: 全部条件满足且数据完整;
- 待核组: 关键数据缺失 / 窗口边界敏感 / 传输方法未验证(合金);
- 未通过组: 记录具体未通过条件;
- 排序: 通过组内两层同为 BVSE 方法按 E_3D,max = max(E_3D^A, E_3D^B) 升序,
  再比较三界面反应能、窗口交叠余量; 含合金组合单独排序。
- 输出: `05_结果/组合筛选_Li2ZrCl6.csv`, `组合筛选_Li3YCl6.csv`

### 通过组多性质排名(step6_排名评分.py)

方法: Aykol et al., Nat. Commun. 2016(电池正极涂层高通量筛选)的
**Weighted sum** 一路——各指标 min–max 归一化到 0–1, 统一"越大越好",
等权合成总分 S = (f_E + f_R + f_V)/3, 总分越高越优。

指标(通过组数据最完整的三项):

| 指标 | 定义 | 方向 |
|---|---|---|
| E_max | max(E3D_A, E3D_B) | 越低越优 |
| E_worst | min(Er_Li_A, Er_A_B, Er_B_SSE) | 越接近 0 越优 |
| dV_min | min(dV_AB, dV_B_SSE) | 越大越优 |

权重敏感性: 在 (E_max, E_worst, dV_min) 权重单形上以 0.25 步长采样 15 组
正权重, 记录每个组合名次的最好/最差区间(rank_min_w/rank_max_w)。

结果(等权基线):
- 前 30 名中权重扰动下始终留在前 30 的: LZC 15/30, LYC 14/30;
- 双 SSE 交集(去重化学式对 4,993 对)按两路线名次最大值升序,
  前三名 LiI|LiCl、LiBr|LiCl、Li₂O|LiCl 在两个 SSE 上均在前 3,
  且权重扰动下名次区间不越出前 3 —— 为最稳健的双层候选。
- 输出: `05_结果/排名评分_Li2ZrCl6.csv`, `排名评分_Li3YCl6.csv`, `双SSE交集排名.csv`

## 已知局限

- PBE 带隙系统性低估, Eg 阈值 1.0 eV 留了余量但仍可能误判边界材料;
- BVSE/softBV 势垒相对 NEB 系统性偏高 ~0.2–0.3 eV, E_3D ≤ 1.0 为宽档;
- −0.10 eV/atom 为动力学容忍余量, 非严格热力学稳定;
- 11 点混合线扫描对极强非对称最低 kink 可能有 ≤ 数 meV 残差(对照验证 ≤2 meV);
- 窗口扫描上限 7.0 V, 超高氧化电位材料 V_ox 记 7.0 封顶;
- 排名评分只做等权基线 + 权重敏感性; E_max/E_worst/dV_min 三项相关性低
  (|ρ| ≤ 0.16), 指标独立性好;
- 本阶段不含缺陷化学、显式界面、成膜/厚度、动力学(界面迁移)评价——留终验(VASP/实验)。

## 筛选结果(2026-10-16 完成)

| 分组 | Li₂ZrCl₆ 路线 | Li₃YCl₆ 路线 |
|---|---|---|
| 组合总数 | 269,602 (340×793) | 363,715 (340×1,070) |
| **通过组** | **6,501** | **8,074** |
| 待核组 | 8,854 | 15,114 |
| 动力学线索 | 6,983 | 9,911 |
| 未通过组 | 247,264 | 330,616 |

- **两个 SSE 均通过的 A|B 对: 6,495(去重化学式对 4,993)**, 已导出
  `05_结果/双SSE均通过_去重化学式对.csv`。
- 通过组全部为非合金组合(合金组合因规则6 BVSE 不适用一律进待核组)。
- 未通过主因: 三界面反应能低于 −0.10 eV/atom(LZC 19.0 万 / LYC 23.4 万对),
  B 层 E_3D > 1.0 eV(LZC 12.6 万 / LYC 16.1 万), 导体/导体组合(LZC 3.0 万 / LYC 5.8 万)。
- 待核主因: A 层为金属间化合物(BVSE 不适用)、B/A 层 E_3D 缺失、窗口边界敏感。
- E_3D,max 最低(0.264 eV, LiCl/LiI 的 BVSE 值)的代表双通过组合:
  LiCl|LiBF₄、LiCl|LiI、LiI|LiGaTe₂、LiCl|Na₃Li₃Sc₂F₁₂、LiCl|LiZnPS₄ 等
  (Er_AB=0 两相平衡或弱放热, 窗口交叠 0.6–3.8 V)。
- 层序注意(导体在 SSE 侧, 电子阻隔靠 A 层): LZC 366 对 / LYC 375 对, 已标记。

### 已知组合验证(与旧 63 体系结论一致)
LiCl|LiF、LiCl|LiBr、LiCl|LiI 均落入通过组, Er_AB=0(两相平衡),
dV_AB 3.4–4.8 V, 与旧报告"全卤素双层最优"结论吻合。

## 优先级分档与短名单(step6_优先级分档.py, 试运行阈值)

不覆盖宽筛通过组, 只在其上加优先档(阈值尚未认定为最佳, 可调整):

| 指标 | 宽筛通过条件 | 优先档 | 理由 |
|---|---|---|---|
| E_3D,max | ≤1.0 eV | ≤0.6 eV | 传输潜力优先 |
| E_worst = min(Er_Li,A, Er_A,B, Er_B,SSE) | ≥−0.10 eV/atom | ≥−0.03 eV/atom | 减少对反应钝化的依赖 |
| ΔV_min = min(ΔV_AB, ΔV_B,SSE) | >0 且避开边界敏感 | ≥0.2 V | 窗口匹配余量 |
| 电子阻隔 | 至少一层 Eg≥1.0 | 保持原要求 | 不重新要求两层都绝缘 |
| Ehull | ≤5 meV/atom | 保持原要求 | 已很严格 |

综合分(0–100, 权重为试运行值): 0.45×传输 + 0.35×反应能 + 0.20×窗口。

分档结果:

| | Li₂ZrCl₆ | Li₃YCl₆ |
|---|---|---|
| 通过组 | 6,501 | 8,074 |
| 优先档 | 302 | 464 |
| 优先档去重化学式对 | 263 | 397 |
| 两 SSE 均入优先档 | — | 247 对 |

短名单(两 SSE 均入优先档 + 每个 A 化学式限额 5): **81 对, 覆盖 18 种 A 材料**
(`05_结果/短名单_双SSE_每A限额5.csv`)。A 层覆盖 LiCl/LiBr/LiI/LiF 卤素、
Li₂O/Li₂S/Li₂Te 硫族、Li₃N 及 Li₃BN₂/Li₃ScN₂/Li₃YN₂ 氮化物、
LiHoO₂/LiErO₂/LiYO₂/LiTmO₂/LiLuO₂ 稀土氧化物、Li₅AlO₄/Li₈TiO₆;
B 层以 LiCl/LiBr/LiI + Li₂SO₄/Li₂SeO₄ + LiBF₄/Li₃AlF₆/LiAlCl₄ 等为主。

排名前列(综合分, LZC/LYC):
1. LiCl|LiBF₄ (80.2/53.7) —— E3D=0.264, E_worst 0/−0.023
2. LiI|LiCl (80.2/53.7) —— 层序相反的另一结构, 分别保留
3. Li₂Te|LiCl (78.8/52.3)
4. Li₂O|LiCl (77.3/50.8)
5. LiBr|LiCl (76.6/50.1)

注意 LiCl/LiBr/LiI 互搭(A 与 B 对调)在表中各占一行——层序不同属不同结构。
优先档全量(含晶型备选): `05_结果/优先档_Li2ZrCl6.csv`, `优先档_Li3YCl6.csv`。

## 复现

```bash
cd AB池双层筛选/06_脚本
python step1_母池.py          # 母池(秒级)
python step2_热力学数据.py    # 窗口+对Li/SSE反应能(~40 min, 旧结果自动复用)
python step2b_SSE窗口.py      # SSE 自身窗口(分钟级)
python step3_AB池.py          # A/B 池(秒级)
python step4_AB界面矩阵.py Li2ZrCl6   # A×B 矩阵(~2h, 可中断重跑)
python step4_AB界面矩阵.py Li3YCl6
python step5_组合筛选.py Li2ZrCl6
python step5_组合筛选.py Li3YCl6
python step6_优先级分档.py     # 优先档 + 排名 + 短名单(秒级)
```

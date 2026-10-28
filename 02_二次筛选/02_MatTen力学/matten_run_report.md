# MatTen 运行报告

日期: 2026-10-21

## 版本与环境
- 模型: wengroup/matten GitHub main (克隆于 2026-10-21), model_identifier=20230627, model_final.ckpt
- torch 2.14.1+cpu, e3nn 0.5.x, pymatgen 2026.5.4, CPU 运行
- 单位: GPa; Voigt 约定: pymatgen ElasticTensor (IEEE, 工程剪切)
- 对官方代码的两处本地修补(均不改变模型权重与数值):
  1. torch.load weights_only 兼容 (torch>=2.6 默认改为 True, 需加安全白名单);
  2. processed 缓存目录名过长(含方括号与等号, 超过 260 字符)在 Windows 触发 PyTorchFileWriter error 3 与 fsspec glob 字符错误, 改为短目录名 proc_*。

## 结果
- 母池 1860 个唯一结构全部完成预测, 张量存于 elastic_tensors.json (6x6 Voigt, GPa)。
- 质检: 正常 1842; 带标记 18 (非正定/病态/泊松比超界/AU异常); 失败 0。
- GVRH 分布: min -8.6 / 中位 41.2 / max 410.4 GPa。
- 负值与非正定条目仅标记, 未截断, 未据此判定真实材料不稳定(预测可能失真)。

## 使用约束
- MatTen 原论文测试集剪切模量 MAE 约 8.38 GPa, 不能作为单样本置信区间。
- 25-45 GPa 待核区间为工程分档, 无概率置信含义。
- GVRH>=35 GPa 仅作"高模量路线"标签, 不作淘汰条件; 不用 max(GA,GB) 代表双层有效模量。

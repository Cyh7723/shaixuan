# 复现说明

## 环境

- Python 3.11
- pymatgen 2026.x, torch 2.x CPU, e3nn, matten (GitHub main)
- Windows 10/11 (MatTen 需两处 Windows 兼容补丁，见下)

## 安装

```bash
pip install -r requirements.txt
pip install git+https://github.com/wengroup/matten.git
```

## MatTen Windows 兼容补丁

matten 官方代码在 Windows 有两个问题（不影响 Linux/Mac）：

1. **torch.load weights_only**：torch≥2.6 默认 `weights_only=True`，需加白名单。
   修改 `matten/predict.py` 的 `get_pretrained_model`：
   ```python
   torch.serialization.add_safe_globals([TensorRegressionTask])
   ```

2. **缓存目录名过长**：processed 目录名含 `[]` 和 `=`，超过 Windows 260 字符路径限制。
   修改 `matten/dataset/structure_scalar_tensor.py`：
   ```python
   processed_dirname = f"proc_{self.tensor_target_name[:16]}_{self.tensor_target_format[:8]}"
   ```

## 数据依赖

第一轮和第二轮脚本需要以下外部数据（未包含在本仓库，需自行准备）：

| 数据 | 来源 | 用途 |
|---|---|---|
| Honrao2021 全体含锂池 (21,712 项) | S3 Open Data | 母池 |
| MP GGA raw entries (153,167) | Materials Project | 相图/反应能 |
| Li₂ZrCl₆ / Li₃YCl₆ VASP 校正条目 | 自算 | 目标 SSE |
| BVEL13k BVSE 数据库 | LiTraj | 迁移势垒查表 |
| 晶体结构 CIF (21,712 个) | MP | MatTen 输入 |

准备好后按 README 的复现命令顺序执行即可。

## 已知问题

- pymatgen 的 `coord_cython.pyd` 在某些 Windows 企业策略下会被 Code Integrity 间歇性阻止，重试通常可恢复
- 反应产物计算 (stepCD) 需 14 进程约 17 分钟；Li₃InCl₆ AB 矩阵约 52 分钟

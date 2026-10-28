# -*- coding: utf-8 -*-
"""
工作包B: MatTen 批量弹性张量预测 (母池 1860 唯一结构)

- 使用官方预训练模型 (wengroup/matten, model_identifier=20230627), 不训练。
- 小批量先试 8 个结构验证接口, 然后分批 200 处理全池, checkpoint 逐批落盘。
- 派生指标: KV/KR/KVRH, GV/GR/GVRH, E, nu, AU (GPa)。
- 质检: 张量有限值/对称性/正定性(刚度矩阵特征值)/异常模量与泊松比标记。
  不正定不截断为正值, 只标记。

输出(02_MatTen力学/):
  elastic_tensors.json     mp_id -> 6x6 Voigt 张量 (GPa, pymatgen IEEE 约定)
  material_elasticity.csv  派生指标 + 质检标记
  matten_run_report.md
"""
import os, sys, csv, json, time, tempfile, warnings
warnings.filterwarnings("ignore")

# MatTen 的 PyG InMemoryDataset 会在当前工作目录写 processed/ 缓存;
# 项目路径含中文导致 PyTorchFileWriter 打开失败, 切到纯 ASCII 临时目录运行。
_ASCII_CWD = tempfile.mkdtemp(prefix="matten_")
os.chdir(_ASCII_CWD)

# 单结构独立预测封装: matten.predict 对整批共享一个 processed 缓存目录,
# 批内任一结构数据损坏会让整批失败; 逐结构调用并在失败时回退到 None。
_MODEL = None
_MODEL_CFG = None

def _init_model():
    """加载一次模型并驻留内存, 避免每个结构重复读权重。"""
    global _MODEL, _MODEL_CFG
    if _MODEL is None:
        from matten.predict import get_pretrained_model, get_pretrained_config
        _MODEL = get_pretrained_model(identifier="20230627", checkpoint="model_final.ckpt")
        _MODEL.eval()
        _MODEL_CFG = get_pretrained_config("20230627")
    return _MODEL, _MODEL_CFG

def predict_one(structure):
    """单结构预测 -> pymatgen ElasticTensor; 失败抛异常由上层记录。"""
    import torch
    from pymatgen.analysis.elasticity import ElasticTensor
    from matten.predict import get_data_loader, evaluate
    model, cfg = _init_model()
    loader = get_data_loader([structure], "20230627", batch_size=1)
    preds = evaluate(model, loader,
                     tensor_target_name=cfg["data"]["tensor_target_name"],
                     tensor_target_formula=cfg["data"]["tensor_target_formula"])
    if not preds:
        raise RuntimeError("empty_predictions")
    return ElasticTensor(preds[0])

ROOT = os.path.dirname(os.path.abspath(__file__))
NEW = os.path.dirname(ROOT)
TOP = os.path.dirname(NEW)
WS = os.path.join(TOP, "AB池双层筛选")
OUT = os.path.join(NEW, "02_MatTen力学")
CKPT_T = os.path.join(OUT, "elastic_tensors.json")
os.makedirs(OUT, exist_ok=True)

CIF_DIR = os.path.join(TOP, "开源数据集", "structures_all_Li")

import numpy as np


def derive(C):
    """C: 6x6 Voigt (GPa) -> 派生指标与质检"""
    from pymatgen.analysis.elasticity import ElasticTensor
    flags = []
    if not np.all(np.isfinite(C)):
        return None, ["非有限值"]
    asym = float(np.abs(C - C.T).max())
    if asym > 1e-6 * max(1.0, float(np.abs(C).max())):
        flags.append(f"非对称({asym:.2e})")
    et = ElasticTensor.from_voigt(C)
    ev = np.linalg.eigvalsh(C)
    if ev.min() <= 0:
        flags.append(f"非正定(minEig={ev.min():.2f})")
    try:
        kv, kr = float(et.k_voigt), float(et.k_reuss)
        gv, gr = float(et.g_voigt), float(et.g_reuss)
    except Exception as ex:
        return None, flags + [f"求逆失败:{type(ex).__name__}"]
    k, g = (kv + kr) / 2, (gv + gr) / 2
    E = 9 * k * g / (3 * k + g) if (3 * k + g) != 0 else float("nan")
    nu = (3 * k - 2 * g) / (2 * (3 * k + g)) if (3 * k + g) != 0 else float("nan")
    au = 5 * gv / gr + kv / kr - 6 if (gr != 0 and kr != 0) else float("nan")
    if gr < 0 or kr < 0:
        flags.append("Reuss模量为负(病态)")
    if abs(nu) > 0.5 + 1e-6:
        flags.append(f"泊松比超界({nu:.2f})")
    if au < -0.5 or au > 100:
        flags.append(f"AU异常({au:.1f})")
    return {"KV": kv, "KR": kr, "KVRH": k, "GV": gv, "GR": gr, "GVRH": g,
            "E_VRH": E, "poisson": nu, "AU": au}, flags


def main():
    from pymatgen.core import Structure
    _init_model()   # 预热模型
    import torch, e3nn

    pool = list(csv.DictReader(open(os.path.join(WS, "01_母池", "母池_无带隙门槛.csv"),
                                    encoding="utf-8-sig")))
    ids = [r["mp_id"] for r in pool]

    tensors = {}
    if os.path.exists(CKPT_T):
        tensors = json.load(open(CKPT_T, encoding="utf-8"))
    todo = [m for m in ids if m not in tensors]
    if not todo:
        # 全部张量已算完, 只做派生指标
        pass
    print(f"[MatTen] 总数 {len(ids)}, 已完成 {len(tensors)}, 待算 {len(todo)}", flush=True)

    # 记录版本
    ver = {"matten_repo": "https://github.com/wengroup/matten (main, cloned 2026-10-21)",
           "model_identifier": "20230627", "checkpoint": "model_final.ckpt",
           "torch": torch.__version__, "e3nn": e3nn.__version__,
           "device": "cuda" if torch.cuda.is_available() else "cpu",
           "unit": "GPa", "voigt_convention": "pymatgen ElasticTensor (IEEE, 工程剪切)"}

    t0 = time.time()
    B = 200
    n_err = 0
    for i in range(0, len(todo), B):
        batch_ids = todo[i:i + B]
        structs, ok_ids = [], []
        for m in batch_ids:
            p = os.path.join(CIF_DIR, f"{m}.cif")
            if not os.path.exists(p):
                tensors[m] = {"error": "cif_missing"}
                n_err += 1
                continue
            try:
                s = Structure.from_file(p)
                if not s.is_ordered:
                    tensors[m] = {"error": "disordered"}
                    n_err += 1
                    continue
                structs.append(s); ok_ids.append(m)
            except Exception as ex:
                tensors[m] = {"error": f"parse:{type(ex).__name__}"}
                n_err += 1
        if structs:
            # 批量预测: 同一批共享一个 processed 目录(短名), 批内损坏结构由 matten
            # 的 failed_entries 机制返回 None; 批级异常再逐结构重试定位。
            try:
                from matten.predict import predict
                ets = predict(structs, batch_size=min(64, len(structs)))
                if not isinstance(ets, list):
                    ets = [ets]
                for m, et in zip(ok_ids, ets):
                    if et is None:
                        tensors[m] = {"error": "predict_failed_entry"}
                        n_err += 1
                        continue
                    C = np.array(et.voigt, dtype=float)
                    tensors[m] = {"C_voigt_GPa": np.round(C, 4).tolist()}
            except Exception as ex:
                # 批级失败 → 逐结构重试
                for m, s in zip(ok_ids, structs):
                    try:
                        et = predict_one(s)
                        C = np.array(et.voigt, dtype=float)
                        tensors[m] = {"C_voigt_GPa": np.round(C, 4).tolist()}
                    except Exception as ex2:
                        tensors[m] = {"error": f"predict:{type(ex2).__name__}:{str(ex2)[:60]}"}
                        n_err += 1
        json.dump(tensors, open(CKPT_T, "w", encoding="utf-8"))
        el = time.time() - t0
        print(f"[MatTen] {min(i + B, len(todo))}/{len(todo)} 批 累计错误 {n_err} 用时 {el:.0f}s", flush=True)

    # 派生指标
    rows = []
    for m in ids:
        rec = tensors.get(m, {})
        if "C_voigt_GPa" not in rec:
            rows.append({"mp_id": m, "状态": rec.get("error", "未算"), "KVRH": "", "GVRH": "",
                         "E_VRH": "", "poisson": "", "AU": "", "质量标记": rec.get("error", "")})
            continue
        d, flags = derive(np.array(rec["C_voigt_GPa"]))
        if d is None:
            rows.append({"mp_id": m, "状态": "质检失败", "KVRH": "", "GVRH": "",
                         "E_VRH": "", "poisson": "", "AU": "", "质量标记": ";".join(flags)})
            continue
        rows.append({"mp_id": m, "状态": "ok" if not flags else "标记",
                     **{k: round(v, 3) for k, v in d.items()}, "质量标记": ";".join(flags)})
    with open(os.path.join(OUT, "material_elasticity.csv"), "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    n_ok = sum(1 for r in rows if r["状态"] == "ok")
    n_flag = sum(1 for r in rows if r["状态"] == "标记")
    n_bad = len(rows) - n_ok - n_flag
    with open(os.path.join(OUT, "matten_run_report.md"), "w", encoding="utf-8") as fo:
        fo.write(f"# MatTen 运行报告\n\n日期: 2026-10-21\n\n")
        fo.write("## 版本与环境\n")
        for k, v in ver.items():
            fo.write(f"- {k}: {v}\n")
        fo.write(f"\n## 结果\n- 总数 {len(rows)}; 正常 {n_ok}; 带质检标记 {n_flag}; 失败/异常 {n_bad}。\n")
        fo.write("- MatTen 原论文测试集剪切模量 MAE≈8.38 GPa, 不构成单样本置信区间。\n")
        fo.write("- 预测不正定的条目仅标记'待核', 未截断为正值, 未据此判定材料真实不稳定。\n")
    print(f"[MatTen] DONE ok={n_ok} flagged={n_flag} bad={n_bad}")


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""
规则4(数据, v5 分桶+混合线扫描): A×B 界面反应能矩阵 E_rxn(A,B)

每桶(并集元素集合 <=6 元)建一次 PhaseDiagram。对每对 (A,B):
  沿 A-B 混合线在 11 个原子分数 x∈[0.05..0.95] 采样,
  Er(x) = E_hull((1-x)A + xB) - [(1-x)E_A + xE_B], 取最小值。
  该值与 pymatgen InterfacialReactivity(norm=True, use_hull_energy=False)
  的最低 kink 一致(10 组对照最大偏差 0.002 eV/atom); Er=0 即 A、B 两相平衡不反应。

方法记录: 封闭组成相图两相混合反应, 归一化 eV/atom(混合原子)。
化学空间 >6 元跳过 -> chemspace>6(组合进待核)。

输出: 04_组合筛选/AB反应能_checkpoint_<sse>.jsonl
"""
import pickle, os, csv, json, time, warnings
from collections import defaultdict
warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.abspath(__file__))
WS = os.path.dirname(ROOT)
TOP = os.path.dirname(WS)

from pymatgen.core import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram

XS = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]


def main(sse="Li2ZrCl6"):
    raw = pickle.load(open(os.path.join(TOP, "multi_screen_Jeon2025", "02_数据获取", "pd_cache", "mp_gga_raw_entries.pkl"), "rb"))
    hal = pickle.load(open(os.path.join(TOP, "halide_vasp", "halide_entries_corrected.pkl"), "rb"))
    entries = raw + [hal["Li2ZrCl6_beta"], hal["Li3YCl6"]]
    by_formula = {}
    for e in entries:
        by_formula.setdefault(e.composition.reduced_formula, []).append(e)
    by_els = {}
    for e in entries:
        key = tuple(sorted(el.symbol for el in e.composition.elements))
        by_els.setdefault(key, []).append(e)

    a_rows = [r for r in csv.DictReader(open(os.path.join(WS, "03_AB池", "A池_靠Li.csv"), encoding="utf-8-sig"))
              if r["A池标记"] == "入池"]
    b_rows = [r for r in csv.DictReader(open(os.path.join(WS, "03_AB池", f"B池_靠SSE_{sse}.csv"), encoding="utf-8-sig"))
              if r["B池标记"] == "入池"]

    ckpt = os.path.join(WS, "04_组合筛选", f"AB反应能_checkpoint_{sse}.jsonl")
    done = set()
    if os.path.exists(ckpt):
        for line in open(ckpt, encoding="utf-8"):
            try:
                d = json.loads(line)
                done.add((d["A"], d["B"]))
            except Exception:
                pass

    def tgt(r):
        cands = by_formula.get(r["formula"])
        return min(cands, key=lambda e: e.energy_per_atom) if cands else None
    ta_map = {r["mp_id"]: tgt(r) for r in a_rows}
    tb_map = {r["mp_id"]: tgt(r) for r in b_rows}

    buckets = defaultdict(list)
    presolved = []
    for ra in a_rows:
        ta = ta_map[ra["mp_id"]]
        if ta is None:
            continue
        for rb in b_rows:
            if (ra["mp_id"], rb["mp_id"]) in done:
                continue
            tb = tb_map[rb["mp_id"]]
            if tb is None:
                continue
            if ra["mp_id"] == rb["mp_id"]:
                presolved.append((ra, rb, 0.0, "same_material"))
                continue
            els = {el.symbol for el in ta.composition.elements} | {el.symbol for el in tb.composition.elements}
            if len(els) > 6:
                presolved.append((ra, rb, None, "chemspace>6"))
                continue
            buckets[tuple(sorted(els))].append((ra, rb))

    total = sum(len(v) for v in buckets.values()) + len(presolved)
    print(f"[{sse}] A池 {len(a_rows)} × B池 {len(b_rows)} = {len(a_rows)*len(b_rows)} 对, "
          f"已完成 {len(done)}, 待算 {total} (预解 {len(presolved)}, 分桶 {len(buckets)})", flush=True)

    fo = open(ckpt, "a", encoding="utf-8")
    def emit(ra, rb, er, note, rxn=None):
        res = {"A": ra["mp_id"], "B": rb["mp_id"],
               "A_formula": ra["formula"], "B_formula": rb["formula"],
               "method": "closed_PD_mix(norm=eV/atom,mixing-line-11pt-scan)",
               "Er_AB": er, "note": note}
        if rxn is not None:
            res["rxn_AB"] = rxn
        fo.write(json.dumps(res, ensure_ascii=False) + "\n")

    for ra, rb, er, note in presolved:
        emit(ra, rb, er, note)
    fo.flush()

    t0 = time.time(); n = 0; n_zero = 0; n_neg = 0
    for els_key, pairs in buckets.items():
        s = set(els_key)
        rel = [e for k, v in by_els.items() if set(k) <= s for e in v]
        try:
            pd = PhaseDiagram(rel)
        except Exception:
            for ra, rb in pairs:
                emit(ra, rb, None, "error:PD_build")
                n += 1
            continue
        for ra, rb in pairs:
            ta = ta_map[ra["mp_id"]]; tb = tb_map[rb["mp_id"]]
            ca, cb = ta.composition, tb.composition
            na, nb = ca.num_atoms, cb.num_atoms
            da, db = ca.get_el_amt_dict(), cb.get_el_amt_dict()
            ea, eb = ta.energy_per_atom, tb.energy_per_atom
            try:
                best = 0.0
                comps = []
                for x in XS:
                    amt = {el: (1 - x) * da.get(el, 0) / na + x * db.get(el, 0) / nb for el in els_key}
                    comps.append((x, Composition(amt)))
                for x, comp in comps:
                    er = pd.get_hull_energy(comp) - ((1 - x) * ea + x * eb)
                    if er < best:
                        best = er
                if abs(best) < 1e-9:
                    best = 0.0
                if best == 0.0:
                    n_zero += 1
                    emit(ra, rb, 0.0, "tie_line_or_no_rxn")
                else:
                    n_neg += 1
                    emit(ra, rb, round(float(best), 4), "")
            except Exception as ex:
                emit(ra, rb, None, f"error:{type(ex).__name__}:{str(ex)[:60]}")
            n += 1
            if n % 10000 == 0:
                fo.flush()
                el = time.time() - t0
                print(f"[{sse}] {n}/{total} {el:.0f}s ({el/n*1000:.1f}ms/对) Er=0:{n_zero} 反应:{n_neg} ETA {(total-n)*el/n/60:.0f}min", flush=True)
        fo.flush()
    fo.close()
    print(f"[{sse}] DONE n={n} Er=0:{n_zero} 反应:{n_neg} {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else "Li2ZrCl6")

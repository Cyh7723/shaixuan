# -*- coding: utf-8 -*-
"""
目标 SSE 自身电化学窗口(规则5 需要 W_SSE)

方法同母池窗口: 巨势相图, 开放 Li, 0.25 V 粗扫 + 二分细化,
tol = 自身Ehull + 50 meV。SSE 条目为 VASP 自算 + Cl 修正(halide_entries_corrected.pkl),
相图其余条目为 MP S3 manifest(已含 MP2020 校正, 不二次校正)。
输出: 02_热力学数据/SSE窗口.json
"""
import pickle, os, json, warnings
warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.abspath(__file__))
WS = os.path.dirname(ROOT)
TOP = os.path.dirname(WS)

from pymatgen.core import Element
from pymatgen.analysis.phase_diagram import PhaseDiagram, GrandPotentialPhaseDiagram, GrandPotPDEntry


def main():
    raw = pickle.load(open(os.path.join(TOP, "multi_screen_Jeon2025", "02_数据获取", "pd_cache", "mp_gga_raw_entries.pkl"), "rb"))
    hal = pickle.load(open(os.path.join(TOP, "halide_vasp", "halide_entries_corrected.pkl"), "rb"))
    sses = {"Li2ZrCl6": hal["Li2ZrCl6_beta"], "Li3YCl6": hal["Li3YCl6"]}
    entries = raw + list(sses.values())
    mu0 = min((e for e in entries if e.composition.reduced_formula == "Li"),
              key=lambda e: e.energy_per_atom).energy_per_atom

    def stable_at(rel, target, mu, tol):
        gpd = GrandPotentialPhaseDiagram(rel, {Element("Li"): mu})
        ge = GrandPotPDEntry(target, {Element("Li"): mu})
        _, ea = gpd.get_decomp_and_e_above_hull(ge, allow_negative=True)
        return ea <= tol

    out = {}
    for name, target in sses.items():
        els = {el.symbol for el in target.composition.elements}
        rel = [e for e in raw if all(el.symbol in els | {"Li"} for el in e.composition.elements)]
        rel = rel + [target]
        pd = PhaseDiagram(rel)
        _, own = pd.get_decomp_and_e_above_hull(target)
        tol = own + 0.0501
        Vs = [round(-0.5 + i * 0.25, 3) for i in range(31)]
        flags = [stable_at(rel, target, mu0 - V, tol) for V in Vs]
        if not any(flags):
            out[name] = {"status": "never_stable"}; continue
        first = flags.index(True)
        last = len(Vs) - 1 - flags[::-1].index(True)

        def refine(lo, hi):
            for _ in range(7):
                mid = (lo + hi) / 2
                if stable_at(rel, target, mu0 - mid, tol):
                    lo = mid
                else:
                    hi = mid
            return (lo + hi) / 2

        vred = Vs[first] if first == 0 else refine(Vs[first - 1], Vs[first])
        vox = Vs[last] if last == len(Vs) - 1 else refine(Vs[last], Vs[last + 1])
        out[name] = {"V_red": round(vred, 2), "V_ox": round(vox, 2),
                     "Ehull_own_meV": round(own * 1000, 1),
                     "method": "GPD_openLi_mu=mu0-V_scan[-0.5,7.0]_tol=Ehull+50meV_SSE_VASP-Cl-corrected"}
        print(name, out[name])

    with open(os.path.join(WS, "02_热力学数据", "SSE窗口.json"), "w", encoding="utf-8") as fo:
        json.dump(out, fo, ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()

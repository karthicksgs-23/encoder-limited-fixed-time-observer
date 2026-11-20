"""Part A: the test as fixed in advance (criteria in arm_margin_test.py).
Part B: exploratory sensitivity sweep, added AFTER the pilot showed no violations.
Part B results must never be presented as the fixed test.
"""
import json
import multiprocessing as mp
import pathlib
import sys
import time

import numpy as np

import arm_margin_test as A

ENC12 = 2 * np.pi / 4096
D_FULL = None


def d_guaranteed():
    Jlo = A.J0 * (1 - A.SPREAD["J"])
    return (1 / Jlo - A.B_NOM) * A.U_MAX + (A.B0 * 1.3 * max(A.W_U, -A.W_L) * 1.5 + A.G0 * 1.2 + A.FC0 * 1.5 + A.D_SIN + A.PULSE_MAX) / Jlo + A.F_MAX


def job(j):
    A.RHO = j["rho"]
    sc = A.scenarios(j["runs"])
    t0 = time.time()
    s, _ = A.simulate(sc, A.margin_fn(j["kind"], j["mp"]), j["lam"], dt=j["dt"], quant=j["quant"], noise_std=j["noise"])
    r = A.summary(j["name"], s, j["runs"])
    r.update({k: j[k] for k in ("part", "lam", "rho", "quant", "noise", "dt")}, **j["mp"], fits=j["fits"], secs=time.time() - t0)
    print(j["part"], j["name"], "lam", j["lam"], "q", round(j["quant"], 5), "viol", r["violating_runs"], "excess", round(r["worst_excess"], 3),
          "rms", round(r["rms_median"], 4), "e2", round(r["max_e2"], 3), "breach", r["margin_breaches"], f'{r["secs"]:.0f}s', flush=True)
    return r


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


def main():
    D = d_guaranteed()
    room = min(A.W_U, -A.W_L)
    jobs = []
    margins = {}
    for lam in (10.0, 20.0, 40.0, 60.0):
        tt, env = A.error_envelope(lam, D)
        margins[lam] = A.margin_from_envelope(tt, env)
        print("margin", lam, margins[lam], flush=True)
    # Part A
    for name, kind in (("no margin", "none"), ("constant margin", "constant"), ("time-varying margin", "tv")):
        jobs.append(dict(part="A", name=name, kind=kind, lam=60.0, rho=0.05, quant=0.0, noise=0.0, dt=1e-4, runs=1000, mp=margins[60.0], fits=True))
    # Part B
    for lam in margins:
        fits = margins[lam]["m0"] < 0.9 * room and margins[lam]["m_ss"] < 0.9 * room
        for q in (0.0, ENC12):
            for name, kind in (("no margin", "none"), ("constant margin", "constant"), ("time-varying margin", "tv")):
                if kind != "none" and not fits:
                    continue
                jobs.append(dict(part="B", name=name, kind=kind, lam=lam, rho=0.01, quant=q, noise=0.0, dt=2e-4, runs=300, mp=margins[lam], fits=fits))
    # dt check: same setting at both step sizes
    jobs.append(dict(part="dtcheck", name="no margin", kind="none", lam=40.0, rho=0.01, quant=ENC12, noise=0.0, dt=1e-4, runs=300, mp=margins[40.0], fits=True))
    print(len(jobs), "jobs", flush=True)
    with mp.Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as p:
        rows = p.map(job, jobs, chunksize=1)
    out = dict(d_guaranteed=D, margins={str(k): v for k, v in margins.items()}, rows=rows)
    (A.OUT / "all_results.json").write_text(json.dumps(clean(out), indent=1))
    print("saved", flush=True)


if __name__ == "__main__":
    main()

"""Test 4: check predictions of the derived gain rule. Written before any of these runs.

Rule:  lam_c = sqrt( b0 * u_max / (kappa * q)^a2 ),  kappa = 1/2 + sigma/q
where q is the encoder count, sigma the noise standard deviation and a2 the injection exponent.
kappa is a heuristic for the peak output error: half a count from quantization plus one noise
standard deviation.

Base design, 12-bit, 1 kHz, 3 seeds of 100 runs. Ceiling = lowest gain, on a grid with 8% steps,
with median saturated share above 10% (median over seeds).

Pass rule, fixed now: every measured ceiling within 15% of its prediction. The rule ignores the
torque needed for the motion itself, so somewhat lower measured ceilings would not be surprising.
"""
import json
import multiprocessing as mp
import time

import numpy as np

import arm_margin_test as A
import test3_realistic as T
from run_all import clean

Q12 = 2 * np.pi / 4096
CASES = [dict(group="noise", sigma=s, a2=0.5, b0=2.0) for s in (0.0, 0.25, 0.5, 1.0, 2.0)]
CASES += [dict(group="exponent", sigma=0.5, a2=a, b0=2.0) for a in (0.25, 0.75)]
CASES += [dict(group="b0", sigma=0.5, a2=0.5, b0=b) for b in (1.5, 3.0)]


def predict(c):
    return float(np.sqrt(c["b0"] * A.U_MAX / (((0.5 + c["sigma"]) * Q12) ** c["a2"])))


def job(j):
    c = j["case"]
    A.A1 = (c["a2"] + 1) / 2
    A.A2 = c["a2"]
    T.B_NOM = c["b0"]
    T.SENSORS["X"] = (c["sigma"], 0.0, 0.0)
    with np.errstate(all="ignore"):
        sc = T.setup(A.scenarios(T.RUNS, seed=j["seed"]), "std")
        r = T.simulate(sc, j["lam"], "frac", "full", 0.0, 12, T.SENSORS["X"], 100 + j["seed"])
    r.update(lam=j["lam"], seed=j["seed"], **c)
    return r


def main():
    jobs = []
    for c in CASES:
        p = predict(c)
        print(c, "predicted ceiling", round(p, 1), flush=True)
        for g in p * 1.08 ** np.arange(-7, 8):
            for s in (7, 8, 9):
                jobs.append(dict(case=c, lam=float(g), seed=s))
    print(len(jobs), flush=True)
    t0 = time.time()
    with mp.Pool(4, maxtasksperchild=4) as p:
        rows = p.map(job, jobs, chunksize=1)
    (A.OUT / "test4.json").write_text(json.dumps(clean(rows), indent=0))
    print("case, predicted, measured, ratio", flush=True)
    for c in CASES:
        rs = [r for r in rows if all(r[k] == c[k] for k in ("sigma", "a2", "b0"))]
        gains = sorted({r["lam"] for r in rs})
        med = [(g, float(np.median([r["sat"] for r in rs if r["lam"] == g]))) for g in gains]
        m = next((g for g, s in med if s > 0.10), None)
        print(c, round(predict(c), 1), None if m is None else round(m, 1), None if m is None else round(m / predict(c), 3), [round(100 * s) for _, s in med], flush=True)
    print("saved", round(time.time() - t0), flush=True)


if __name__ == "__main__":
    main()

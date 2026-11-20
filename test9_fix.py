"""Test 9: fixes for the two-link violations. Rules fixed before any run.

Diagnosis (Test 8, seeds 30 to 34): all 11 angle violations had a reference within 0.02 rad of, or past, a limit.
All 18 velocity violations were in start-up (0.15 to 0.44 s), with no saturation and no pulse nearby.

Changes applied to all cases: references kept 0.1 rad inside the angle limits. Gains unchanged (24, 29), low-pass.
Clearance from the steady-error formula plus twice the estimate RMS, using values from seeds 30 to 34:
  rho_j = 2 D_j^1.5 / lam_j^2 + 2 E_j, with D = (10.5, 21.9), E = (0.0487, 0.0591).
Cases: X0 5% clearance; X1 rho_j at all times; X2 rho_j until 0.6 s, falling to 5% by 0.8 s.
New seeds 40 to 44 (500 runs), sensors S2 and S3.
  Y1: X0 shows whether the reference fix alone meets <= 1% violating runs.
  Y2: a case passes safety with <= 1% velocity AND <= 1% angle violations under S2 and S3.
  Y3: X2 passes cost if tracking RMS is within 25% of X0 on both joints under S2.
All outcomes are reported. No further tuning follows this test.
"""
import json
import multiprocessing as mp

import numpy as np

import test7_two_link as T
from run_all import clean

LAM = np.array([24.0, 29.0])
D = np.array([10.5, 21.9])
E = np.array([0.0487, 0.0591])
RHO = 2 * D ** 1.5 / LAM ** 2 + 2 * E
R5 = 0.05 * (T.W_U - T.W_L)
CASES = {"X0": dict(), "X1": dict(RHO_ABS=RHO), "X2": dict(RHO_T=(RHO, R5, 0.6, 0.8))}


def job(j):
    T.REF_SAFE = True
    T.RHO_ABS = CASES[j["case"]].get("RHO_ABS")
    T.RHO_T = CASES[j["case"]].get("RHO_T")
    sc = T.scen(T.RUNS, j["seed"])
    n = T.RUNS
    d = {k: np.full((n, 2), np.nan) for k in ("t", "ta", "e2", "x2", "vc", "sat", "th", "yd")}
    for k in ("wmax", "amax"):
        d[k] = np.full((n, 2), -np.inf)
    d["dmax"], d["emax"] = np.zeros((n, 2)), np.zeros((n, 2))
    with np.errstate(all="ignore"):
        r = T.simulate(sc, LAM.tolist(), ["frac"] * 2, ["lp"] * 2, T.SENS[j["sensor"]], 200 + j["seed"], diag=d)
    r.update(j, tfirst=[float(x) for x in d["t"][~np.isnan(d["t"])]], wexc=float(d["wmax"].max()), aexc=float(d["amax"].max()))
    return r


if __name__ == "__main__":
    print("rho", RHO.round(3), "5%", R5.round(3), flush=True)
    jobs = [dict(case=c, sensor=s, seed=x) for c in CASES for s in ("S2", "S3") for x in (40, 41, 42, 43, 44)]
    with mp.Pool(4, maxtasksperchild=4) as p:
        rows = p.map(job, jobs, chunksize=1)
    json.dump(clean(dict(rho=RHO.tolist(), rows=rows)), open("results/test9.json", "w"), indent=0)
    for c in CASES:
        for s in ("S2", "S3"):
            rs = [r for r in rows if r["case"] == c and r["sensor"] == s]
            tf = sum((r["tfirst"] for r in rs), [])
            print(c, s, "viol", sum(r["viol"] for r in rs), "ang", sum(r["ang"] for r in rs), "/", sum(r["runs"] for r in rs), "sat", np.round(100 * np.median([r["sat"] for r in rs], 0), 1), "rms", np.round(np.mean([r["rms"] for r in rs], 0), 5),
                  "e2", np.round(np.mean([r["e2"] for r in rs], 0), 4), "wexc", round(max(r["wexc"] for r in rs), 3), "aexc", round(max(r["aexc"] for r in rs), 4), "first", (round(min(tf), 3), round(max(tf), 3)) if tf else None, flush=True)

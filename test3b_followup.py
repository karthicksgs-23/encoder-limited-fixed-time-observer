"""Follow-up to Test 3, added after seeing its results. Exploratory, not fixed in advance,
except that the torque prediction below was written before this was run.

T: torque dependence of the gain rule. Base design, 12-bit S2, torque limit 4, 8, 16 N m.
   Rule with no fitted constant: lam_c = sqrt(b0 * u_max / q^a2) -> 14.3, 20.2, 28.6.
N: 'no cancel' design, omitted from Test 3, under S2.
C: clearance sensitivity. All designs under S2 with command clearance 5% instead of 1%.
"""
import json
import multiprocessing as mp
import time

import numpy as np

import arm_margin_test as A
import test3_realistic as T
from run_all import clean

T.DESIGNS["no cancel"] = ("frac", "none", 0.0)
GT = (10.0, 12.0, 14.0, 17.0, 20.0, 24.0, 28.0, 34.0, 40.0)


def job(j):
    T.U_MAX = j["umax"]
    T.CLEAR = j["clear"]
    with np.errstate(all="ignore"):
        sc = T.setup(A.scenarios(T.RUNS, seed=j["seed"]), "std")
        obs, cancel, dz = T.DESIGNS[j["design"]]
        r = T.simulate(sc, j["lam"], obs, cancel, dz, 12, T.SENSORS["S2"], 100 + j["seed"])
    r.update(j)
    return r


def main():
    jobs = []
    for s in (7, 8, 9):
        for um in (4.0, 8.0, 16.0):
            for lam in GT:
                jobs.append(dict(part="T", design="base", lam=lam, seed=s, umax=um, clear=0.01))
        for lam in T.GAINS:
            jobs.append(dict(part="N", design="no cancel", lam=lam, seed=s, umax=8.0, clear=0.01))
        for d in T.DESIGNS:
            for lam in (14.0, 20.0, 28.0, 40.0, 57.0):
                jobs.append(dict(part="C", design=d, lam=lam, seed=s, umax=8.0, clear=0.05))
    print(len(jobs), flush=True)
    t0 = time.time()
    with mp.Pool(4, maxtasksperchild=1) as p:
        rows = p.map(job, jobs, chunksize=1)
    (A.OUT / "test3b.json").write_text(json.dumps(clean(rows), indent=0))
    print("saved", round(time.time() - t0), flush=True)


if __name__ == "__main__":
    main()

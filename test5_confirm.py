"""Test 5: confirmation on unseen scenario seeds.

Settings were chosen from Tests 3 and 3b, then frozen here before this run:
gains and designs below, 12-bit encoder, 1 kHz, 5% command clearance, seeds 20 to 24
(100 runs each, never used before). Sensor S2 and S3. All outcomes are reported as they come.
Time responses for one scenario are saved for the figures.
"""
import json
import multiprocessing as mp

import numpy as np

import arm_margin_test as A
import test3_realistic as T
from run_all import clean

D = [("base", 14.0), ("low-pass", 20.0), ("dead zone 1", 20.0), ("linear", 57.0), ("base", 40.0)]


def job(j):
    T.CLEAR = 0.05
    with np.errstate(all="ignore"):
        sc = T.setup(A.scenarios(T.RUNS, seed=j["seed"]), "std")
        obs, cancel, dz = T.DESIGNS[j["design"]]
        r = T.simulate(sc, j["lam"], obs, cancel, dz, 12, T.SENSORS[j["sensor"]], 100 + j["seed"])
    r.update(j)
    return r


def traces():
    T.CLEAR = 0.05
    sc = T.setup(A.scenarios(T.RUNS, seed=20), "std")
    ok = np.where((sc["pt"] > 1.2) & (sc["pt"] < 3.2) & (np.abs(sc["pa"]) > 2.0) & (np.abs(sc["w0"]) > 0.3) & (sc["pd"] > 0.1))[0]
    idx = int(ok[0])
    out = dict(idx=idx, pulse=[float(sc["pt"][idx]), float(sc["pd"][idx]), float(sc["pa"][idx])], w0=float(sc["w0"][idx]))
    for d, lam in (("low-pass", 20.0), ("base", 40.0), ("linear", 57.0)):
        tr = []
        o, c, z = T.DESIGNS[d]
        with np.errstate(all="ignore"):
            T.simulate({k: (v.copy() if hasattr(v, "copy") else v) for k, v in sc.items()}, lam, o, c, z, 12, T.SENSORS["S2"], 120, trace=tr, idx=idx)
        out[f"{d} {lam:g}"] = tr[::2]
    return out


def main():
    jobs = [dict(design=d, lam=l, sensor=s, seed=x) for d, l in D for s in ("S2", "S3") for x in (20, 21, 22, 23, 24)]
    with mp.Pool(4, maxtasksperchild=4) as p:
        rows = p.map(job, jobs, chunksize=1)
    (A.OUT / "test5.json").write_text(json.dumps(clean(rows), indent=0))
    (A.OUT / "traces.json").write_text(json.dumps(clean(traces())))
    for s in ("S2", "S3"):
        for d, l in D:
            rs = [r for r in rows if r["design"] == d and r["lam"] == l and r["sensor"] == s]
            n = sum(r["runs"] for r in rs)
            print(s, d, l, "viol", sum(r["viol"] for r in rs), "/", n, "ang", sum(r["ang"] for r in rs), "sat", round(100 * float(np.median([r["sat"] for r in rs])), 1),
                  "rms", round(float(np.mean([r["rms"] for r in rs])), 5), "e2", round(float(np.mean([r["e2"] for r in rs])), 4), "utv", round(float(np.mean([r["utv"] for r in rs]))),
                  "exc", round(max(r["excess"] for r in rs), 3), flush=True)


if __name__ == "__main__":
    main()

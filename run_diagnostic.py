"""Part C: one bounded diagnostic after Parts A and B. Exploratory only.

Q1 Is there a setting where the margin both fits and is needed, if the disturbance
   bound uses D = 8, set from the pilot maximum of about 7.4, and not the guaranteed 28.8?
Q2 Does the encoder breakdown come from estimation error or from something else?
"""
import json
import multiprocessing as mp
import sys

import numpy as np

import arm_margin_test as A
from run_all import ENC12, clean

ENC16 = 2 * np.pi / 65536
D_OBS = 8.0


def job(j):
    A.RHO = 0.01
    sc = A.scenarios(300)
    s, _ = A.simulate(sc, A.margin_fn(j["kind"], j["mp"]), j["lam"], dt=2e-4, quant=j["quant"])
    r = A.summary(j["name"], s, 300)
    first = s["first"][~np.isnan(s["first"])]
    r.update(lam=j["lam"], quant=j["quant"], **j["mp"], first_min=float(first.min()) if first.size else None,
             first_max=float(first.max()) if first.size else None)
    print("lam", j["lam"], "q", round(j["quant"], 6), j["name"], "viol", r["violating_runs"], "excess", round(r["worst_excess"], 3), "estviol", r["est_violating_runs"],
          "breach", r["margin_breaches"], "rms", round(r["rms_median"], 4), "sat", round(r["sat_time_median"], 3), "first", r["first_min"], r["first_max"], flush=True)
    return r


def main():
    jobs, margins = [], {}
    for lam in (10.0, 15.0, 20.0):
        tt, env = A.error_envelope(lam, D_OBS)
        margins[lam] = A.margin_from_envelope(tt, env)
        print("margin D=8", lam, margins[lam], flush=True)
        for q in (0.0, ENC12):
            for name, kind in (("no margin", "none"), ("constant margin", "constant"), ("time-varying margin", "tv")):
                if kind != "none" and margins[lam]["m0"] >= 0.81:
                    continue
                jobs.append(dict(name=name, kind=kind, lam=lam, quant=q, mp=margins[lam]))
    for lam in (40.0, 60.0):
        tt, env = A.error_envelope(lam, D_OBS)
        mpm = A.margin_from_envelope(tt, env)
        for q in (ENC12, ENC16):
            jobs.append(dict(name="no margin", kind="none", lam=lam, quant=q, mp=mpm))
    print(len(jobs), "jobs", flush=True)
    with mp.Pool(4) as p:
        rows = p.map(job, jobs, chunksize=1)
    (A.OUT / "diagnostic.json").write_text(json.dumps(clean(dict(margins={str(k): v for k, v in margins.items()}, rows=rows)), indent=1))
    print("saved", flush=True)


if __name__ == "__main__":
    main()

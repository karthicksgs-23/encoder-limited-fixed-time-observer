"""Test 2: is the encoder breakdown a real resolution limit or a design artefact?

Fixed before any result was seen.

Part F, simple fixes at 12-bit, gains 20, 40, 60. Variants: base, no injection
cancellation in u, low-pass injection in u (30 Hz), one-count dead zone, linear observer.
  F1 artefact rule: if 'no cancel' or 'low-pass' gives <= 1% violating runs and
     < 10% saturated time at gains 40 AND 60, with median tracking RMS no worse than
     2x the base at gain 20, the breakdown is a design artefact. Then no
     resolution-limit claim may be made for that design.

Part L, gain ceiling across resolutions (10, 12, 14, 16 bit; gains 10 to 113).
Ceiling = lowest tested gain with median saturated time above 10% of the run.
  Prediction: the injection step equals the torque authority at
     fractional observer: lam_c^2 * q^a2 ~ c * b0 * u_max  -> slope d(log lam_c)/d(log q) = -a2/2 = -0.25
     linear observer:     2 lam_c^2 * q  ~ c * b0 * u_max  -> slope = -0.50
  F2 law rule: fitted slope within 0.10 of the prediction for each observer.
  F3 value rule: at the highest gain below each ceiling, report the time for |e2| < 0.05
     and the estimate error after t = 1 s. Fixed-time terms are worthwhile only if the
     fractional observer is faster at equal resolution.
"""
import json
import multiprocessing as mp
import sys
import time

import numpy as np

import arm_margin_test as A
from run_all import clean

BITS = (10, 12, 14, 16)
GAINS = (10.0, 14.0, 20.0, 28.0, 40.0, 57.0, 80.0, 113.0)
TT = 4.0


def q(bits):
    return 2 * np.pi / 2 ** bits


def job(j):
    A.RHO = 0.01
    if j["obs"] == "linear":
        A.A1 = A.B1 = A.A2 = A.B2 = 1.0
    sc = A.scenarios(j["runs"])
    s, _ = A.simulate(sc, A.margin_fn("none", {}), j["lam"], T=TT, dt=2e-4, quant=q(j["bits"]), cancel=j["cancel"], dz=j["dz"] * q(j["bits"]))
    r = A.summary(j["variant"], s, j["runs"])
    r.update(part=j["part"], variant=j["variant"], obs=j["obs"], lam=j["lam"], bits=j["bits"], sat_share=float(np.median(s["sat"])) / TT,
             t_e2=float(np.nanmedian(s["t_e2"])), e2_rms=float(np.median(np.sqrt(s["e2rms"] / max(s["e2cnt"], 1)))))
    print(j["part"], j["variant"], "lam", j["lam"], "bits", j["bits"], "viol", r["violating_runs"], "sat", round(r["sat_share"], 3), "rms", round(r["rms_median"], 4),
          "te2", round(r["t_e2"], 4), "e2rms", round(r["e2_rms"], 4), flush=True)
    return r


def main():
    jobs = []
    V = [("base", "frac", "full", 0.0), ("no cancel", "frac", "none", 0.0), ("low-pass", "frac", "lp", 0.0), ("dead zone", "frac", "full", 1.0), ("linear", "linear", "full", 0.0)]
    for lam in (20.0, 40.0, 60.0):
        for v, o, c, d in V:
            jobs.append(dict(part="F", variant=v, obs=o, cancel=c, dz=d, lam=lam, bits=12, runs=200))
    for o in ("frac", "linear"):
        for b in BITS:
            for lam in GAINS:
                jobs.append(dict(part="L", variant=o, obs=o, cancel="full", dz=0.0, lam=lam, bits=b, runs=60))
    print(len(jobs), "jobs", flush=True)
    t0 = time.time()
    with mp.Pool(4, maxtasksperchild=1) as p:
        rows = p.map(job, jobs, chunksize=1)
    (A.OUT / "test2.json").write_text(json.dumps(clean(rows), indent=1))
    print("saved", round(time.time() - t0), flush=True)


if __name__ == "__main__":
    main()

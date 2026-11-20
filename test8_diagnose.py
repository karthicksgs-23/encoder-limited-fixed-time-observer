"""Test 8: diagnosis of two-link violations. Same design, seeds and sensor as Test 7 (S2, seeds 30 to 34)."""
import json
import multiprocessing as mp

import numpy as np

import test7_two_link as T


def job(seed):
    sc = T.scen(T.RUNS, seed)
    n = T.RUNS
    d = {k: np.full((n, 2), np.nan) for k in ("t", "ta", "e2", "x2", "vc", "sat", "th", "yd")}
    d["wmax"] = np.full((n, 2), -np.inf)
    d["dmax"] = np.zeros((n, 2))
    d["emax"] = np.zeros((n, 2))
    d["amax"] = np.full((n, 2), -np.inf)
    with np.errstate(all="ignore"):
        T.simulate(sc, [24.0, 29.0], ["frac"] * 2, ["lp"] * 2, T.SENS["S2"], 200 + seed, diag=d)
    tt = np.arange(0, T.T_END, 1e-3)[None, :, None]
    yd = sc["off"][:, None, :] + sc["amp"][:, None, :] * np.sin(sc["w"][:, None, :] * tt + sc["ph"][:, None, :])
    d["ref_lo"] = yd.min(1)
    d["ref_hi"] = yd.max(1)
    d["pk"] = sc["amp"] * sc["w"]
    d["pt"], d["pd"], d["pa"] = sc["pt"], sc["pd"], sc["pa"]
    d["w0"], d["th0"] = sc["w0"], sc["th0"]
    return {k: v.tolist() for k, v in d.items()}


if __name__ == "__main__":
    with mp.Pool(4) as p:
        R = p.map(job, (30, 31, 32, 33, 34))
    D = {k: np.concatenate([np.array(r[k]) for r in R]) for k in R[0]}
    json.dump({k: v.tolist() for k, v in D.items()}, open("results/test8.json", "w"))
    bad_ref = (D["ref_lo"] <= T.TH_L) | (D["ref_hi"] >= T.TH_U)
    print("runs", len(D["t"]))
    print("reference outside angle limits, per joint", bad_ref.sum(0), "runs any", bad_ref.any(1).sum())
    vel = ~np.isnan(D["t"]); ang = ~np.isnan(D["ta"])
    print("velocity violations per joint", vel.sum(0), "any", vel.any(1).sum())
    print("angle violations per joint", ang.sum(0), "any", ang.any(1).sum())
    for name, m in (("velocity", vel.any(1)), ("angle", ang.any(1))):
        print(name, "violating runs", m.sum(), "of which reference infeasible", (m & bad_ref.any(1)).sum(), "| infeasible runs that violated", (m & bad_ref.any(1)).sum(), "of", bad_ref.any(1).sum())
    both = vel.any(1) & ~bad_ref.any(1)
    print("velocity violations with feasible reference:", both.sum())
    print("peak |delta| per joint: max", D["dmax"].max(0).round(1), "p95", np.percentile(D["dmax"], 95, 0).round(1), "median", np.median(D["dmax"], 0).round(1))
    print("peak |e2| after 0.1 s: max", D["emax"].max(0).round(3), "p95", np.percentile(D["emax"], 95, 0).round(3), "median", np.median(D["emax"], 0).round(3))
    print("clearance now", (0.05 * (T.W_U - T.W_L)).round(3))
    for i in np.where(ang.any(1) & ~bad_ref.any(1))[0]:
        j = int(np.nanargmin(np.where(ang[i], D["ta"][i], np.inf)))
        print(f"angle run {i} joint {j+1} t {D['ta'][i,j]:.3f} excess {D['amax'][i,j]:.4f} ref range {D['ref_lo'][i,j]:.3f} to {D['ref_hi'][i,j]:.3f} limit {T.TH_L[j]} to {T.TH_U[j]}")
    for i in []:
        j = int(np.nanargmin(np.where(vel[i], D["t"][i], np.inf)))
        pul = D["pt"][i, 0] <= D["t"][i, j] <= D["pt"][i, 0] + D["pd"][i, 0] + .3 or D["pt"][i, 1] <= D["t"][i, j] <= D["pt"][i, 1] + D["pd"][i, 1] + .3
        print(f"run {i} joint {j+1} t {D['t'][i,j]:.3f} excess {D['wmax'][i,j]:.3f} w_est {D['x2'][i,j]:.3f} e2 {D['e2'][i,j]:.3f} vc {D['vc'][i,j]:.3f} sat {int(D['sat'][i,j])} badref {bool(bad_ref[i].any())} th {D['th'][i,j]:.3f} yd {D['yd'][i,j]:.3f} near_pulse {bool(pul)} w0 {D['w0'][i,j]:.2f}")

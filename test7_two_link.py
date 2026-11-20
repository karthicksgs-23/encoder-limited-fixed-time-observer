"""Test 7: second benchmark, a geared two-link arm in a vertical plane with coupled joints.

Each joint has its own 12-bit encoder, observer, fuzzy system and controller (decentralised).
Inertial, Coriolis and gravity coupling is unknown to the controller. 1 kHz, one-sample delay, 5% clearance.

Fixed before any run. Predictions come from the single-link rule with no change:
  lam_c,j = sqrt(b0_j * u_max,j / (kappa q)^a2), kappa = 1 under S2.
  W1 gain rule: with the other joint held at a safe gain, each joint's ceiling (first gain on an 8% grid
     with median saturated share of that joint above 10%) lies within 15% of its prediction.
  W2 safety: the proposed design, with gains set from the rule and NOT tuned, gives <= 1% runs with any
     velocity-limit violation on seeds 30 to 34 (500 runs) under S2 and S3.
  W3 value: on the same seeds under S2, the proposed design has mean tracking RMS no worse than a linear
     observer whose gain is selected on seed 7 (lowest estimate RMS with <= 1% violations and < 10% saturation).
  W4 contrast: full cancellation at twice the rule's gain is reported as it comes.
All outcomes are reported whether they pass or fail.
"""
import json
import multiprocessing as mp
import sys
import time

import numpy as np

import arm_margin_test as A
from arm_margin_test import C_FILT, K1, K1B, K2, K2B, P, Q, h, h_x, sig
from run_all import clean

TS, SUB, T_END, RUNS, CLEAR = 1e-3, 10, 4.0, 100, 0.05
G9 = 9.81
M1, M2, L1, L2 = 1.2, 0.8, 0.35, 0.30
LC1, LC2 = L1 / 2, L2 / 2
I1, I2 = M1 * L1 ** 2 / 12, M2 * L2 ** 2 / 12
JM = np.array([0.30, 0.10])         # reflected rotor inertia
FV = np.array([0.30, 0.15])
FC = np.array([0.15, 0.10])
U = np.array([12.0, 4.0])
TH_L, TH_U = np.array([-0.9, -0.7]), np.array([1.3, 1.1])
W_L, W_U = np.array([-1.0, -1.1]), np.array([1.4, 1.5])
REF_SAFE = False                    # keep references 0.1 rad inside the angle limits
RHO_T = None                        # optional (start-up clearance, final clearance, t1, t2)
RHO_ABS = None                      # optional clearance per joint, rad/s
COS0 = 0.9                          # nominal elbow cosine for b0


def m_nom():
    m11 = I1 + I2 + M1 * LC1 ** 2 + M2 * (L1 ** 2 + LC2 ** 2 + 2 * L1 * LC2 * COS0) + JM[0]
    m22 = I2 + M2 * LC2 ** 2 + JM[1]
    return np.array([m11, m22])


B0 = 1 / m_nom()
Q12 = 2 * np.pi / 4096
LAMC = np.sqrt(B0 * U / Q12 ** 0.5)
FMAX = np.array([14.0, 12.0])
GAMMA = 2.0
CY = np.stack([np.linspace(TH_L[j], TH_U[j], 5) for j in range(2)])
CW = np.stack([np.linspace(W_L[j], W_U[j], 5) for j in range(2)])
SY = 0.6 * (CY[:, 1] - CY[:, 0])
SW = 0.6 * (CW[:, 1] - CW[:, 0])
DESIGN = {"base": ("frac", "full"), "low-pass": ("frac", "lp"), "linear": ("linear", "full")}
SENS = {"S2": (0.5, 0.0, 0.0), "S3": (0.5, 2.0, 0.5)}


def scen(n, seed):
    r = np.random.default_rng(seed)
    u = lambda s, *sh: r.uniform(1 - s, 1 + s, (n, *sh))
    w = r.uniform(1.2, 2.0, (n, 2))
    pk = r.uniform(0.5, 0.8, (n, 2)) * np.minimum(W_U, -W_L)
    off = r.uniform(-.1, .1, (n, 2))
    ph = r.uniform(0, 2 * np.pi, (n, 2))
    amp = pk / w
    if REF_SAFE:
        amp = np.minimum(amp, np.minimum(TH_U - 0.1 - off, off - TH_L - 0.1))
    return dict(m1=M1 * u(.2), m2=M2 * u(.2), jm=JM * u(.1, 2), fv=FV * u(.3, 2), fc=FC * u(.5, 2),
                th0=r.uniform(-.15, .15, (n, 2)), w0=r.uniform(-.4, .4, (n, 2)), amp=amp, w=w, ph=ph, off=off,
                pa=r.uniform(-1, 1, (n, 2)) * np.array([2.0, 0.8]), pt=r.uniform(1, 3.5, (n, 2)), pd=r.uniform(.05, .2, (n, 2)))


def accel(sc, th, w, tau):
    m1, m2 = sc["m1"], sc["m2"]
    i1, i2 = m1 * L1 ** 2 / 12, m2 * L2 ** 2 / 12
    c2, s2 = np.cos(th[:, 1]), np.sin(th[:, 1])
    a11 = i1 + i2 + m1 * LC1 ** 2 + m2 * (L1 ** 2 + LC2 ** 2 + 2 * L1 * LC2 * c2) + sc["jm"][:, 0]
    a12 = i2 + m2 * (LC2 ** 2 + L1 * LC2 * c2)
    a22 = i2 + m2 * LC2 ** 2 + sc["jm"][:, 1]
    hh = m2 * L1 * LC2 * s2
    c1 = -hh * (2 * w[:, 0] * w[:, 1] + w[:, 1] ** 2)
    c2_ = hh * w[:, 0] ** 2
    g1 = (m1 * LC1 + m2 * L1) * G9 * np.cos(th[:, 0] + np.pi / 2) + m2 * LC2 * G9 * np.cos(th[:, 0] + th[:, 1] + np.pi / 2)
    g2 = m2 * LC2 * G9 * np.cos(th[:, 0] + th[:, 1] + np.pi / 2)
    r = tau - sc["fv"] * w - sc["fc"] * np.tanh(20 * w)
    r1, r2 = r[:, 0] - c1 - g1, r[:, 1] - c2_ - g2
    det = a11 * a22 - a12 ** 2
    return np.stack([(a22 * r1 - a12 * r2) / det, (a11 * r2 - a12 * r1) / det], 1)


def simulate(sc, lam, obs, cancel, sensor, seed, trace=None, idx=0, diag=None):
    """lam, obs and cancel each hold one entry per joint."""
    n = sc["th0"].shape[0]
    bits_q = Q12
    ns, bs, dr = sensor
    rng = np.random.default_rng(seed)
    lin = np.array([o == "linear" for o in obs])
    a1 = np.where(lin, 1.0, A.A1); b1 = np.where(lin, 1.0, A.B1)
    a2 = np.where(lin, 1.0, A.A2); b2 = np.where(lin, 1.0, A.B2)
    lam = np.asarray(lam, float)
    l1, l2 = 2 * lam, lam ** 2
    lp = np.array([c == "lp" for c in cancel])
    bias = rng.uniform(-bs, bs, (n, 2)) * bits_q
    drift = rng.uniform(-dr, dr, (n, 2)) * bits_q
    th, w = sc["th0"].copy(), sc["w0"].copy()
    x1, x2 = th.copy(), np.zeros((n, 2))
    theta = np.zeros((n, 2, 25))
    vh, i2f, ua, ul = (np.zeros((n, 2)) for _ in range(4))
    steps = int(T_END / TS)
    o1h, o1l, oh, ol = TH_U, -TH_L, W_U, -W_L
    rho = CLEAR * (oh + ol) if RHO_ABS is None else RHO_ABS
    lo, hi = -ol + 1e-3, oh - 1e-3
    wv, tv = np.zeros(n), np.zeros(n)
    sq, e2s = np.zeros((n, 2)), np.zeros((n, 2))
    cnt = cnt2 = 0
    tvu, sat = np.zeros((n, 2)), np.zeros((n, 2))
    alp = TS * 2 * np.pi * 30.0
    dts = TS / SUB
    for k in range(steps):
        t = k * TS
        yd = sc["off"] + sc["amp"] * np.sin(sc["w"] * t + sc["ph"])
        dyd = sc["amp"] * sc["w"] * np.cos(sc["w"] * t + sc["ph"])
        y = th + bias + drift * t + ns * bits_q * rng.standard_normal((n, 2))
        y = np.round(y / bits_q) * bits_q
        yc = np.clip(y, TH_L + 1e-3, TH_U - 1e-3)
        e1 = y - x1
        i1 = l1 * (sig(e1, a1) + sig(e1, b1))
        i2 = l2 * (sig(e1, a2) + sig(e1, b2))
        psi = np.exp(-((yc[:, :, None, None] - CY[None, :, :, None]) / SY[None, :, None, None]) ** 2 - ((x2[:, :, None, None] - CW[None, :, None, :]) / SW[None, :, None, None]) ** 2).reshape(n, 2, 25)
        psi /= psi.sum(2, keepdims=True) + 1e-9
        fh = (theta * psi).sum(2)
        z1 = h(yc, o1h, o1l) - h(yd, o1h, o1l)
        v = (h_x(yd, o1h, o1l) * dyd - K1 * sig(z1, P) - K1B * sig(z1, Q)) / h_x(yc, o1h, o1l)
        if RHO_T is not None:
            a_, b_, t1_, t2_ = RHO_T
            rho = a_ + (b_ - a_) * min(max((t - t1_) / (t2_ - t1_), 0.0), 1.0)
        vc = np.clip(v, -ol + rho, oh - rho)
        dvh = C_FILT * (vc - vh)
        vhc, xc = np.clip(vh, lo, hi), np.clip(x2, lo, hi)
        z2 = h(xc, oh, ol) - h(vhc, oh, ol)
        dx = (h_x(vhc, oh, ol) * dvh - K2 * sig(z2, P) - K2B * sig(z2, Q)) / h_x(xc, oh, ol)
        i2f = i2f + alp * (i2 - i2f)
        iu = np.where(lp, i2f, i2)
        ur = (dx - fh - iu) / B0
        un = np.clip(ur, -U, U)
        sat += np.abs(ur) >= U
        tvu += np.abs(un - ul)
        ul = un
        for j in range(SUB):
            tj = t + j * dts
            d = 0.2 * np.sin(2 * tj) + np.where((tj >= sc["pt"]) & (tj < sc["pt"] + sc["pd"]), sc["pa"], 0.0)
            dw = accel(sc, th, w, ua + d)
            th = th + dts * w
            w = w + dts * dw
            wv = np.maximum(wv, np.maximum(w - W_U, W_L - w).max(1))
            tv = np.maximum(tv, np.maximum(th - TH_U, TH_L - th).max(1))
        if diag is not None:
            ov = np.maximum(w - W_U, W_L - w)
            new = (ov > 0) & np.isnan(diag["t"])
            if new.any():
                diag["t"][new] = t
                diag["e2"][new] = (w - x2)[new]
                diag["x2"][new] = x2[new]
                diag["vc"][new] = vc[new]
                diag["sat"][new] = (np.abs(ur) >= U)[new]
                diag["th"][new] = th[new]
                diag["yd"][new] = yd[new]
            ao = np.maximum(th - TH_U, TH_L - th)
            na = (ao > 0) & np.isnan(diag["ta"])
            diag["ta"][na] = t
            diag["wmax"] = np.maximum(diag["wmax"], ov)
            diag["dmax"] = np.maximum(diag["dmax"], np.abs(dw - B0 * ua - fh))
            diag["emax"] = np.maximum(diag["emax"], np.abs(w - x2) if t >= 0.1 else 0.0)
            diag["amax"] = np.maximum(diag["amax"], ao)
        if trace is not None:
            trace.append((t, *th[idx], *yd[idx], *w[idx], *x2[idx], *un[idx]))
        if t >= 2.0:
            sq += (th - yd) ** 2
            cnt += 1
        if t >= 1.0:
            e2s += (w - x2) ** 2
            cnt2 += 1
        theta = np.clip(theta + TS * GAMMA * psi * i2[:, :, None], -FMAX[None, :, None], FMAX[None, :, None])
        x1 = x1 + TS * (x2 + i1)
        x2 = x2 + TS * (B0 * ua + fh + i2)
        vh = vh + TS * dvh
        ua = un
    return dict(viol=int((wv > 0).sum()), ang=int((tv > 0).sum()), excess=float(wv.max()), runs=n, finite=bool(np.isfinite(th).all()),
                sat=[float(x) for x in np.median(sat / steps, 0)], rms=[float(x) for x in np.median(np.sqrt(sq / cnt), 0)], e2=[float(x) for x in np.median(np.sqrt(e2s / cnt2), 0)],
                utv=[float(x) for x in np.median(tvu / T_END, 0)])


def job(j):
    with np.errstate(all="ignore"):
        r = simulate(scen(RUNS, j["seed"]), j["lam"], j["obs"], j["cancel"], SENS[j["sensor"]], 200 + j["seed"])
    r.update(j)
    return r


def main():
    print("b0", B0.round(3), "predicted ceilings", LAMC.round(1), flush=True)
    jobs = []
    grid = 1.08 ** np.arange(-7, 8)
    for jt in (0, 1):
        for g in grid:
            for s in (7, 8, 9):
                lam = [15.0, 15.0]
                lam[jt] = float(LAMC[jt] * g)
                cancel = ["lp", "lp"]
                cancel[jt] = "full"
                jobs.append(dict(part="W1", joint=jt, lam=lam, obs=["frac", "frac"], cancel=cancel, sensor="S2", seed=s))
    lin = [28.0, 40.0, 57.0, 68.0, 80.0, 96.0]
    for g in lin:
        jobs.append(dict(part="sel", lam=[g, g * float(LAMC[1] / LAMC[0])], obs=["linear"] * 2, cancel=["full"] * 2, sensor="S2", seed=7))
    t0 = time.time()
    with mp.Pool(4, maxtasksperchild=4) as p:
        rows = p.map(job, jobs, chunksize=1)
    print("stage 1", round(time.time() - t0), flush=True)
    c = [r for r in rows if r["part"] == "sel" and r["viol"] / r["runs"] <= .01 and max(r["sat"]) < .10]
    # No candidate met the rule in the first run. Deviation, added after seeing that: use the fewest violations.
    sel = min(c, key=lambda r: sum(r["e2"])) if c else min((r for r in rows if r["part"] == "sel"), key=lambda r: r["viol"])
    print("linear candidates", [(round(r["lam"][0]), r["viol"], [round(x, 3) for x in r["sat"]], [round(x, 4) for x in r["e2"]]) for r in rows if r["part"] == "sel"], flush=True)
    lamp = [float(round(x)) for x in LAMC]
    D = [("proposed", lamp, ["frac"] * 2, ["lp"] * 2), ("base 2x", [2 * x for x in lamp], ["frac"] * 2, ["full"] * 2)]
    if sel:
        D.append(("linear", sel["lam"], ["linear"] * 2, ["full"] * 2))
    jobs = [dict(part="final", name=n, lam=l, obs=o, cancel=c_, sensor=s, seed=x) for n, l, o, c_ in D for s in ("S2", "S3") for x in (30, 31, 32, 33, 34)]
    with mp.Pool(4, maxtasksperchild=4) as p:
        rows += p.map(job, jobs, chunksize=1)
    tr = {}
    sc = scen(RUNS, 30)
    idx = int(np.where((np.abs(sc["w0"]).max(1) > .3) & (np.abs(sc["pa"][:, 0]) > 1.2) & (sc["pt"][:, 0] > 1.2) & (sc["pt"][:, 0] < 3) & (sc["pd"][:, 0] > .1))[0][0])
    for n, l, o, c_ in D[:2]:
        x = []
        with np.errstate(all="ignore"):
            simulate(sc, l, o, c_, SENS["S2"], 230, trace=x, idx=idx)
        tr[n] = x[::2]
    tr["pulse"] = [float(sc["pt"][idx, 0]), float(sc["pd"][idx, 0])]
    (A.OUT / "test7.json").write_text(json.dumps(clean(dict(b0=B0.tolist(), lamc=LAMC.tolist(), linear_eligible=bool(c), rows=rows)), indent=0))
    (A.OUT / "traces7.json").write_text(json.dumps(clean(tr)))
    print("W1")
    for jt in (0, 1):
        rs = [r for r in rows if r["part"] == "W1" and r["joint"] == jt]
        g = sorted({r["lam"][jt] for r in rs})
        med = [(x, float(np.median([r["sat"][jt] for r in rs if r["lam"][jt] == x]))) for x in g]
        m = next((x for x, s in med if s > .10), None)
        print("joint", jt + 1, "pred", round(LAMC[jt], 1), "measured", m and round(m, 1), "ratio", m and round(m / LAMC[jt], 3), [round(100 * s) for _, s in med], flush=True)
    print("final")
    for n, *_ in D:
        for s in ("S2", "S3"):
            rs = [r for r in rows if r["part"] == "final" and r["name"] == n and r["sensor"] == s]
            print(n, rs[0]["lam"], s, "viol", sum(r["viol"] for r in rs), "/", sum(r["runs"] for r in rs), "ang", sum(r["ang"] for r in rs), "sat", np.round(100 * np.median([r["sat"] for r in rs], 0), 1),
                  "rms", np.round(np.mean([r["rms"] for r in rs], 0), 5), "e2", np.round(np.mean([r["e2"] for r in rs], 0), 4), "utv", np.round(np.mean([r["utv"] for r in rs], 0)), "exc", round(max(r["excess"] for r in rs), 3), "finite", all(r["finite"] for r in rs), flush=True)


if __name__ == "__main__":
    main()

"""Test 3: realistic digital implementation. Rules fixed before any result was seen.

Digital control at 1 kHz with zero-order hold and a one-sample delay. The plant is integrated
at 10 kHz. Observer, fuzzy update and command filter use forward Euler at the sample rate.

Sensor conditions (12-bit unless stated)
  S1 quantization and sampling only
  S2 S1 + Gaussian noise of 0.5 count before quantization
  S3 S2 + constant bias within +-2 counts + drift within +-0.5 count/s

Designs: base (full cancellation), low-pass, dead zone 1 count, dead zone 2 counts, linear.
Scenario seeds: 7 (development), 8 and 9 (test). 100 runs per job.

  R1 robustness: under S2 at gains 40 AND 57, a dead zone (1 or 2 counts, the same one at
     both gains) gives <= 1% violating runs and < 10% saturated time, pooled over 3 seeds,
     with tracking RMS no more than 2x the best design at that gain.
  R2 value of fixed-time terms: choose each design's gain on seed 7 under S2: the gain with the
     lowest velocity-estimate RMS among gains with <= 1% violating runs and < 10% saturation.
     Evaluate on seeds 8 and 9. Pass if the best dead-zone design has velocity-estimate RMS at
     least 20% below the linear observer AND tracking RMS no worse than linear.
     If R2 fails, the practical recommendation is a linear observer, and the fixed-time
     framing has no practical support from these tests.
  R3 small motions: for slow and hold scenarios under S2 at gains 28, 40 and 57, the chosen dead zone
     has tracking RMS <= 1.5x the linear observer at the same gain, and peak-to-peak hold
     error in the last 2 s <= 4 counts.
  R4 gain law: base design, S2-type noise in counts, 10/12/14/16 bits, 3 seeds, finer gains.
     Ceiling = lowest gain with pooled median saturated share > 10%. Fitted slope within 0.10 of -0.25.
"""
import json
import multiprocessing as mp
import sys
import time

import numpy as np

import arm_margin_test as A
from arm_margin_test import B_NOM, C_FILT, CW, CY, GAMMA, F_MAX, K1, K1B, K2, K2B, P, Q, SW, SY, TH_L, TH_U, U_MAX, W_L, W_U, D_SIN, h, h_x, sig
from run_all import clean

TS = 1e-3
SUB = 10
T_END = 4.0
RUNS = 100
CLEAR = 0.01                    # command clearance from the limits, as a share of their width
GAINS = (14.0, 20.0, 24.0, 28.0, 34.0, 40.0, 48.0, 57.0, 68.0, 80.0, 113.0)
DESIGNS = {"base": ("frac", "full", 0.0), "low-pass": ("frac", "lp", 0.0), "dead zone 1": ("frac", "full", 1.0),
           "dead zone 2": ("frac", "full", 2.0), "linear": ("linear", "full", 0.0)}
SENSORS = {"S1": (0.0, 0.0, 0.0), "S2": (0.5, 0.0, 0.0), "S3": (0.5, 2.0, 0.5)}


def setup(sc, kind):
    n = sc["J"].size
    r = np.random.default_rng(1234)
    if kind == "slow":
        sc["w"] = r.uniform(0.3, 0.6, n)
        sc["amp"] = r.uniform(0.02, 0.05, n) / sc["w"]
        sc["pa"] = np.zeros(n)
    elif kind == "hold":
        sc["amp"] = np.zeros(n)
        sc["pa"] = np.zeros(n)
        sc["w0"] = np.zeros(n)
        sc["th0"] = sc["off"].copy()
    return sc


def simulate(sc, lam, obs, cancel, dz_c, bits, sensor, seed, trace=None, idx=0):
    n = sc["J"].size
    q = 2 * np.pi / 2 ** bits
    ns, bs, dr = sensor
    a1, b1 = (1.0, 1.0) if obs == "linear" else (A.A1, A.B1)
    a2, b2 = (1.0, 1.0) if obs == "linear" else (A.A2, A.B2)
    l1, l1b, l2, l2b = A.obs_gains(lam)
    rng = np.random.default_rng(seed)
    bias = rng.uniform(-bs, bs, n) * q
    drift = rng.uniform(-dr, dr, n) * q
    th, w = sc["th0"].copy(), sc["w0"].copy()
    xh1, xh2 = th.copy(), np.zeros(n)
    theta = np.zeros((n, 25))
    vh, i2f, u_app, u_last = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
    steps = int(T_END / TS)
    o1h, o1l, oh, ol = TH_U, -TH_L, W_U, -W_L
    rho = CLEAR * (oh + ol)
    lo, hi = -ol + 1e-3, oh - 1e-3
    wv, tv, xv = np.zeros(n), np.zeros(n), np.zeros(n)
    sq, e2s, cnt, cnt2 = np.zeros(n), np.zeros(n), 0, 0
    tvu, sat, upk = np.zeros(n), np.zeros(n), np.zeros(n)
    mx, mn = np.full(n, -np.inf), np.full(n, np.inf)
    a_lp = TS * 2 * np.pi * 30.0
    dts = TS / SUB
    for k in range(steps):
        t = k * TS
        yd = sc["off"] + sc["amp"] * np.sin(sc["w"] * t + sc["ph"])
        dyd = sc["amp"] * sc["w"] * np.cos(sc["w"] * t + sc["ph"])
        y = th + bias + drift * t + ns * q * rng.standard_normal(n)
        y = np.round(y / q) * q
        yc = np.clip(y, TH_L + 1e-3, TH_U - 1e-3)
        e1 = y - xh1
        if dz_c > 0:
            e1 = np.sign(e1) * np.maximum(np.abs(e1) - dz_c * q, 0.0)
        i1 = l1 * sig(e1, a1) + l1b * sig(e1, b1)
        i2 = l2 * sig(e1, a2) + l2b * sig(e1, b2)
        psi = np.exp(-((yc[:, None, None] - CY[None, :, None]) / SY) ** 2 - ((xh2[:, None, None] - CW[None, None, :]) / SW) ** 2).reshape(n, 25)
        psi /= psi.sum(1, keepdims=True) + 1e-9
        fh = (theta * psi).sum(1)
        z1 = h(yc, o1h, o1l) - h(yd, o1h, o1l)
        v = (h_x(yd, o1h, o1l) * dyd - K1 * sig(z1, P) - K1B * sig(z1, Q)) / h_x(yc, o1h, o1l)
        vc = np.clip(v, -ol + rho, oh - rho)
        dvh = C_FILT * (vc - vh)
        vhc, xc = np.clip(vh, lo, hi), np.clip(xh2, lo, hi)
        alpha = h(vhc, oh, ol)
        dalpha = h_x(vhc, oh, ol) * dvh
        z2 = h(xc, oh, ol) - alpha
        ds2 = dalpha - K2 * sig(z2, P) - K2B * sig(z2, Q)
        dx2 = ds2 / h_x(xc, oh, ol)
        i2f = i2f + a_lp * (i2 - i2f)
        i2u = {"full": i2, "none": 0.0, "lp": i2f}[cancel]
        u_raw = (dx2 - fh - i2u) / B_NOM
        u_new = np.clip(u_raw, -U_MAX, U_MAX)
        sat += np.abs(u_raw) >= U_MAX
        tvu += np.abs(u_new - u_last)
        upk = np.maximum(upk, np.abs(u_new))
        u_last = u_new
        # plant over one sample with the delayed torque
        for j in range(SUB):
            tj = t + j * dts
            pulse = np.where((tj >= sc["pt"]) & (tj < sc["pt"] + sc["pd"]), sc["pa"], 0.0)
            d = D_SIN * np.sin(2 * tj) + pulse
            dw = (u_app - sc["B"] * w - sc["G"] * np.sin(th) - sc["FC"] * np.tanh(20 * w) + d) / sc["J"]
            th = th + dts * w
            w = w + dts * dw
            wv = np.maximum(wv, np.maximum(w - W_U, W_L - w))
            tv = np.maximum(tv, np.maximum(th - TH_U, TH_L - th))
        xv = np.maximum(xv, np.maximum(xh2 - oh, -ol - xh2))
        if trace is not None:
            trace.append((t, float(th[idx]), float(yd[idx]), float(w[idx]), float(xh2[idx]), float(u_new[idx])))
        if t >= 2.0:
            err = th - yd
            sq += err * err
            cnt += 1
            mx, mn = np.maximum(mx, err), np.minimum(mn, err)
        if t >= 1.0:
            e2s += (w - xh2) ** 2
            cnt2 += 1
        # discrete observer, fuzzy update and filter, using the torque applied in this sample
        theta = np.clip(theta + TS * GAMMA * psi * i2[:, None], -F_MAX, F_MAX)
        xh1 = xh1 + TS * (xh2 + i1)
        xh2 = xh2 + TS * (B_NOM * u_app + fh + i2)
        vh = vh + TS * dvh
        u_app = u_new
    return dict(viol=int((wv > 0).sum()), est=int((xv > 0).sum()), ang=int((tv > 0).sum()), excess=float(wv.max()), runs=n,
                sat=float(np.median(sat / steps)), rms=float(np.median(np.sqrt(sq / cnt))), rms95=float(np.percentile(np.sqrt(sq / cnt), 95)),
                e2=float(np.median(np.sqrt(e2s / cnt2))), p2p=float(np.median(mx - mn)), utv=float(np.median(tvu / T_END)), upk=float(np.median(upk)), finite=bool(np.isfinite(th).all()))


def job(j):
    with np.errstate(all="ignore"):
        sc = setup(A.scenarios(RUNS, seed=j["seed"]), j["scen"])
        obs, cancel, dz = DESIGNS[j["design"]]
        r = simulate(sc, j["lam"], obs, cancel, dz, j["bits"], SENSORS[j["sensor"]], 100 + j["seed"])
    r.update({k: j[k] for k in ("part", "design", "lam", "bits", "sensor", "seed", "scen")})
    return r


def main():
    jobs = []
    add = lambda **k: jobs.append(k)
    for d in DESIGNS:
        for lam in GAINS:
            add(part="S1", design=d, lam=lam, bits=12, sensor="S1", seed=7, scen="std")
            for s in (7, 8, 9):
                add(part="S2", design=d, lam=lam, bits=12, sensor="S2", seed=s, scen="std")
        for lam in (28.0, 40.0, 57.0):
            for s in (7, 8, 9):
                add(part="S3", design=d, lam=lam, bits=12, sensor="S3", seed=s, scen="std")
                add(part="small", design=d, lam=lam, bits=12, sensor="S2", seed=s, scen="slow")
                add(part="small", design=d, lam=lam, bits=12, sensor="S2", seed=s, scen="hold")
    for b in (10, 14, 16):
        for lam in GAINS:
            for s in (7, 8, 9):
                add(part="law", design="base", lam=lam, bits=b, sensor="S2", seed=s, scen="std")
    print(len(jobs), "jobs", flush=True)
    t0 = time.time()
    rows = []
    with mp.Pool(4) as p:
        for i, r in enumerate(p.imap_unordered(job, jobs, chunksize=2)):
            rows.append(r)
            if i % 50 == 0:
                print(i, round(time.time() - t0), flush=True)
    (A.OUT / "test3.json").write_text(json.dumps(clean(rows), indent=0))
    print("saved", round(time.time() - t0), flush=True)


if __name__ == "__main__":
    main()

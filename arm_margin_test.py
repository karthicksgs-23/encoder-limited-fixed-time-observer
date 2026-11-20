"""First feasibility test for the observer-aware margin.

Single-link arm, only the joint angle is measured. The same controller and observer
are run with different margins on the velocity limits. All cases see identical
randomised scenarios.

Criteria fixed before any result was seen:
  C1 motive : with no margin, true velocity breaks its limit in >= 1% of runs.
  C2 safety : with the time-varying margin, 0 true-state violations in all runs.
  C3 cost   : with the time-varying margin, tracking RMS (t >= 2 s) is within 25% of
              the no-margin case, peak torque rises by no more than 10%, and the
              margin fits inside the limits (m(0) < smallest limit distance).
  C4 value  : the time-varying margin tracks better than a constant margin of equal
              safety. If it does not, a constant back-off is enough and the time-varying
              margin adds nothing.
"""
import argparse
import json
import pathlib
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

# Nominal plant and uncertainty ranges
J0, B0, G0, FC0 = 0.5, 0.4, 2.5, 0.2
SPREAD = dict(J=0.2, B=0.3, G=0.2, FC=0.5)
U_MAX = 8.0
TH_L, TH_U = -1.0, 1.4          # angle limits (rad), measured
W_L, W_U = -0.9, 1.3            # velocity limits (rad/s), not measured
E0 = 0.5                        # known bound on initial velocity error
PULSE_MAX = 3.0                 # torque pulse bound (N m)
D_SIN = 0.2                     # 0.2 sin(2t) disturbance
B_NOM = 1.0 / J0

# Controller
P, Q = 0.8, 1.2
K1, K1B = 4.0, 1.0
K2, K2B = 25.0, 2.0
C_FILT = 40.0
RHO = 0.05                      # command clearance from the limits, as a share of their width
# Observer exponents, a2 = 2 a1 - 1 and b2 = 2 b1 - 1
A1, B1 = 0.75, 1.25
A2, B2 = 2 * A1 - 1, 2 * B1 - 1
# Fuzzy system: 5 x 5 Gaussian rules on (y, velocity estimate)
CY = np.linspace(TH_L, TH_U, 5)
CW = np.linspace(W_L, W_U, 5)
SY, SW = 0.35, 0.33
GAMMA = 2.0
F_MAX = 6.0                     # projection box on each rule consequent


def sig(x, p):
    return np.sign(x) * np.abs(x) ** p


def h(x, oh, ol):
    return oh * ol * x / ((oh - x) * (ol + x))


def h_x(x, oh, ol):
    return oh * ol * (oh * ol + x * x) / ((oh - x) * (ol + x)) ** 2


def h_t(x, oh, ol, doh, dol):
    return (ol * x / (ol + x)) * (-x / (oh - x) ** 2) * doh + (oh * x / (oh - x)) * (x / (ol + x) ** 2) * dol


def obs_gains(lam):
    return 2.0 * lam, 2.0 * lam, lam * lam, lam * lam


def error_envelope(lam, D, dt=1e-4, T=6.0):
    """Sampled worst case of |e2| for the observer error system under |delta| <= D.

    This is an empirical envelope from adversarial bang-bang inputs, not a proof.
    """
    l1, l1b, l2, l2b = obs_gains(lam)
    n = int(T / dt)
    tt = np.arange(n) * dt
    periods = np.array([np.inf, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0])
    phases = np.linspace(0, 1, 5, endpoint=False)
    e20 = np.array([-E0, 0.0, E0])
    Pm, Ph, Ee = np.meshgrid(periods, phases, e20, indexing="ij")
    Pm, Ph, Ee = Pm.ravel(), Ph.ravel(), Ee.ravel()
    sgn0 = np.array([1.0, -1.0])
    Pm, Ph, Ee = np.tile(Pm, 2), np.tile(Ph, 2), np.tile(Ee, 2)
    s0 = np.repeat(sgn0, Pm.size // 2)
    e1 = np.zeros_like(Ee)
    e2 = Ee.copy()
    env = np.zeros(n)
    for k in range(n):
        env[k] = np.abs(e2).max()
        ph = np.where(np.isinf(Pm), 0.0, tt[k] / np.where(np.isinf(Pm), 1.0, Pm) + Ph)
        delta = D * s0 * np.where((ph % 1.0) < 0.5, 1.0, -1.0)
        de1 = e2 - l1 * sig(e1, A1) - l1b * sig(e1, B1)
        de2 = delta - l2 * sig(e1, A2) - l2b * sig(e1, B2)
        e1 = e1 + dt * de1
        e2 = e2 + dt * de2
    return tt, env


def margin_from_envelope(tt, env, inflate=1.1):
    """Non-increasing smooth margin above the sampled envelope."""
    top = np.maximum.accumulate(env[::-1])[::-1] * inflate
    m_ss = top[-1]
    m0 = top[0]
    idx = np.argmax(top <= m_ss * 1.02 + 1e-12)
    t_o = max(tt[idx], 1e-3)
    return dict(m0=float(m0), m_ss=float(m_ss), t_o=float(t_o))


def margin_fn(kind, mp):
    """Returns m(t) and its derivative for scalar t."""
    if kind == "none":
        return lambda t: (0.0, 0.0)
    if kind == "constant":
        return lambda t: (mp["m0"], 0.0)
    m0, ms, to = mp["m0"], mp["m_ss"], mp["t_o"]

    def f(t):
        if t >= to:
            return ms, 0.0
        r = t / to
        return m0 + (ms - m0) * (3 * r * r - 2 * r ** 3), (ms - m0) * (6 * r - 6 * r * r) / to
    return f


def scenarios(n, seed=7):
    r = np.random.default_rng(seed)
    u = lambda s: r.uniform(1 - s, 1 + s, n)
    a_w = r.uniform(0.5, 0.8, n)        # peak reference speed
    w = r.uniform(1.2, 2.0, n)
    return dict(
        J=J0 * u(SPREAD["J"]), B=B0 * u(SPREAD["B"]), G=G0 * u(SPREAD["G"]), FC=FC0 * u(SPREAD["FC"]),
        th0=r.uniform(-0.2, 0.2, n), w0=r.uniform(-E0, E0, n),
        amp=a_w / w, w=w, ph=r.uniform(0, 2 * np.pi, n), off=r.uniform(-0.1, 0.3, n),
        pa=r.uniform(-PULSE_MAX, PULSE_MAX, n), pt=r.uniform(1.0, 5.0, n), pd=r.uniform(0.05, 0.2, n),
        noise=r.standard_normal,
    )


def simulate(sc, margin, lam, T=6.0, dt=1e-4, noise_std=0.0, quant=0.0, keep=None, seed=11, cancel="full", dz=0.0, lp_hz=30.0):
    """cancel: how the observer injection enters u ('full', 'none' or 'lp' for low-pass).
    dz: dead zone on the output error, in rad."""
    n = sc["J"].size
    i2f = np.zeros(n)
    l1, l1b, l2, l2b = obs_gains(lam)
    rng = np.random.default_rng(seed)
    th, w = sc["th0"].copy(), sc["w0"].copy()
    xh1, xh2 = th.copy(), np.zeros(n)
    theta = np.zeros((n, 25))
    vh = np.zeros(n)
    yd0 = sc["off"] + sc["amp"] * np.sin(sc["ph"])
    steps = int(T / dt)
    o1h, o1l = TH_U, -TH_L
    stat = dict(wv=np.zeros(n), tv=np.zeros(n), xv=np.zeros(n), first=np.full(n, np.nan), e2=np.zeros(n), env_ok=np.ones(n, bool),
                t_e2=np.full(n, np.nan), e2rms=np.zeros(n), e2cnt=0, sq=np.zeros(n), cnt=0, upk=np.zeros(n), usq=np.zeros(n), sat=np.zeros(n), dmax=np.zeros(n))
    tr = {k: [] for k in ("t", "th", "w", "xh2", "yd", "u", "m")} if keep is not None else None
    for k in range(steps):
        t = k * dt
        m, dm = margin(t)
        oh, ol = W_U - m, -W_L - m
        yd = sc["off"] + sc["amp"] * np.sin(sc["w"] * t + sc["ph"])
        dyd = sc["amp"] * sc["w"] * np.cos(sc["w"] * t + sc["ph"])
        y = th
        if noise_std > 0:
            y = y + noise_std * rng.standard_normal(n)
        if quant > 0:
            y = np.round(y / quant) * quant
        yc = np.clip(y, TH_L + 1e-3, TH_U - 1e-3)
        e1 = y - xh1
        if dz > 0:
            e1 = np.sign(e1) * np.maximum(np.abs(e1) - dz, 0.0)
        i1 = l1 * sig(e1, A1) + l1b * sig(e1, B1)
        i2 = l2 * sig(e1, A2) + l2b * sig(e1, B2)
        psi = np.exp(-((yc[:, None, None] - CY[None, :, None]) / SY) ** 2 - ((xh2[:, None, None] - CW[None, None, :]) / SW) ** 2).reshape(n, 25)
        psi /= psi.sum(1, keepdims=True) + 1e-9
        fh = (theta * psi).sum(1)
        z1 = h(yc, o1h, o1l) - h(yd, o1h, o1l)
        v = (h_x(yd, o1h, o1l) * dyd - K1 * sig(z1, P) - K1B * sig(z1, Q)) / h_x(yc, o1h, o1l)
        rho = RHO * (oh + ol)
        vc = np.clip(v, -ol + rho, oh - rho)
        dvh = C_FILT * (vc - vh)
        lo, hi = -ol + 1e-3, oh - 1e-3
        vhc = np.clip(vh, lo, hi)
        xc = np.clip(xh2, lo, hi)
        alpha = h(vhc, oh, ol)
        dalpha = h_x(vhc, oh, ol) * dvh + h_t(vhc, oh, ol, -dm, -dm)
        z2 = h(xc, oh, ol) - alpha
        ds2 = dalpha - K2 * sig(z2, P) - K2B * sig(z2, Q)
        dx2 = (ds2 - h_t(xc, oh, ol, -dm, -dm)) / h_x(xc, oh, ol)
        i2f = i2f + dt * 2 * np.pi * lp_hz * (i2 - i2f)
        i2u = {"full": i2, "none": 0.0, "lp": i2f}[cancel]
        u_raw = (dx2 - fh - i2u) / B_NOM
        u = np.clip(u_raw, -U_MAX, U_MAX)
        pulse = np.where((t >= sc["pt"]) & (t < sc["pt"] + sc["pd"]), sc["pa"], 0.0)
        d = D_SIN * np.sin(2 * t) + pulse
        f_true = (-sc["B"] * w - sc["G"] * np.sin(th) - sc["FC"] * np.tanh(20 * w) + d) / sc["J"]
        dw = u / sc["J"] + f_true
        delta = (1 / sc["J"] - B_NOM) * u + f_true - fh

        # metrics use true states and true limits
        over = np.maximum(w - W_U, W_L - w)
        new = (over > 0) & np.isnan(stat["first"])
        stat["first"][new] = t
        stat["wv"] = np.maximum(stat["wv"], over)
        stat["tv"] = np.maximum(stat["tv"], np.maximum(th - TH_U, TH_L - th))
        stat["xv"] = np.maximum(stat["xv"], np.maximum(xh2 - oh, -ol - xh2))
        e2 = np.abs(w - xh2)
        stat["e2"] = np.maximum(stat["e2"], e2)
        hit = (e2 < 0.05) & np.isnan(stat["t_e2"])
        stat["t_e2"][hit] = t
        if t >= 1.0:
            stat["e2rms"] += e2 * e2
            stat["e2cnt"] += 1
        stat["env_ok"] &= e2 <= m + 1e-9
        stat["dmax"] = np.maximum(stat["dmax"], np.abs(delta))
        stat["upk"] = np.maximum(stat["upk"], np.abs(u))
        stat["usq"] += u * u * dt
        stat["sat"] += (np.abs(u_raw) >= U_MAX) * dt
        if t >= 2.0:
            stat["sq"] += (th - yd) ** 2
            stat["cnt"] += 1
        if tr is not None and k % 50 == 0:
            for a, b in zip(("t", "th", "w", "xh2", "yd", "u", "m"), (t, th[keep], w[keep], xh2[keep], yd[keep], u[keep], m)):
                tr[a].append(float(b))

        # integrate
        dth1 = xh2 + i1
        dth2 = B_NOM * u + fh + i2
        theta = np.clip(theta + dt * GAMMA * psi * i2[:, None], -F_MAX, F_MAX)
        xh1 = xh1 + dt * dth1
        xh2 = xh2 + dt * dth2
        vh = vh + dt * dvh
        th = th + dt * w
        w = w + dt * dw
    stat["rms"] = np.sqrt(stat["sq"] / max(stat["cnt"], 1))
    return stat, tr


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    s = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(0.0, c - s), min(1.0, c + s)


def summary(name, s, n):
    k = int((s["wv"] > 0).sum())
    lo, hi = wilson(k, n)
    v = s["wv"][s["wv"] > 0]
    return dict(case=name, runs=n, violating_runs=k, pct=100 * k / n, ci_lo=100 * lo, ci_hi=100 * hi,
                worst_excess=float(s["wv"].max()), median_excess=float(np.median(v)) if k else 0.0,
                median_first_t=float(np.nanmedian(s["first"])) if k else None,
                late_violations=int((s["first"] > 1.0).sum()),
                angle_violating_runs=int((s["tv"] > 0).sum()), est_violating_runs=int((s["xv"] > 0).sum()),
                margin_breaches=int((~s["env_ok"]).sum()),
                rms_median=float(np.median(s["rms"])), rms_p95=float(np.percentile(s["rms"], 95)),
                peak_u_median=float(np.median(s["upk"])), peak_u_max=float(s["upk"].max()),
                energy_median=float(np.median(s["usq"])), sat_time_median=float(np.median(s["sat"])),
                max_e2=float(s["e2"].max()), max_delta=float(s["dmax"].max()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", type=int, default=100)
    ap.add_argument("--lam", type=float, nargs="+", default=[60.0])
    ap.add_argument("--tag", default="pilot")
    ap.add_argument("--dmargin", type=float, default=None, help="delta bound used for the margin; default is the guaranteed bound")
    ap.add_argument("--noise", type=float, default=0.0)
    ap.add_argument("--quant", type=float, default=0.0)
    ap.add_argument("--dt", type=float, default=1e-4)
    a = ap.parse_args()
    # guaranteed bound: gain error x torque + worst f + worst fuzzy output
    Jlo = J0 * (1 - SPREAD["J"])
    d_full = (1 / Jlo - B_NOM) * U_MAX + (B0 * 1.3 * max(W_U, -W_L) * 1.5 + G0 * 1.2 + FC0 * 1.5 + D_SIN + PULSE_MAX) / Jlo + F_MAX
    sc = scenarios(a.runs)
    rows, meta = [], dict(d_guaranteed=d_full, args=vars(a))
    for lam in a.lam:
        D = a.dmargin if a.dmargin else d_full
        tt, env = error_envelope(lam, D)
        mp = margin_from_envelope(tt, env)
        room = min(W_U, -W_L)
        print(f"lam={lam} D={D:.1f} margin {mp}, room {room}", flush=True)
        meta[f"margin_lam{lam:g}"] = mp
        cases = [("no margin", "none")]
        if mp["m0"] < 0.9 * room:
            cases += [("constant margin", "constant"), ("time-varying margin", "tv")]
        else:
            print("  margin does not fit inside the limits; margin cases skipped", flush=True)
        for name, kind in cases:
            t0 = time.time()
            s, _ = simulate(sc, margin_fn(kind, mp), lam, dt=a.dt, noise_std=a.noise, quant=a.quant)
            r = summary(name, s, a.runs)
            r.update(lam=lam, D=D, **mp)
            rows.append(r)
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k not in ("D", "m0", "m_ss", "t_o", "lam")}, f"{time.time()-t0:.0f}s", flush=True)
    (OUT / f"{a.tag}.json").write_text(json.dumps(dict(meta=meta, rows=rows), indent=1))


if __name__ == "__main__":
    main()

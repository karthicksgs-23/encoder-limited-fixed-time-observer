"""Test 6: does the observer error settle within a bound that stops growing with the initial error?

Undisturbed error system with lam = 1, then lam = 20 to check the scaling T(lam) = T1/lam.
Settling = first time max(|e1|, |e2|/lam) < 1e-3, after which the error stayed below it to the
end of the integration. A sampled check, not a proof of fixed-time stability.
"""
import json

import numpy as np
from scipy.integrate import solve_ivp

A1, B1 = 0.75, 1.25
A2, B2 = 0.5, 1.5
sig = lambda e, p: np.sign(e) * np.abs(e) ** p


def f(t, x, lam):
    e1, e2 = x
    return [e2 - 2 * lam * (sig(e1, A1) + sig(e1, B1)), -lam ** 2 * (sig(e1, A2) + sig(e1, B2))]


def settle(e10, e20, lam, T=60.0):
    ev = lambda t, x, l: max(abs(x[0]), abs(x[1]) / l) - 1e-3
    ev.terminal = True
    ev.direction = -1
    s = solve_ivp(f, (0, T / lam), [e10, e20], args=(lam,), method="LSODA", rtol=1e-7, atol=1e-10, events=ev)
    return float(s.t_events[0][0]) if s.t_events[0].size else None


rows = []
for lam in (1.0, 20.0):
    for mag in (1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6):
        ts = []
        for th in np.linspace(0, 2 * np.pi, 16, endpoint=False):
            ts.append(settle(mag * np.cos(th), lam * mag * np.sin(th), lam))
        ok = [t for t in ts if t is not None]
        rows.append(dict(lam=lam, mag=mag, n=len(ts), settled=len(ok), tmax=max(ok) if ok else None))
        print(lam, f"{mag:g}", len(ok), "/", len(ts), None if not ok else round(max(ok), 4), flush=True)
json.dump(rows, open("results/test6.json", "w"), indent=0)

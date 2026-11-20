"""Figures, all drawn from saved simulation results."""
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

R = pathlib.Path(__file__).parent / "results"
F = pathlib.Path(__file__).parent / "figures"
F.mkdir(exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman"], "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#E3E3E3", "grid.linewidth": .6, "lines.linewidth": 1.5, "savefig.dpi": 600, "savefig.bbox": "tight"})
C = dict(prop="#0B6E4F", base="#C8501E", lin="#2A5CAA", dz="#7A4FA3", ink="#222222", lim="#8A8A8A")


@matplotlib.rc_context({"mathtext.fontset": "stix"})
def diagram():
    fig, ax = plt.subplots(figsize=(6.9, 3.9))
    ax.set_xlim(0, 106)
    ax.set_ylim(0, 60)
    ax.axis("off")
    ax.grid(False)
    ink = "#333333"

    def box(x, y, w, h, t, c):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.2,rounding_size=1", fc=c, ec=ink, lw=.9, zorder=3))
        ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=8.6, zorder=4)

    def line(pts, head=True):
        """Polyline with an arrow head on its last segment."""
        for a_, b_ in zip(pts[:-2], pts[1:-1]):
            ax.plot([a_[0], b_[0]], [a_[1], b_[1]], color=ink, lw=.9, solid_capstyle="round", zorder=2)
        ax.add_patch(FancyArrowPatch(pts[-2], pts[-1], arrowstyle="-|>" if head else "-", mutation_scale=9, lw=.9, color=ink, shrinkA=0, shrinkB=0, zorder=2))

    def lab(x, y, t, ha="center", va="bottom"):
        ax.text(x, y, t, ha=ha, va=va, fontsize=8.6, zorder=4, bbox=dict(fc="white", ec="none", pad=.4))

    def dot(x, y):
        ax.plot([x], [y], "o", ms=3, color=ink, zorder=3)

    box(13, 42, 15, 10, "Limit\nmapping", "#EAF2FB")
    box(35, 42, 17, 10, "Backstepping\ncontrol law", "#EAF2FB")
    box(59, 42, 10, 10, "Torque\nlimit", "#FDF0EA")
    box(76, 42, 11, 10, "Robot\narm", "#F1F1F1")
    box(93, 42, 12, 10, "Encoder\n(step $q$)", "#FDF0EA")
    box(64, 6, 18, 10, "Fixed-time\nobserver", "#E8F7F1")
    box(38, 6, 17, 10, "Fuzzy\napproximator", "#E8F7F1")
    box(13, 6, 15, 10, "Low-pass\nfilter", "#E8F7F1")

    # Input and forward path
    line([(1, 47), (12.8, 47)])
    lab(6.5, 47.7, "$y_d,\\ \\dot{y}_d$\nlimits")
    line([(28.2, 47), (34.8, 47)])
    lab(31.5, 47.7, "$z_1, z_2$")
    line([(52.2, 47), (58.8, 47)])
    lab(55.5, 47.7, "$u$")
    line([(69.2, 47), (75.8, 47)])
    lab(72.5, 47.7, "sat$(u)$")
    line([(87.2, 47), (92.8, 47)])
    lab(90, 47.7, "$x_1$")
    # External effects
    line([(81.5, 58), (81.5, 52.4)])
    lab(80.5, 54.2, "$d(t)$", ha="right", va="bottom")
    line([(99, 58), (99, 52.4)])
    lab(98, 54.2, "$n,\\ \\beta,\\ \\gamma t$", ha="right", va="bottom")
    # Measured output to the observer, applied torque to the observer
    line([(99, 41.8), (99, 11), (82.2, 11)])
    lab(90.5, 11.6, "$y_q$")
    line([(72.5, 47), (72.5, 16.2)])
    dot(72.5, 47)
    # Feedback bus: y_q and the velocity estimate
    line([(99, 30), (20.5, 30), (20.5, 41.8)])
    dot(99, 30)
    line([(78, 16.2), (78, 30)], head=False)
    dot(78, 30)
    lab(79, 21, "$\\hat{x}_2$", ha="left", va="center")
    lab(88, 30.6, "$y_q$")
    lab(60, 30.6, "$y_q,\\ \\hat{x}_2$")
    line([(48, 30), (48, 41.8)])
    line([(48, 30), (48, 16.2)])
    dot(48, 30)
    # Fuzzy output
    line([(42.5, 16.2), (42.5, 41.8)])
    lab(43.3, 36, "$\\hat{f}$", ha="left", va="center")
    line([(55.2, 11), (63.8, 11)])
    lab(59.5, 11.6, "$\\hat{f}$")
    # Injection goes directly from the observer to the fuzzy update and the filter
    line([(73, 5.8), (73, 2), (20.5, 2), (20.5, 5.8)])
    line([(46.5, 2), (46.5, 5.8)])
    dot(46.5, 2)
    lab(60, 2.5, "$\\iota$")
    # Filtered injection
    line([(20.5, 16.2), (20.5, 23), (37.5, 23), (37.5, 41.8)])
    lab(29, 23.6, "$\\iota_f$")
    fig.savefig(F / "fig1_block_diagram.png")
    plt.close(fig)


def responses():
    T = json.load(open(R / "traces.json"))
    fig, ax = plt.subplots(4, 1, figsize=(6.5, 6.6), sharex=True)
    S = [("low-pass 20", C["prop"], "Proposed, gain 20"), ("base 40", C["base"], "Full cancellation, gain 40")]
    a = np.array(T["low-pass 20"])
    ax[0].plot(a[:, 0], a[:, 2], color=C["ink"], ls="--", lw=1.1, label="Reference")
    for k, c, l in S:
        a = np.array(T[k])
        ax[0].plot(a[:, 0], a[:, 1], color=c, label=l)
        ax[2 if k == "low-pass 20" else 3].plot(a[:, 0], a[:, 5], color=c, lw=.5 if k == "base 40" else 1.0, label=l)
    a = np.array(T["low-pass 20"])
    ax[1].plot(a[:, 0], a[:, 3], color=C["prop"], label="True velocity")
    ax[1].plot(a[:, 0], a[:, 4], color=C["lin"], lw=.9, label="Estimated velocity")
    for y in (-0.9, 1.3):
        ax[1].axhline(y, color=C["lim"], ls=":", lw=1.2)
    ax[1].text(3.98, 1.22, "Velocity limits", ha="right", va="top", fontsize=8.5, color="#555555")
    ax[0].set_ylabel("Joint angle (rad)")
    ax[1].set_ylabel("Velocity (rad/s)")
    ax[2].set_ylabel("Torque (N·m)")
    ax[3].set_ylabel("Torque (N·m)")
    ax[3].set_xlabel("Time (s)")
    ax[2].set_ylim(-8.8, 8.8)
    ax[3].set_ylim(-8.8, 8.8)
    for x in ax:
        x.legend(frameon=False, fontsize=8.5, ncol=3, loc="upper center", bbox_to_anchor=(.5, 1.2))
        p = T["pulse"]
        x.axvspan(p[0], p[0] + p[1], color="#DDDDDD", zorder=0, lw=0)
    ax[0].set_xlim(0, 4)
    fig.tight_layout()
    fig.savefig(F / "fig2_time_responses.png")
    plt.close(fig)


def ceilings():
    q = lambda b: 2 * np.pi / 2 ** b
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.7), sharey=True)
    b = np.linspace(9.5, 16.5, 50)
    ax[0].plot(b, np.sqrt(16 / q(b) ** .5), color=C["prop"], label="Gain rule")
    ax[0].scatter([10, 12, 14, 16], [14, 20, 28, 40], color=C["base"], zorder=3, s=26, label="Measured")
    ax[0].set_xlabel("Encoder resolution (bits)")
    ax[0].set_ylabel("Highest usable gain")
    u = np.linspace(3, 18, 50)
    ax[1].plot(u, np.sqrt(2 * u / q(12) ** .5), color=C["prop"], label="Gain rule")
    ax[1].scatter([4, 8, 16], [14, 20, 28], color=C["base"], zorder=3, s=26, label="Measured")
    ax[1].set_xlabel("Torque limit (N·m)")
    ax[0].legend(frameon=False, fontsize=9)
    ax[0].set_ylim(0, 48)
    fig.tight_layout()
    fig.savefig(F / "fig3_gain_ceiling.png")
    plt.close(fig)


def saturation():
    Rw = json.load(open(R / "test3.json"))
    fig, ax = plt.subplots(figsize=(5.6, 2.9))
    for d, c, l, m in (("base", C["base"], "Full cancellation", "o"), ("low-pass", C["prop"], "Low-pass", "s"), ("dead zone 1", C["dz"], "Dead zone", "^"), ("linear", C["lin"], "Linear observer", "D")):
        g = sorted({r["lam"] for r in Rw if r["part"] == "S2" and r["design"] == d})
        y = [100 * float(np.median([r["sat"] for r in Rw if r["part"] == "S2" and r["design"] == d and r["lam"] == x])) for x in g]
        ax.plot(g, y, color=c, marker=m, ms=4, label=l)
    ax.axhline(10, color=C["lim"], ls=":", lw=1.2)
    ax.text(113, 12, "10% threshold", ha="right", fontsize=8.5, color="#555555")
    ax.set_xscale("log")
    ax.set_xticks([14, 20, 28, 40, 57, 80, 113])
    ax.set_xticklabels(["14", "20", "28", "40", "57", "80", "113"])
    ax.minorticks_off()
    ax.set_xlabel("Observer gain")
    ax.set_ylabel("Time in saturation (%)")
    ax.legend(frameon=False, fontsize=8.5, ncol=2, loc="upper left")
    ax.set_ylim(0, 118)
    fig.savefig(F / "fig4_saturation.png")
    plt.close(fig)


def two_link():
    T = json.load(open(R / "traces7.json"))
    a = np.array(T["proposed"])
    fig, ax = plt.subplots(3, 1, figsize=(6.5, 5.0), sharex=True)
    for j, c, n in ((0, C["prop"], "Joint 1"), (1, C["lin"], "Joint 2")):
        ax[0].plot(a[:, 0], a[:, 1 + j], color=c, label=n)
        ax[0].plot(a[:, 0], a[:, 3 + j], color=c, ls="--", lw=1, label=n + " reference")
        ax[1].plot(a[:, 0], a[:, 5 + j], color=c, label=n)
        ax[2].plot(a[:, 0], a[:, 9 + j], color=c, lw=.8, label=n)
    ax[0].set_ylabel("Angle (rad)")
    ax[1].set_ylabel("Velocity (rad/s)")
    ax[2].set_ylabel("Torque (N·m)")
    ax[2].set_xlabel("Time (s)")
    ax[0].legend(frameon=False, fontsize=8.5, ncol=4, loc="upper center", bbox_to_anchor=(.5, 1.22))
    for x in ax:
        x.axvspan(T["pulse"][0], T["pulse"][0] + T["pulse"][1], color="#DDDDDD", zorder=0, lw=0)
    ax[0].set_xlim(0, 4)
    fig.tight_layout()
    fig.savefig(F / "fig5_two_link.png")
    plt.close(fig)
    print("two-link velocity range", a[:, 5].min().round(3), a[:, 5].max().round(3), a[:, 6].min().round(3), a[:, 6].max().round(3))


diagram()
two_link()
responses()
ceilings()
saturation()
print("ok")

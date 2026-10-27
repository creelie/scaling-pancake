"""Figures of the paper "Weak-Type Endpoint Bounds for Bochner-Riesz Means in Every Dimension".

Writes fig_fibre.png, fig_cylinder.png and fig_sparse.png (200 dpi) next to this script.
Run from any directory:  python3 figures/make_figures.py
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10})

BLUE = "#2b5d8a"
RED = "#b23a3a"
GREEN = "#3f7f4f"
GREY = "#6b6b6b"


def save(fig, name):
    fig.savefig(os.path.join(HERE, name), dpi=200, bbox_inches="tight")
    plt.close(fig)


def fibre():
    """The fibre in the modulus coordinates (|W|, |Q|), with w0 = 1 and gamma = 1/10."""
    w0, gam = 1.0, 0.1
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    q = np.linspace(0, 0.1, 400)
    w = w0 * np.sqrt(np.clip(1 - 100 * q**2, 0, None))
    ax.fill_betweenx(q, 0, w, color=BLUE, alpha=0.08)
    ax.plot(w, q, color=BLUE, lw=1.2)
    heavy = q <= 0.05
    ax.plot(w[heavy], q[heavy], color=RED, lw=3.2, solid_capstyle="round")
    ax.axhline(0.05, color=GREY, lw=0.6, ls=":")
    ax.axvline(np.sqrt(3) / 2 * w0, color=GREY, lw=0.6, ls=":")
    ax.add_patch(plt.Rectangle((0, 0), 2 * gam, 0.01, fc=GREEN, alpha=0.18, ec=GREEN, lw=0.8))
    ax.add_patch(plt.Rectangle((0, 0), gam, 0.0025, fc=GREEN, alpha=0.45, ec=GREEN, lw=0.8))
    ax.annotate(r"wall $r=0$ carrying the flux, $|Q|\leq 1/20$", xy=(0.975, 0.02),
                xytext=(0.45, 0.012), color=RED,
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))
    ax.annotate(r"centre: $|W|\leq\gamma$, $|Q|\leq 1/400$", xy=(0.06, 0.0018),
                xytext=(0.24, 0.026), color=GREEN,
                arrowprops=dict(arrowstyle="->", color=GREEN, lw=0.8))
    ax.annotate(r"patch: $|W|\leq 2\gamma$, $|Q|\leq 1/100$, $g(Q)=1$", xy=(0.17, 0.009),
                xytext=(0.24, 0.042), color=GREEN,
                arrowprops=dict(arrowstyle="->", color=GREEN, lw=0.8))
    ax.text(np.sqrt(3) / 2 * w0 + 0.01, 0.091, r"$|W|=\frac{\sqrt{3}}{2}w_0$", color=GREY, fontsize=9)
    ax.text(0.01, 0.052, r"$|Q|=\frac{1}{20}$", color=GREY, fontsize=9)
    ax.set_xlim(0, 1.08)
    ax.set_ylim(0, 0.105)
    ax.set_xlabel(r"$|W|/w_0$")
    ax.set_ylabel(r"$|Q|$")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    save(fig, "fig_fibre.png")


def shell(ax, r, t0, t1, color, alpha, n=60):
    th = np.linspace(0, 2 * np.pi, n)
    t = np.linspace(t0, t1, 2)
    TH, TT = np.meshgrid(th, t)
    ax.plot_surface(r * np.cos(TH), r * np.sin(TH), TT, color=color, alpha=alpha, linewidth=0, shade=True)


def cylinder():
    """The two comparisons in the region b <= |tau| <= 3 (drawn to 3/2), with time T vertical."""
    b = 0.12
    fig = plt.figure(figsize=(7.2, 3.6))
    titles = [r"inner wall $s=b$ (red) to band $\frac{3}{4}\leq s\leq\frac{5}{4}$ (blue)",
              r"outer wall $s=1$ (blue) back to $s=2b$ (green)"]
    for k, title in enumerate(titles):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        ax.view_init(elev=18, azim=-58)
        shell(ax, 1.5, -1, 2, GREY, 0.06)
        if k == 0:
            shell(ax, b, 1, 2, RED, 0.85)
            shell(ax, 0.75, 1, 2, BLUE, 0.18)
            shell(ax, 1.25, 1, 2, BLUE, 0.18)
        else:
            shell(ax, 1.0, -1, 1, BLUE, 0.35)
            shell(ax, 2 * b, -1, 1, GREEN, 0.30)
            ax.scatter([2 * b], [0], [0], color="k", s=14)
            ax.text(2 * b + 0.1, 0.1, -0.35, r"$(0,2b)$", color="k")
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)
        ax.set_zlim(-1, 2.4)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([-1, 0, 1, 2])
        ax.set_xlabel(r"$\tau_1$", labelpad=-12)
        ax.set_ylabel(r"$\tau_2$", labelpad=-12)
        ax.set_zlabel(r"$T$", labelpad=2)
        ax.set_title(title, fontsize=9)
        ax.set_box_aspect((1, 1, 1.15))
    fig.subplots_adjust(wspace=0.02)
    save(fig, "fig_cylinder.png")


def sparse(n=3):
    """The critical sides of the sparse trapezoids for n = 3: (1/p, 1/q, lambda)."""
    fig = plt.figure(figsize=(6.4, 4.6))
    ax = fig.add_subplot(projection="3d")
    ax.view_init(elev=22, azim=-125)
    lam = np.linspace(0.0, (n - 1) / 2, 60)
    x = (2 * lam + n + 1) / (2 * n)  # 1/p with lambda = lambda(p)
    y1 = (n - 2 * lam - 1) / (2 * n)  # P1: 1/q = 1 - 1/p
    y2 = (n - 1) / (2 * n) + lam * (n + 1) / (n * (n - 1))  # P2
    # surface swept by P1P2 (Theorem B), and its dual P4P3
    for (xa, ya, xb, yb, col) in [(x, y1, x, y2, BLUE), (y1, x, y2, x, GREEN)]:
        polys = []
        for i in range(len(lam) - 1):
            polys.append([(xa[i], ya[i], lam[i]), (xb[i], yb[i], lam[i]),
                          (xb[i + 1], yb[i + 1], lam[i + 1]), (xa[i + 1], ya[i + 1], lam[i + 1])])
        ax.add_collection3d(Poly3DCollection(polys, facecolor=col, alpha=0.28, edgecolor="none"))
    ax.plot(x, y2, lam, color=RED, lw=2.4, label=r"corner $P_2$ (open)")
    ax.plot(y2, x, lam, color=RED, lw=2.4, ls="--")
    for l0 in [0.25, 0.5, 0.75]:
        xx = (2 * l0 + n + 1) / (2 * n)
        a = (n - 2 * l0 - 1) / (2 * n)
        b = (n - 1) / (2 * n) + l0 * (n + 1) / (n * (n - 1))
        P = np.array([[xx, a], [xx, b], [b, xx], [a, xx], [xx, a]])
        ax.plot(P[:, 0], P[:, 1], l0, color=GREY, lw=0.8)
    ax.plot(x, y1, lam, color=BLUE, lw=1.2, label=r"side $P_1P_2$: Theorem B")
    ax.plot(y1, x, lam, color=GREEN, lw=1.2, label=r"dual side $P_4P_3$")
    ax.set_xlabel(r"$1/p$")
    ax.set_ylabel(r"$1/q$")
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel(r"$\lambda$")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 1)
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    save(fig, "fig_sparse.png")


if __name__ == "__main__":
    fibre()
    cylinder()
    sparse()

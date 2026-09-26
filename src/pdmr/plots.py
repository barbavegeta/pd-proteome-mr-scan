"""Static figures for the scan report (matplotlib, PNG)."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .ranking import TIERS  # noqa: E402

INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#b9b8b2", "#e7e6e2", "#fcfcfb"
# Only three tiers get a hue (scatter charts cap at three validated categorical slots);
# the rest are neutral and identified by legend + marker shape.
TIER_STYLE = {
    TIERS[0]: ("#2a78d6", "o", 9, 3),
    TIERS[1]: ("#eb6834", "D", 7, 3),
    TIERS[2]: ("#1baf7a", "s", 7, 3),
    TIERS[3]: ("#8f8e88", "^", 6, 2),
    TIERS[4]: (MUTED, "o", 4, 1),
}


def _style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def volcano(res: pd.DataFrame, path, label_tiers=3, highlight=()):
    d = res[res.mr_p.notna()].copy()
    d["nlp"] = -np.log10(d.mr_p.clip(lower=1e-300))
    d["lb"] = np.log(d.mr_OR)
    fig, ax = plt.subplots(figsize=(8, 5.2), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    _style(ax)
    for tier in reversed(list(TIER_STYLE)):
        g = d[d.tier == tier]
        if g.empty:
            continue
        col, mk, sz, z = TIER_STYLE[tier]
        ax.scatter(g.lb, g.nlp, s=sz ** 2, c=col, marker=mk, edgecolors=SURFACE, linewidths=1.2,
                   zorder=z, label=f"{tier} (n={len(g)})")
    ax.axvline(0, color=MUTED, lw=1)
    lab = d[d.tier.isin(TIERS[:label_tiers]) | d.gene.isin(highlight)]
    for r in lab.itertuples():
        ax.annotate(r.gene, (r.lb, r.nlp), xytext=(5, 4), textcoords="offset points", fontsize=8.5, color=INK)
    ax.set_xlabel("MR effect on Parkinson’s disease, log OR per SD protein", color=INK2, fontsize=10)
    ax.set_ylabel("−log$_{10}$ MR p-value", color=INK2, fontsize=10)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="upper center",
              bbox_to_anchor=(0.5, -0.14), ncol=3, handletextpad=0.3, columnspacing=1.2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def forest(res: pd.DataFrame, path, tiers=TIERS[:3], max_rows=40):
    d = res[res.tier.isin(tiers)].head(max_rows).iloc[::-1]
    if d.empty:
        return False
    h = 1.2 + 0.32 * len(d)
    fig, ax = plt.subplots(figsize=(8, h), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    _style(ax)
    ax.grid(axis="y", visible=False)
    y = np.arange(len(d))
    cols = [TIER_STYLE[t][0] for t in d.tier]
    ax.hlines(y, d.mr_OR_lci, d.mr_OR_uci, colors=cols, lw=2)
    ax.scatter(d.mr_OR, y, c=cols, s=40, zorder=3, edgecolors=SURFACE, linewidths=1.2)
    ax.axvline(1, color=MUTED, lw=1)
    ax.set_xscale("log")
    from matplotlib.ticker import FuncFormatter, LogLocator
    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1, 1.25, 1.5, 2, 3, 5, 7)))
    ax.xaxis.set_minor_locator(LogLocator(base=10, subs="auto"))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))
    ax.set_yticks(y, [f"{g}  (PP.H4 {h4:.2f})" for g, h4 in zip(d.gene, d["coloc_PP.H4"])], fontsize=8.5, color=INK)
    ax.set_xlabel("Odds ratio for Parkinson’s disease per SD protein (95% CI)", color=INK2, fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return True


def regional(region: pd.DataFrame, instruments: pd.DataFrame, title, path, subtitle=""):
    fig, axes = plt.subplots(2, 1, figsize=(8, 4.9), dpi=150, sharex=True)
    fig.patch.set_facecolor(SURFACE)
    x = region.pos37 / 1e6
    for ax, col, lab in ((axes[0], "p_exp", "pQTL (protein)"), (axes[1], "p_out", "Parkinson’s disease")):
        _style(ax)
        y = -np.log10(region[col].clip(lower=1e-300))
        ax.scatter(x, y, s=10, c=MUTED, edgecolors="none")
        if len(instruments):
            m = region.pos37.isin(instruments.pos37)
            ax.scatter(x[m], y[m], s=64, c="#2a78d6", edgecolors=SURFACE, linewidths=1.5, zorder=3,
                       label="MR instrument")
        ax.set_ylabel(f"−log$_{{10}}$ p\n{lab}", color=INK2, fontsize=9)
    axes[0].set_title(title + ("\n" + subtitle if subtitle else ""), loc="left", fontsize=10, color=INK)
    if len(instruments):
        axes[0].legend(frameon=False, fontsize=8, labelcolor=INK2)
    axes[1].set_xlabel("Position (GRCh37, Mb)", color=INK2, fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)

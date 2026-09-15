"""Shared plotting theme and data loaders for the gallrna notebooks.

Kept in one place so every notebook draws from the same palette: a host keeps
its colour in every figure, whatever subset is on screen.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- paths ------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
if not (ROOT / "samples.tsv").exists():                # fall back to cwd
    p = Path.cwd().resolve()
    ROOT = next(c for c in (p, *p.parents) if (c / "samples.tsv").exists())
FIGDIR = ROOT / "figures"

# --- palette ----------------------------------------------------------------
# Fixed slot order; hues are assigned in order and never cycled.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
          "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
BLUES = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
         "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
SEQ = mpl.colors.LinearSegmentedColormap.from_list("seq_blue", BLUES)

# Roles that recur across notebooks.
AGRO, PLANT = SERIES[1], SERIES[2]
STRAIN_COLOR = {"1416": SERIES[0], "29": SERIES[1]}
HOST_ORDER = ["Euonymus japonicus", "Citrus sinensis", "Poncirus trifoliata",
              "Brassica juncea", "Carica papaya", "Solanum lycopersicum"]
HOST_COLOR = {h: SERIES[i] for i, h in enumerate(HOST_ORDER)}
REPL = ["circ", "lin", "pAt1", "pAt2", "pTi"]
REPL_LABEL = {"circ": "circular chromosome", "lin": "linear chromosome",
              "pAt1": "pAt1", "pAt2": "pAt2", "pTi": "Ti plasmid"}
REPL_COLOR = {r: SERIES[i] for i, r in enumerate(REPL)}


def use_theme():
    mpl.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "figure.dpi": 110, "savefig.dpi": 150, "savefig.bbox": "tight",
        "font.size": 9, "text.color": INK,
        "axes.edgecolor": GRID, "axes.labelcolor": INK2, "axes.titlesize": 11,
        "axes.titleweight": "semibold", "axes.titlecolor": INK, "axes.titlelocation": "left",
        "axes.titlepad": 10, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.frameon": False, "legend.fontsize": 8,
        "lines.linewidth": 2, "lines.markersize": 8,
    })
    FIGDIR.mkdir(exist_ok=True)


# --- chart furniture --------------------------------------------------------
def clean(ax, axis="x"):
    """Recessive frame: grid on the measure axis only, no top/right spines."""
    hide = {"x": ("top", "right", "left"), "y": ("top", "right", "bottom"),
            "both": ("top", "right")}[axis]
    for side in hide:
        ax.spines[side].set_visible(False)
    ax.grid(axis="both" if axis == "both" else axis, alpha=0.9)
    ax.set_axisbelow(True)
    return ax


def note(ax, text, y=-0.34):
    ax.text(0, y, text, transform=ax.transAxes, fontsize=7.5, color=INK2,
            va="top", ha="left", wrap=True)


def barh(ax, labels, values, colors, fmt="{:,.0f}", pad=0.02, log=False):
    values = np.asarray(values, dtype=float)
    y = np.arange(len(labels))
    ax.barh(y, values, height=0.62, color=colors)      # height<1 leaves a gap between bars
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    span = values.max() if len(values) else 1
    for yi, v in zip(y, values):
        # offset has to be multiplicative on a log axis, additive on a linear one
        ax.text(v * 1.12 if log else v + span * pad, yi, fmt.format(v),
                va="center", fontsize=7.5, color=INK2)
    if not log:
        ax.set_xlim(0, span * 1.18)
    return clean(ax, "x")


def host_legend(ax, hosts=None, y=-0.17, ncols=3):
    """Legend below the plot: with 14 sorted bars there is no free space inside."""
    hosts = hosts if hosts is not None else HOST_ORDER
    handles = [mpl.patches.Patch(facecolor=HOST_COLOR[h], label=h) for h in hosts]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, y),
              ncols=ncols, borderaxespad=0, columnspacing=1.4)


def save(fig, name):
    fig.savefig(FIGDIR / f"{name}.png")
    return fig


# --- schematic furniture ----------------------------------------------------
def box(ax, x, y, w, h, text, face=SURFACE, edge=None, fontsize=8.5,
        weight="normal", text_color=None, radius=0.02):
    """Rounded box with centred text, in axes coordinates."""
    edge = edge or GRID
    ax.add_patch(mpl.patches.FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0,rounding_size={radius}",
        facecolor=face, edgecolor=edge, linewidth=1.4, transform=ax.transAxes,
        clip_on=False, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, transform=ax.transAxes, ha="center",
            va="center", fontsize=fontsize, color=text_color or INK,
            weight=weight, zorder=3, linespacing=1.5)


def arrow(ax, xy_from, xy_to, color=None, style="-|>", lw=1.6):
    ax.add_patch(mpl.patches.FancyArrowPatch(
        xy_from, xy_to, arrowstyle=style, mutation_scale=13,
        color=color or INK2, linewidth=lw, transform=ax.transAxes,
        clip_on=False, zorder=1, shrinkA=2, shrinkB=2))


def blank(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis("off")
    return ax


# --- data loaders -----------------------------------------------------------
def load_meta():
    m = pd.read_csv(ROOT / "samples.tsv", sep="\t")
    m["gall_age_days"] = pd.to_numeric(m["gall_age_days"], errors="coerce")
    m["gall_mass_g"] = pd.to_numeric(m["gall_mass_g"], errors="coerce")
    m["strain"] = m["strain"].astype(str)
    m = m.set_index("sample")
    order = (m.assign(_h=m["host_species"].map({h: i for i, h in enumerate(HOST_ORDER)}))
              .sort_values(["_h", "strain", "gall_age_days"]).index.tolist())
    return m.loc[order]


def load_idxstats(sample):
    return pd.read_csv(ROOT / "02_align" / f"{sample}.idxstats", sep="\t",
                       names=["ref", "length", "mapped", "unmapped"]).query("ref != '*'")


def load_gff_genes(path):
    rec = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or "\tgene\t" not in line:
            continue
        f = line.split("\t")
        attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        rec.append(dict(gene=attr.get("ID"), contig=f[0], start=int(f[3]), end=int(f[4]),
                        strand=f[6], name=attr.get("Name", "")))
    return pd.DataFrame(rec)

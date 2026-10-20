"""Figures for the final segmentation validation.

  output/segmentation_validation/Figure_S_seg_overlays.png
      Representative overlays: one slice of a low-density (3 h) and a
      high-density (7 h) + CAA image per carbon, StarDist outlines on the raw
      image. Once points/ exist, outlines are coloured by outcome (matched,
      merged = 2+ clicks, unmatched = no click) and clicks are shown.
  output/segmentation_validation/Figure_S_seg_metrics.png   (needs results/)
      a manual vs automated counts, b count error and c recall by timepoint
      (grey band = 3-7 h growth-rate window).

Run from the repository root (after score_final.py):
  python validation/final/figures_final.py
"""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from skimage.segmentation import find_boundaries

sys.path.insert(0, str(Path("validation")))
import score_validation as sv  # noqa: E402

FINAL = Path("validation/final")
RES = FINAL / "results"
OUT = Path("output/segmentation_validation")
XY_UM = 0.0916
CARBONS = ["glucose", "glycerol", "succinate"]
COL = {"glucose": "#E69F00", "glycerol": "#56B4E9", "succinate": "#009E73"}  # Okabe-Ito, as in the R figures
STATUS = {"matched": (0.25, 0.85, 1.0), "merged": (1.0, 0.6, 0.0), "unmatched": (1.0, 0.2, 0.2)}

plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 9, "axes.titleweight": "bold", "axes.titlelocation": "left"})


def read_csv(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def overlay(entry):
    vid = entry["id"]
    raw = tifffile.imread(FINAL / "images" / f"{vid}.tif").astype(float)
    lab = tifffile.imread(FINAL / "_labels" / f"{vid}.tif")
    eligible = set(map(int, (FINAL / "_eligible" / f"{vid}.txt").read_text().split()))
    pts_path = FINAL / "points" / f"{vid}.csv"
    points = sv.read_points(pts_path) if pts_path.exists() else []
    hits = {}
    for p in points:
        o = sv.match_point(lab, *p)
        hits[o] = hits.get(o, 0) + 1
    z = int(np.argmax((lab > 0).sum(axis=(1, 2))))
    sl, l2 = raw[z], lab[z]
    # zoom whole stacks to the cells
    ys, xs = np.nonzero(lab.max(0))
    pad = 25
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad, sl.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad, sl.shape[1])
    side = max(y1 - y0, x1 - x0)
    y1, x1 = min(y0 + side, sl.shape[0]), min(x0 + side, sl.shape[1])
    sl, l2 = sl[y0:y1, x0:x1], l2[y0:y1, x0:x1]
    lo, hi = np.percentile(sl, [1, 99.7])
    rgb = np.dstack([np.clip((sl - lo) / (hi - lo + 1e-9), 0, 1)] * 3)
    edges = find_boundaries(l2, mode="inner")
    for obj in np.unique(l2[edges]):
        if obj == 0:
            continue
        if not points:
            colour = STATUS["matched"]
        elif obj not in eligible:
            colour = (0.55, 0.55, 0.55)
        else:
            n = hits.get(obj, 0)
            colour = STATUS["unmatched" if n == 0 else "merged" if n > 1 else "matched"]
        rgb[edges & (l2 == obj)] = colour
    clicks = [(py - y0, px - x0) for pz, py, px in points
              if abs(pz - z) <= 1 and y0 <= py < y1 and x0 <= px < x1]
    return rgb, clicks, bool(points)


def overlay_figure():
    manifest = read_csv(FINAL / "_manifest.csv")

    def pick(carbon, t):
        return next(e for e in manifest if e["carbon"] == carbon and e["medium"] == "caa"
                    and int(e["timepoint_h"]) == t)

    fig, axes = plt.subplots(2, 3, figsize=(7.1, 5.0))
    scored = False
    for col, carbon in enumerate(CARBONS):
        for row, t in enumerate((3, 7)):
            e = pick(carbon, t)
            rgb, clicks, scored = overlay(e)
            ax = axes[row, col]
            ax.imshow(rgb, interpolation="nearest")
            for y, x in clicks:
                ax.plot(x, y, "+", color="yellow", ms=5, mew=1)
            bar = 2 / XY_UM
            h, w = rgb.shape[:2]
            ax.plot([w - 6 - bar, w - 6], [h - 6, h - 6], color="white", lw=2)
            ax.set_xticks([]), ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
            where = "whole stack" if e["kind"] == "whole" else "crop of densest layer"
            ax.set_title(f"{carbon.capitalize()}, {t} h ({where})",
                         fontsize=8, fontweight="normal")
    note = ("Outlines: cyan matched, orange merged (2+ clicks), red no click, grey "
            "mostly outside crop; + manual clicks (±1 slice). Scale bar 2 µm."
            if scored else "StarDist outlines on the densest slice. Scale bar 2 µm.")
    fig.text(0.01, 0.005, note, fontsize=7)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "Figure_S_seg_overlays.png", dpi=300)
    plt.close(fig)


def metrics_figure():
    a = read_csv(RES / "partA_per_image.csv")
    for r in a:
        for k in ("timepoint_h", "manual_cells", "automated_objects", "precision", "recall", "count_error_pct"):
            r[k] = float(r[k])
    fig, ax = plt.subplots(1, 3, figsize=(7.1, 2.5))
    jitter = {"glucose": -0.15, "glycerol": 0, "succinate": 0.15}

    def mk(r):
        return dict(marker="o", ms=4.5, mec=COL[r["carbon"]], mew=1,
                    mfc=COL[r["carbon"]] if r["medium"] == "caa" else "white", ls="none")

    for r in a:
        ax[0].plot(r["manual_cells"], r["automated_objects"], **mk(r))
    lim = [1, max(max(r["manual_cells"], r["automated_objects"]) for r in a) * 1.3]
    ax[0].plot(lim, lim, color="0.5", lw=0.8, zorder=0)
    ax[0].set(xscale="log", yscale="log", xlim=lim, ylim=lim,
              xlabel="Manual count", ylabel="Automated count", title="a")

    for i, (key, label, title) in enumerate((("count_error_pct", "Count error (%)", "b"),
                                             ("recall", "Recall", "c")), start=1):
        ax[i].axvspan(2.7, 7.3, color="0.93", zorder=0)
        for c in CARBONS:
            rows = sorted((r for r in a if r["carbon"] == c and r["medium"] == "caa"),
                          key=lambda r: r["timepoint_h"])
            ax[i].plot([r["timepoint_h"] + jitter[c] for r in rows], [r[key] for r in rows],
                       color=COL[c], lw=1, alpha=0.7)
            for r in a:
                if r["carbon"] == c:
                    ax[i].plot(r["timepoint_h"] + jitter[c], r[key], **mk(r))
        ax[i].set(xlabel="Time (h)", ylabel=label, title=title, xticks=range(2, 9))
        if key == "count_error_pct":
            ax[i].axhline(0, color="0.5", lw=0.8, zorder=0)
        else:
            ax[i].set_ylim(0, 1.05)

    for c in CARBONS:
        ax[0].plot([], [], "o", color=COL[c], label=c.capitalize())
    ax[0].plot([], [], "o", mfc="white", mec="0.3", label="−CAA")
    ax[0].legend(fontsize=6.5, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "Figure_S_seg_metrics.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    overlay_figure()
    if (RES / "partA_per_image.csv").exists():
        metrics_figure()
    else:
        print("No results yet: made the overlay figure only.")

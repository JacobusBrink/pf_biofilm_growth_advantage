"""Set up the segmentation-validation subset (Reviewer comment 5).

Reads StarDist label stacks (ch2) and raw stacks (ch1) from the processed
image folders on E: (READ-ONLY - nothing there is modified) and writes COPIES
of a stratified subset into this folder:

  validation/images/Vxx.tif      raw image to annotate in Fiji (blinded name)
  validation/_labels/Vxx.tif     matching StarDist labels (do not open while
                                 annotating)
  validation/_manifest.csv       which sample each Vxx is, crop box, counts

Subset:
  - low density, whole stacks: 2 h and 4 h for glucose, glycerol and
    succinate, with and without CAA (12 stacks)
  - high density, crops holding ~30 cells: a 9-slice band around the
    densest layer and a square xy window, 6 h and 8 h for glucose, glycerol
    and succinate with CAA (6 crops)

For every sample the label file used is the one whose number of objects
equals the BiofilmQ Ncells in input/biofilm/, i.e. the segmentation that
produced the published numbers.

Run from the repository root:  python validation/setup_validation.py
"""

import csv
import random
import re
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage

ROOTS = [
    Path("E:/early-stage_growth_carbons_processed"),
    Path("E:/early-stage_growth_carbons_controls_processed"),
]
CSV_DIR = Path("input/biofilm")
OUT = Path("validation")

TARGET_CELLS_PER_CROP = 30
Z_SLICES_PER_CROP = 9  # ~2.8 um: one to two cell layers
SEED = 365  # fixed so the blinded order is reproducible

NAME_RE = re.compile(r"^([a-z]+?)(\d+)hbr(\d+)tr(\d+)_pos1_ch2_time1\.tif$")


def biofilmq_ncells(cond, t, br, tr):
    """Ncells from the BiofilmQ summary CSV used in the paper (None if absent)."""
    path = CSV_DIR / f"{cond}_{t}h_br{br}_tr{tr}.csv"
    if not path.exists():
        return None
    with path.open(newline="") as fh:
        # BiofilmQ writes an empty units row before the data row
        rows = [r for r in csv.DictReader(fh) if r["Ncells"]]
    return int(float(rows[0]["Ncells"]))


def index_label_files():
    """All ch2 label stacks on E:, grouped by sample name."""
    index = {}
    for root in ROOTS:
        for path in root.rglob("*_ch2_time1.tif"):
            m = NAME_RE.match(path.name)
            if not m:
                continue
            cond, t, br, tr = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
            index.setdefault((cond, t, br, tr), []).append(path)
    return index


def matching_label_file(candidates, ncells):
    """The candidate whose object count equals the published Ncells."""
    for path in candidates:
        labels = tifffile.imread(path)
        if len(np.unique(labels)) - 1 == ncells:
            return path, labels
    return None, None


def pick_sample(index, cond, t, br, max_cells=None):
    """First technical replicate with a label file matching BiofilmQ."""
    for tr in range(1, 16):
        ncells = biofilmq_ncells(cond, t, br, tr)
        if ncells is None or (max_cells and ncells > max_cells):
            continue
        path, labels = matching_label_file(index.get((cond, t, br, tr), []), ncells)
        if path is not None:
            return tr, ncells, path, labels
    raise RuntimeError(f"no matching label file for {cond} {t}h br{br}")


def dense_crop_box(labels, target, z_slices=Z_SLICES_PER_CROP):
    """Thin z band around the densest slice, with a square xy window around
    the colony centre grown until it holds ~target objects.

    An object counts toward the window when at least half of its volume lies
    inside the box (the same rule the scoring script uses).
    """
    zdim, ydim, xdim = labels.shape
    z_peak = int(np.argmax((labels > 0).sum(axis=(1, 2))))
    z0 = max(z_peak - z_slices // 2, 0)
    z1 = min(z0 + z_slices, zdim)
    total = np.bincount(labels.ravel())

    ids = np.unique(labels[z0:z1])
    ids = ids[ids > 0]
    centroids = np.array(ndimage.center_of_mass(labels > 0, labels, ids))
    cy, cx = np.median(centroids[:, 1]), np.median(centroids[:, 2])
    for half in range(20, ydim // 2, 2):
        y0, y1 = max(int(cy - half), 0), min(int(cy + half), ydim)
        x0, x1 = max(int(cx - half), 0), min(int(cx + half), xdim)
        inside = np.bincount(labels[z0:z1, y0:y1, x0:x1].ravel(), minlength=len(total))
        with np.errstate(divide="ignore", invalid="ignore"):
            eligible = np.sum(inside[1:] / total[1:] >= 0.5)
        if eligible >= target:
            return z0, z1, y0, y1, x0, x1
    return z0, z1, 0, ydim, 0, xdim


def main():
    (OUT / "images").mkdir(parents=True, exist_ok=True)
    (OUT / "_labels").mkdir(parents=True, exist_ok=True)
    index = index_label_files()

    plan = []
    reps = [1, 2, 3]
    i = 0
    for cond in ["glucose", "glycerol", "succinate",
                 "glucosecontrol", "glycerolcontrol", "succinatecontrol"]:
        for t in [2, 4]:
            plan.append(("whole", cond, t, reps[i % 3]))
            i += 1
    for cond in ["glucose", "glycerol", "succinate"]:
        for t in [6, 8]:
            plan.append(("crop", cond, t, reps[i % 3]))
            i += 1

    random.Random(SEED).shuffle(plan)

    manifest = []
    for n, (kind, cond, t, br) in enumerate(plan, start=1):
        vid = f"V{n:02d}"
        tr, ncells, label_path, labels = pick_sample(
            index, cond, t, br, max_cells=25 if kind == "whole" else None
        )
        raw_path = label_path.with_name(label_path.name.replace("_ch2_", "_ch1_"))
        raw = tifffile.imread(raw_path)

        if kind == "crop":
            z0, z1, y0, y1, x0, x1 = dense_crop_box(labels, TARGET_CELLS_PER_CROP)
        else:
            z0, z1, y0, y1, x0, x1 = 0, labels.shape[0], 0, labels.shape[1], 0, labels.shape[2]

        raw_crop = raw[z0:z1, y0:y1, x0:x1]
        lab_crop = labels[z0:z1, y0:y1, x0:x1]
        tifffile.imwrite(OUT / "images" / f"{vid}.tif", raw_crop, imagej=True,
                         metadata={"axes": "ZYX"})
        tifffile.imwrite(OUT / "_labels" / f"{vid}.tif", lab_crop, imagej=True,
                         metadata={"axes": "ZYX"})

        manifest.append({
            "id": vid, "kind": kind, "condition": cond, "timepoint_h": t,
            "biological_replicate": br, "technical_replicate": tr,
            "biofilmq_ncells": ncells,
            "objects_in_crop": len(np.unique(lab_crop)) - 1,
            "z0": z0, "z1": z1, "y0": y0, "y1": y1, "x0": x0, "x1": x1,
            "z_slices": raw_crop.shape[0],
            "raw_source": str(raw_path), "label_source": str(label_path),
        })
        print(vid, kind, cond, f"{t}h", f"br{br}tr{tr}", "objects:",
              manifest[-1]["objects_in_crop"])

    with (OUT / "_manifest.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)


if __name__ == "__main__":
    main()

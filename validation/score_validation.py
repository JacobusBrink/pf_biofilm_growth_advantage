"""Score StarDist segmentation against manual point annotations.

For every validation/images/Vxx.tif that has a validation/points/Vxx.csv
(Fiji Results table: X, Y, Slice, one row per clicked cell), each click is
matched to the StarDist object it falls in (or the nearest object within
MATCH_RADIUS_XY pixels and one z-slice), and:

  true positive   object hit by exactly one click
  merged          object hit by 2+ clicks: 1 true positive, the rest
                  counted as missed cells (under-segmentation)
  missed          click matching no object (false negative)
  unmatched       object with no click: false positive or split cell
                  (over-segmentation)

In crops, objects with less than half of their volume inside the crop are
ignored, together with any clicks on them, so cells cut by the crop edge do
not count for or against the segmentation.

Outputs:
  validation/results/per_image.csv
  validation/results/summary.csv   (by density: low = whole stacks,
                                    high = 6-8 h crops, and overall)

Run from the repository root:  python validation/score_validation.py
Self-test (clicks placed on object centroids, must score 100%):
  python validation/score_validation.py --selftest
"""

import csv
import sys
import tempfile
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage

OUT = Path("validation")
MATCH_RADIUS_XY = 3  # pixels (~0.27 um)
MATCH_RADIUS_Z = 1  # slices (~0.31 um)
MIN_FRACTION_IN_CROP = 0.5


def read_manifest():
    with (OUT / "_manifest.csv").open(newline="") as fh:
        return {row["id"]: row for row in csv.DictReader(fh)}


def read_points(path):
    """Clicks as (z, y, x) voxel indices from a Fiji Results CSV."""
    with path.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    if rows and "Slice" not in {k.strip() for k in rows[0] if k is not None}:
        raise ValueError(
            f"{path.name} has no Slice column: in Fiji tick 'Stack position' "
            "under Analyze > Set Measurements, then measure the points again"
        )
    points = []
    for row in rows:
        row = {k.strip(): v for k, v in row.items() if k is not None}
        z = int(float(row.get("Slice", 1))) - 1
        points.append((z, int(float(row["Y"])), int(float(row["X"]))))
    return points


def eligible_objects(labels, entry):
    """Object ids counted in this image (>= half their volume inside)."""
    ids, counts_in_crop = np.unique(labels, return_counts=True)
    counts_in_crop = dict(zip(ids[ids > 0], counts_in_crop[ids > 0]))
    if entry["kind"] != "crop":
        return set(counts_in_crop)
    full = tifffile.imread(entry["label_source"])
    full_ids, full_counts = np.unique(full, return_counts=True)
    full_counts = dict(zip(full_ids, full_counts))
    return {
        i for i, n in counts_in_crop.items()
        if n / full_counts[i] >= MIN_FRACTION_IN_CROP
    }


def match_point(labels, z, y, x):
    """Object id at the click, else the nearest one within the radius."""
    zmax, ymax, xmax = labels.shape
    z, y, x = min(max(z, 0), zmax - 1), min(max(y, 0), ymax - 1), min(max(x, 0), xmax - 1)
    if labels[z, y, x]:
        return int(labels[z, y, x])
    z0, z1 = max(z - MATCH_RADIUS_Z, 0), min(z + MATCH_RADIUS_Z + 1, zmax)
    y0, y1 = max(y - MATCH_RADIUS_XY, 0), min(y + MATCH_RADIUS_XY + 1, ymax)
    x0, x1 = max(x - MATCH_RADIUS_XY, 0), min(x + MATCH_RADIUS_XY + 1, xmax)
    window = labels[z0:z1, y0:y1, x0:x1]
    zz, yy, xx = np.nonzero(window)
    if len(zz) == 0:
        return 0
    d = (zz + z0 - z) ** 2 * 9 + (yy + y0 - y) ** 2 + (xx + x0 - x) ** 2
    k = int(np.argmin(d))
    return int(window[zz[k], yy[k], xx[k]])


def score_image(labels, points, eligible):
    hits = {}
    missed = 0
    edge_clicks = 0
    for z, y, x in points:
        obj = match_point(labels, z, y, x)
        if obj == 0:
            missed += 1
        elif obj not in eligible:
            edge_clicks += 1
        else:
            hits[obj] = hits.get(obj, 0) + 1
    tp = len(hits)
    merged_objects = sum(1 for n in hits.values() if n > 1)
    merged_extra = sum(n - 1 for n in hits.values() if n > 1)
    unmatched = len(eligible - set(hits))
    fn = missed + merged_extra
    return {
        "manual_cells": len(points) - edge_clicks,
        "automated_objects": len(eligible),
        "true_positive": tp,
        "merged_objects": merged_objects,
        "missed_cells": missed,
        "unmatched_objects": unmatched,
        "false_negative": fn,
        "false_positive": unmatched,
    }


def metrics(tp, fp, fn):
    precision = tp / (tp + fp) if tp + fp else float("nan")
    recall = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else float("nan")
    return round(precision, 3), round(recall, 3), round(f1, 3)


def run(points_dir, results_dir):
    manifest = read_manifest()
    per_image = []
    for vid, entry in sorted(manifest.items()):
        pts_path = points_dir / f"{vid}.csv"
        if not pts_path.exists():
            continue
        labels = tifffile.imread(OUT / "_labels" / f"{vid}.tif")
        eligible = eligible_objects(labels, entry)
        row = {
            "id": vid, "kind": entry["kind"], "condition": entry["condition"],
            "timepoint_h": entry["timepoint_h"],
            **score_image(labels, read_points(pts_path), eligible),
        }
        row["precision"], row["recall"], row["f1"] = metrics(
            row["true_positive"], row["false_positive"], row["false_negative"]
        )
        row["count_error_pct"] = round(
            100 * (row["automated_objects"] - row["manual_cells"]) / row["manual_cells"], 1
        ) if row["manual_cells"] else float("nan")
        per_image.append(row)

    if not per_image:
        print("No point files found in", points_dir)
        return []

    results_dir.mkdir(parents=True, exist_ok=True)
    with (results_dir / "per_image.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(per_image[0]))
        writer.writeheader()
        writer.writerows(per_image)

    summary = []
    groups = {
        "low density (2-4 h, whole stacks)": [r for r in per_image if r["kind"] == "whole"],
        "high density (6-8 h, crops)": [r for r in per_image if r["kind"] == "crop"],
        "all": per_image,
    }
    for name, rows in groups.items():
        if not rows:
            continue
        tp = sum(r["true_positive"] for r in rows)
        fp = sum(r["false_positive"] for r in rows)
        fn = sum(r["false_negative"] for r in rows)
        p, rc, f1 = metrics(tp, fp, fn)
        summary.append({
            "group": name, "images": len(rows),
            "manual_cells": sum(r["manual_cells"] for r in rows),
            "automated_objects": sum(r["automated_objects"] for r in rows),
            "true_positive": tp, "false_positive": fp, "false_negative": fn,
            "merged_objects": sum(r["merged_objects"] for r in rows),
            "precision": p, "recall": rc, "f1": f1,
        })
    with (results_dir / "summary.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    for row in summary:
        print(row)
    return summary


def selftest():
    """Place one click on each object centroid; every metric must be 1."""
    manifest = read_manifest()
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for vid, entry in manifest.items():
            labels = tifffile.imread(OUT / "_labels" / f"{vid}.tif")
            eligible = sorted(eligible_objects(labels, entry))
            with (tmp / f"{vid}.csv").open("w", newline="") as fh:
                fh.write(" ,X,Y,Slice\n")
                for n, obj in enumerate(eligible, start=1):
                    zz, yy, xx = np.nonzero(labels == obj)
                    # the object voxel closest to its centroid (always inside it)
                    c = np.array([zz.mean(), yy.mean(), xx.mean()])
                    k = int(np.argmin((zz - c[0]) ** 2 + (yy - c[1]) ** 2 + (xx - c[2]) ** 2))
                    fh.write(f"{n},{xx[k] + 0.5},{yy[k] + 0.5},{zz[k] + 1}\n")
        summary = run(tmp, tmp / "results")
    assert all(r["precision"] == 1 and r["recall"] == 1 for r in summary), "self-test failed"
    print("Self-test passed.")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
    else:
        run(OUT / "points", OUT / "results")

"""Score the final segmentation validation (see setup_final.py / README.md).

Part A - clicks (points/Fxx.csv, Fiji Results: X, Y, Slice) vs 3D StarDist
labels, with the round-1/2 matching rules (validation/score_validation.py,
not modified): precision, recall, F1, merged objects, count error; plus the
mean nearest-neighbour distance from clicks vs from object centroids, and the
bias that the count errors would put on the specific growth rate (slope of
ln(automated / manual) against time, 3-7 h, per carbon, + CAA).

Part B - outlines (outlines/Bxx.zip, Fiji ROI Manager) vs the label slice:
objects matched one-to-one by IoU (Hungarian), F1 at IoU 0.5-0.9, mean IoU,
area ratio, centroid offset, orientation error and the 2D nematic order
parameter of manual vs automated cells.

Optional second annotator: points_rater2/Fxx.csv are scored the same way and
reported next to rater 1 (annotator variability).

Outputs in validation/final/results/.
Run from the repository root:   python validation/final/score_final.py
Self-test (synthetic perfect annotations must score 1):
                                python validation/final/score_final.py --selftest
"""

import csv
import math
import sys
import tempfile
from pathlib import Path

import numpy as np
import roifile
import tifffile
from scipy import ndimage
from scipy.optimize import linear_sum_assignment
from skimage.draw import polygon2mask
from skimage.measure import find_contours, regionprops

sys.path.insert(0, str(Path("validation")))
import score_validation as sv  # noqa: E402

FINAL = Path("validation/final")
RESULTS = FINAL / "results"
XY_UM, Z_UM = 0.0916, 0.31  # voxel size of the CZI stacks
IOU_THRESHOLDS = (0.5, 0.6, 0.7, 0.8, 0.9)
MIN_ELONGATION = 1.5  # major/minor axis ratio for orientation analysis
MIN_FP_AREA_FRACTION = 0.25  # of median manual cell area (out-of-focus tips)
MU_WINDOW = (3, 7)


def read_csv(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def pct(auto, manual):
    return round(100 * (auto - manual) / manual, 1) if manual else float("nan")


# --- Part A ----------------------------------------------------------------

def mean_nn_um(points_zyx):
    if len(points_zyx) < 2:
        return float("nan")
    p = np.asarray(points_zyx, float) * [Z_UM, XY_UM, XY_UM]
    d = np.sqrt(((p[:, None] - p[None]) ** 2).sum(-1))
    np.fill_diagonal(d, np.inf)
    return round(float(d.min(1).mean()), 3)


def score_part_a(points_dir):
    rows = []
    for e in read_csv(FINAL / "_manifest.csv"):
        pts_path = points_dir / f"{e['id']}.csv"
        if not pts_path.exists():
            continue
        labels = tifffile.imread(FINAL / "_labels" / f"{e['id']}.tif")
        text = (FINAL / "_eligible" / f"{e['id']}.txt").read_text().split()
        eligible = set(map(int, text))
        points = sv.read_points(pts_path)
        row = {k: e[k] for k in ("id", "kind", "carbon", "medium", "biological_replicate")}
        row["timepoint_h"] = int(e["timepoint_h"])
        row.update(sv.score_image(labels, points, eligible))
        row["precision"], row["recall"], row["f1"] = sv.metrics(
            row["true_positive"], row["false_positive"], row["false_negative"])
        row["count_error_pct"] = pct(row["automated_objects"], row["manual_cells"])
        # clicks on edge objects are left out, as in the counts
        kept = [p for p in points if sv.match_point(labels, *p) in eligible | {0}]
        ids = sorted(eligible)
        cents = ndimage.center_of_mass(labels > 0, labels, ids) if ids else []
        row["nn_manual_um"] = mean_nn_um(kept)
        row["nn_auto_um"] = mean_nn_um(cents)
        rows.append(row)
    return rows


def pool(name, g):
    tp, fp, fn = (sum(r[k] for r in g) for k in ("true_positive", "false_positive", "false_negative"))
    manual = sum(r["manual_cells"] for r in g)
    auto = sum(r["automated_objects"] for r in g)
    p, rc, f1 = sv.metrics(tp, fp, fn)
    errs = [r["count_error_pct"] for r in g if not math.isnan(r["count_error_pct"])]
    return {"group": name, "images": len(g), "manual_cells": manual, "automated_objects": auto,
            "true_positive": tp, "false_positive": fp, "false_negative": fn,
            "merged_objects": sum(r["merged_objects"] for r in g),
            "missed_cells": sum(r["missed_cells"] for r in g),
            "precision": p, "recall": rc, "f1": f1,
            "count_error_pct_pooled": pct(auto, manual),
            "count_error_pct_min": min(errs) if errs else float("nan"),
            "count_error_pct_max": max(errs) if errs else float("nan")}


def summarise_a(rows):
    groups = {"all": rows,
              "low density (whole stacks)": [r for r in rows if r["kind"] == "whole"],
              "high density (crops)": [r for r in rows if r["kind"] == "crop"],
              "+CAA": [r for r in rows if r["medium"] == "caa"],
              "-CAA": [r for r in rows if r["medium"] == "no_caa"]}
    for c in ("glucose", "glycerol", "succinate"):
        groups[c] = [r for r in rows if r["carbon"] == c]
        for kind, label in (("whole", "low density"), ("crop", "high density")):
            groups[f"{c}, {label}"] = [r for r in rows if r["carbon"] == c and r["kind"] == kind]
    for t in sorted({r["timepoint_h"] for r in rows}):
        groups[f"{t} h"] = [r for r in rows if r["timepoint_h"] == t]
    return [pool(n, g) for n, g in groups.items() if g]


def mu_bias(rows):
    """Slope of ln(auto/manual) vs time = error the segmentation adds to mu (h^-1)."""
    out = []
    win = [r for r in rows if r["medium"] == "caa" and MU_WINDOW[0] <= r["timepoint_h"] <= MU_WINDOW[1]
           and r["manual_cells"] > 0 and r["automated_objects"] > 0]
    for name in ("glucose", "glycerol", "succinate", "pooled"):
        g = win if name == "pooled" else [r for r in win if r["carbon"] == name]
        if len({r["timepoint_h"] for r in g}) < 3:
            continue
        t = np.array([r["timepoint_h"] for r in g], float)
        y = np.log([r["automated_objects"] / r["manual_cells"] for r in g])
        if name == "pooled":  # common slope, separate intercept per carbon
            carbons = sorted({r["carbon"] for r in g})
            X = np.column_stack([t] + [[r["carbon"] == c for r in g] for c in carbons]).astype(float)
        else:
            X = np.column_stack([t, np.ones_like(t)])
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        dof = len(y) - X.shape[1]
        se = math.sqrt(resid @ resid / dof * np.linalg.inv(X.T @ X)[0, 0]) if dof > 0 else float("nan")
        from scipy.stats import t as tdist
        q = tdist.ppf(0.975, dof) if dof > 0 else float("nan")
        out.append({"carbon": name, "images": len(g),
                    "mu_bias_per_h": round(beta[0], 4),
                    "ci95_low": round(beta[0] - q * se, 4), "ci95_high": round(beta[0] + q * se, 4)})
    return out


# --- Part B ----------------------------------------------------------------

def read_outlines(path, shape):
    """Fiji ROI zip -> list of 2D boolean masks (pixel-centre convention)."""
    masks = []
    for roi in roifile.roiread(path):
        coords = roi.coordinates(multi=True)
        mask = np.zeros(shape, bool)
        for xy in coords:
            if len(xy) >= 3:
                # ImageJ vertices sit on pixel corners; skimage uses centres
                mask ^= polygon2mask(shape, np.column_stack([xy[:, 1], xy[:, 0]]) - 0.5)
        if mask.any():
            masks.append(mask)
    return masks


def touches_border(mask):
    return mask[0].any() or mask[-1].any() or mask[:, 0].any() or mask[:, -1].any()


def shape_props(mask):
    p = regionprops(mask.astype(np.uint8))[0]
    elong = p.axis_major_length / p.axis_minor_length if p.axis_minor_length > 0 else np.inf
    return {"area": p.area, "centroid": np.array(p.centroid), "theta": p.orientation, "elong": elong}


def nematic_2d(thetas):
    """|<exp(2i theta)>|: 1 = all aligned, 0 = random (2D)."""
    if len(thetas) < 2:
        return float("nan")
    return round(float(abs(np.mean(np.exp(2j * np.asarray(thetas))))), 3)


def score_part_b(outline_dir):
    per_image, per_cell = [], []
    for e in read_csv(FINAL / "_manifest_trace.csv"):
        path = outline_dir / f"{e['id']}.zip"
        if not path.exists():
            continue
        lab = tifffile.imread(FINAL / "_labels" / f"{e['id']}.tif")
        manual = [m for m in read_outlines(path, lab.shape) if not touches_border(m)]
        auto = [lab == i for i in np.unique(lab) if i > 0]
        auto = [m for m in auto if not touches_border(m)]
        mprops = [shape_props(m) for m in manual]
        aprops = [shape_props(m) for m in auto]
        iou = np.zeros((len(manual), len(auto)))
        for i, m in enumerate(manual):
            for j, a in enumerate(auto):
                inter = np.logical_and(m, a).sum()
                if inter:
                    iou[i, j] = inter / np.logical_or(m, a).sum()
        rr, cc = linear_sum_assignment(-iou) if iou.size else ([], [])
        pairs = [(i, j, iou[i, j]) for i, j in zip(rr, cc) if iou[i, j] > 0]
        med_area = np.median([p["area"] for p in mprops]) if mprops else 0
        n_auto = sum(p["area"] >= MIN_FP_AREA_FRACTION * med_area for p in aprops)
        row = {"id": e["id"], **{k: e[k] for k in ("kind", "carbon", "medium", "timepoint_h")},
               "manual_cells": len(manual), "automated_objects": n_auto}
        for thr in IOU_THRESHOLDS:
            tp = sum(v >= thr for *_, v in pairs)
            fp = n_auto - sum(v >= thr and aprops[j]["area"] >= MIN_FP_AREA_FRACTION * med_area
                              for _, j, v in pairs)
            fn = len(manual) - tp
            row[f"f1_iou{int(thr * 100)}"] = round(2 * tp / (2 * tp + fp + fn), 3) if tp + fp + fn else float("nan")
        good = [(i, j, v) for i, j, v in pairs if v >= 0.5]
        row["matched_iou50"] = len(good)
        row["mean_iou_matched"] = round(float(np.mean([v for *_, v in good])), 3) if good else float("nan")
        row["area_ratio_median"] = round(float(np.median(
            [aprops[j]["area"] / mprops[i]["area"] for i, j, _ in good])), 3) if good else float("nan")
        row["foreground_area_ratio"] = round(
            float(np.sum([p["area"] for p in aprops]) / np.sum([p["area"] for p in mprops])), 3) if mprops else float("nan")
        el_m = [p["theta"] for p in mprops if p["elong"] >= MIN_ELONGATION]
        el_a = [p["theta"] for p in aprops if p["elong"] >= MIN_ELONGATION]
        row["nematic_2d_manual"], row["nematic_2d_auto"] = nematic_2d(el_m), nematic_2d(el_a)
        for i, j, v in good:
            m, a = mprops[i], aprops[j]
            both = m["elong"] >= MIN_ELONGATION and a["elong"] >= MIN_ELONGATION
            dtheta = abs(math.degrees(m["theta"] - a["theta"])) % 180
            per_cell.append({"id": e["id"], "carbon": e["carbon"], "medium": e["medium"],
                             "timepoint_h": e["timepoint_h"], "iou": round(v, 3),
                             "area_manual_um2": round(m["area"] * XY_UM ** 2, 3),
                             "area_auto_um2": round(a["area"] * XY_UM ** 2, 3),
                             "centroid_offset_um": round(float(np.linalg.norm(m["centroid"] - a["centroid"])) * XY_UM, 3),
                             "orientation_error_deg": round(min(dtheta, 180 - dtheta), 1) if both else ""})
        cells = [c for c in per_cell if c["id"] == e["id"]]
        ori = [c["orientation_error_deg"] for c in cells if c["orientation_error_deg"] != ""]
        row["centroid_offset_um_median"] = round(float(np.median([c["centroid_offset_um"] for c in cells])), 3) if cells else float("nan")
        row["orientation_error_deg_median"] = round(float(np.median(ori)), 1) if ori else float("nan")
        per_image.append(row)
    return per_image, per_cell


def summarise_b(per_image, per_cell):
    groups = {"all": (per_image, per_cell)}
    for kind, label in (("whole", "low density"), ("crop", "high density")):
        ids = {r["id"] for r in per_image if r["kind"] == kind}
        groups[label] = ([r for r in per_image if r["id"] in ids], [c for c in per_cell if c["id"] in ids])
    for c in ("glucose", "glycerol", "succinate"):
        groups[c] = ([r for r in per_image if r["carbon"] == c], [x for x in per_cell if x["carbon"] == c])
    out = []
    for name, (imgs, cells) in groups.items():
        if not imgs:
            continue
        ori = [c["orientation_error_deg"] for c in cells if c["orientation_error_deg"] != ""]
        row = {"group": name, "images": len(imgs),
               "manual_cells": sum(r["manual_cells"] for r in imgs),
               "matched_iou50": len(cells)}
        for thr in IOU_THRESHOLDS:
            k = f"f1_iou{int(thr * 100)}"
            row[k] = round(float(np.nanmean([r[k] for r in imgs])), 3)
        row.update({
            "mean_iou_matched": round(float(np.mean([c["iou"] for c in cells])), 3) if cells else float("nan"),
            "area_ratio_median": round(float(np.median([c["area_auto_um2"] / c["area_manual_um2"] for c in cells])), 3) if cells else float("nan"),
            "centroid_offset_um_median": round(float(np.median([c["centroid_offset_um"] for c in cells])), 3) if cells else float("nan"),
            "orientation_error_deg_median": round(float(np.median(ori)), 1) if ori else float("nan"),
            "orientation_cells": len(ori)})
        out.append(row)
    return out


# --- main ------------------------------------------------------------------

def run(points_dir, outline_dir, results_dir, quiet=False):
    a_rows = score_part_a(points_dir)
    a_sum = summarise_a(a_rows) if a_rows else []
    write_csv(results_dir / "partA_per_image.csv", a_rows)
    write_csv(results_dir / "partA_summary.csv", a_sum)
    bias = mu_bias(a_rows)
    write_csv(results_dir / "partA_mu_bias.csv", bias)
    b_img, b_cell = score_part_b(outline_dir)
    b_sum = summarise_b(b_img, b_cell) if b_img else []
    write_csv(results_dir / "partB_per_image.csv", b_img)
    write_csv(results_dir / "partB_per_cell.csv", b_cell)
    write_csv(results_dir / "partB_summary.csv", b_sum)
    rater2 = FINAL / "points_rater2"
    if rater2.exists() and any(rater2.glob("*.csv")) and points_dir == FINAL / "points":
        r2 = score_part_a(rater2)
        r1 = {r["id"]: r for r in a_rows}
        write_csv(results_dir / "rater_comparison.csv", [
            {"id": r["id"], "manual_cells_rater1": r1[r["id"]]["manual_cells"] if r["id"] in r1 else "",
             "manual_cells_rater2": r["manual_cells"],
             "precision_rater1": r1[r["id"]]["precision"] if r["id"] in r1 else "",
             "precision_rater2": r["precision"],
             "recall_rater1": r1[r["id"]]["recall"] if r["id"] in r1 else "",
             "recall_rater2": r["recall"]} for r in r2])
    if quiet:
        return a_sum, b_sum
    print(f"Part A: {len(a_rows)} images scored")
    for s in a_sum:
        print(f"  {s['group']:<30} n={s['images']:2d} manual={s['manual_cells']:4d} "
              f"auto={s['automated_objects']:4d} P={s['precision']:.3f} R={s['recall']:.3f} "
              f"F1={s['f1']:.3f} merged={s['merged_objects']:3d} count error={s['count_error_pct_pooled']:+.1f}%")
    for b in bias:
        print(f"  mu bias {b['carbon']:<9}: {b['mu_bias_per_h']:+.3f} h-1 "
              f"(95% CI {b['ci95_low']:+.3f} to {b['ci95_high']:+.3f})")
    print(f"Part B: {len(b_img)} slices scored")
    for s in b_sum:
        print(f"  {s['group']:<14} cells={s['manual_cells']:4d} F1@0.5={s['f1_iou50']:.3f} "
              f"mean IoU={s['mean_iou_matched']:.3f} area ratio={s['area_ratio_median']:.2f} "
              f"orientation error={s['orientation_error_deg_median']} deg")
    return a_sum, b_sum


def selftest():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        pts, outl = tmp / "points", tmp / "outlines"
        pts.mkdir(), outl.mkdir()
        for e in read_csv(FINAL / "_manifest.csv"):
            labels = tifffile.imread(FINAL / "_labels" / f"{e['id']}.tif")
            ids = map(int, (FINAL / "_eligible" / f"{e['id']}.txt").read_text().split())
            with (pts / f"{e['id']}.csv").open("w", newline="") as fh:
                fh.write(" ,X,Y,Slice\n")
                for n, obj in enumerate(ids, start=1):
                    zz, yy, xx = np.nonzero(labels == obj)
                    c = np.array([zz.mean(), yy.mean(), xx.mean()])
                    k = int(np.argmin((zz - c[0]) ** 2 + (yy - c[1]) ** 2 + (xx - c[2]) ** 2))
                    fh.write(f"{n},{xx[k] + 0.5},{yy[k] + 0.5},{zz[k] + 1}\n")
        for e in read_csv(FINAL / "_manifest_trace.csv"):
            lab = tifffile.imread(FINAL / "_labels" / f"{e['id']}.tif")
            rois = []
            for i in np.unique(lab):
                if i == 0:
                    continue
                padded = np.pad(lab == i, 1)
                contour = max(find_contours(padded.astype(float), 0.5), key=len) - 1
                rois.append(roifile.ImagejRoi.frompoints(contour[:, ::-1] + 0.5))
            if rois:
                roifile.roiwrite(outl / f"{e['id']}.zip", rois)
        a_sum, b_sum = run(pts, outl, tmp / "results", quiet=True)
    assert all(s["precision"] == 1 and s["recall"] == 1 for s in a_sum), a_sum
    allb = next(s for s in b_sum if s["group"] == "all")
    assert allb["f1_iou50"] >= 0.99 and allb["mean_iou_matched"] >= 0.85, allb
    print("Self-test passed. Part A P = R = 1;",
          f"Part B F1@0.5 = {allb['f1_iou50']}, mean IoU = {allb['mean_iou_matched']}",
          "(polygon rasterisation of label outlines, < 1 expected)")


if __name__ == "__main__":
    selftest() if "--selftest" in sys.argv else run(FINAL / "points", FINAL / "outlines", RESULTS)

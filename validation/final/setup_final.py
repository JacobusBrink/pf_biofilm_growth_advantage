"""Final segmentation validation (Reviewer comment on segmentation, Table S4).

Source: <stacks>/<carbon>/<caa|no_caa>/rep*/<t>h/*_stardist_opp_output/
(READ-ONLY), where <stacks> is the folder of re-segmented image stacks
(not in this repository). The StarDist label stack of a sample is the ch2 file in its
*_channel0 folder (object count = BiofilmQ Ncells in the matching
*_Pos1/data/txt_output/Summary*.csv); the raw stack is the ch1 file next to it.

Pre-specified, stratified, blinded design (fixed seed, written before any
annotation):

  Part A - cell counting (one click per cell, 3D)            F01 ... F30
    + CAA : glucose, glycerol, succinate x 2, 3, 4, 5, 6, 7, 8 h   (21)
    - CAA : glucose, glycerol, succinate x 2, 5, 8 h                 (9)
    Biological replicates rotate across timepoints, technical
    replicate random. Low-density stacks (<= MAX_WHOLE_CELLS objects) are
    annotated whole; denser stacks are cropped to ~30 cells in a 9-slice
    band around the densest layer (worst case for merging; same rule as
    rounds 1 and 2).

  Part B - cell outlines (2D, one slice)                       B01 ... B12
    + CAA 3 h, 5 h, 8 h and - CAA 8 h, for each carbon: the same samples as
    Part A, one slice (the densest slice of the crop / stack). Gives IoU,
    cell area, centroid position and orientation (biovolume, spacing and
    nematic order checks).

Technical replicate 1 (round 1) and the round-2 glucose samples are excluded,
so every image is new to the annotator.

Writes COPIES into validation/final/:
  images/Fxx.tif        raw stacks to click (Part A)
  trace/Bxx.tif         raw single slices to outline (Part B)
  _labels/Fxx.tif, _labels/Bxx.tif   StarDist labels (do not open while annotating)
  _eligible/Fxx.txt     object ids counted in each crop (>= half inside)
  _manifest.csv, _manifest_trace.csv   key (do not open until finished)

Run from the repository root:  python validation/final/setup_final.py <stacks>
Source paths in the manifests are written relative to <stacks>.
"""

import csv
import random
import re
import sys
from pathlib import Path

import numpy as np
import tifffile

sys.path.insert(0, str(Path("validation")))
from setup_validation import dense_crop_box  # same crop rule as rounds 1 and 2

SRC = None  # set from the command line in main()
OUT = Path("validation/final")

TARGET_CELLS_PER_CROP = 30
MAX_WHOLE_CELLS = 50
MIN_FRACTION_IN_CROP = 0.5
SEED = 20261010

CARBONS = ["glucose", "glycerol", "succinate"]
CAA_TIMES = [2, 3, 4, 5, 6, 7, 8]
NOCAA_TIMES = [2, 5, 8]
TRACE_STRATA = {("caa", 3), ("caa", 5), ("caa", 8), ("no_caa", 8)}

LABEL_RE = re.compile(r"^([a-z]+?)(caa|control)(\d+)hbr(\d+)tr(\d+)_pos1_ch2_time1\.tif$")


def index_label_files():
    index = {}
    for path in SRC.glob("*/*/rep*/*h/*_stardist_opp_output/*_channel0/*_ch2_time1.tif"):
        m = LABEL_RE.match(path.name)
        if m:
            carbon, medium = m.group(1), "caa" if m.group(2) == "caa" else "no_caa"
            t, br, tr = map(int, m.groups()[2:])
            index[(carbon, medium, t, br, tr)] = path
    return index


def excluded():
    """Samples already seen by the annotator (round 1: tr1; round 2 list)."""
    seen = set()
    with Path("validation/round2/_manifest.csv").open(newline="") as fh:
        for r in csv.DictReader(fh):
            seen.add(("glucose", "caa", int(r["timepoint_h"]),
                      int(r["biological_replicate"]), int(r["technical_replicate"])))
    return seen


def biofilmq_ncells(label_path):
    summary = next(label_path.parent.parent.glob(
        f"{label_path.parent.name}_Pos1/data/txt_output/Summary_*.csv"))
    with summary.open(newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["Ncells"]]
    return int(float(rows[0]["Ncells"]))


def build_plan():
    strata = [(c, "caa", t) for c in CARBONS for t in CAA_TIMES]
    strata += [(c, "no_caa", t) for c in CARBONS for t in NOCAA_TIMES]
    # rotate biological replicates so each carbon x medium uses all three
    plan = []
    for i, (carbon, medium, t) in enumerate(strata):
        offset = CARBONS.index(carbon)
        plan.append((carbon, medium, t, (i + offset) % 3 + 1))
    return plan


def pick(index, seen, carbon, medium, t, br, rng):
    trs = sorted(tr for (c, m, tt, b, tr) in index
                 if (c, m, tt, b) == (carbon, medium, t, br) and tr != 1
                 and (c, m, tt, b, tr) not in seen)
    rng.shuffle(trs)
    for tr in trs:
        path = index[(carbon, medium, t, br, tr)]
        ncells = biofilmq_ncells(path)
        if ncells == 0:
            continue
        labels = tifffile.imread(path)
        assert len(np.unique(labels)) - 1 == ncells, path
        return tr, ncells, path, labels
    raise RuntimeError(f"no sample for {carbon} {medium} {t} h br{br}")


def eligible_ids(lab_crop, labels_full):
    total = np.bincount(labels_full.ravel())
    inside = np.bincount(lab_crop.ravel(), minlength=len(total))
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = inside / total
    return [int(i) for i in np.nonzero(frac >= MIN_FRACTION_IN_CROP)[0] if i > 0]


def write_tif(path, array):
    path.parent.mkdir(parents=True, exist_ok=True)
    axes = "ZYX" if array.ndim == 3 else "YX"
    tifffile.imwrite(path, array, imagej=True, metadata={"axes": axes})


def write_csv(path, rows):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    global SRC
    if len(sys.argv) != 2:
        sys.exit("usage: python validation/final/setup_final.py <stacks folder>")
    SRC = Path(sys.argv[1])
    if (OUT / "_manifest.csv").exists():
        sys.exit("validation/final/_manifest.csv exists: the set is already made "
                 "(delete it only if no annotation has started).")
    index = index_label_files()
    seen = excluded()
    rng = random.Random(SEED)
    plan = build_plan()
    rng.shuffle(plan)

    manifest, trace = [], []
    for n, (carbon, medium, t, br) in enumerate(plan, start=1):
        vid = f"F{n:02d}"
        tr, ncells, label_path, labels = pick(index, seen, carbon, medium, t, br, rng)
        raw_path = label_path.with_name(label_path.name.replace("_ch2_", "_ch1_"))
        raw = tifffile.imread(raw_path)
        if ncells <= MAX_WHOLE_CELLS:
            kind = "whole"
            z0, z1, y0, y1, x0, x1 = 0, labels.shape[0], 0, labels.shape[1], 0, labels.shape[2]
        else:
            kind = "crop"
            z0, z1, y0, y1, x0, x1 = dense_crop_box(labels, TARGET_CELLS_PER_CROP)
        raw_crop = raw[z0:z1, y0:y1, x0:x1]
        lab_crop = labels[z0:z1, y0:y1, x0:x1]
        eligible = eligible_ids(lab_crop, labels)
        write_tif(OUT / "images" / f"{vid}.tif", raw_crop)
        write_tif(OUT / "_labels" / f"{vid}.tif", lab_crop)
        (OUT / "_eligible").mkdir(exist_ok=True)
        (OUT / "_eligible" / f"{vid}.txt").write_text("\n".join(map(str, eligible)))
        entry = {
            "id": vid, "kind": kind, "carbon": carbon, "medium": medium,
            "timepoint_h": t, "biological_replicate": br, "technical_replicate": tr,
            "biofilmq_ncells": ncells, "objects_in_image": len(eligible),
            "z0": z0, "z1": z1, "y0": y0, "y1": y1, "x0": x0, "x1": x1,
            "raw_source": raw_path.relative_to(SRC).as_posix(),
            "label_source": label_path.relative_to(SRC).as_posix(),
        }
        manifest.append(entry)
        print(vid, kind, carbon, medium, f"{t}h", f"br{br}tr{tr}",
              "Ncells", ncells, "objects in image:", len(eligible))

        if (medium, t) in TRACE_STRATA:
            # densest slice of the image (most foreground in the labels)
            z = int(np.argmax((lab_crop > 0).sum(axis=(1, 2))))
            trace.append({"source_id": vid, "slice_in_image": z + 1,
                          "raw": raw_crop[z], "labels": lab_crop[z], **{
                              k: entry[k] for k in ("kind", "carbon", "medium", "timepoint_h",
                                                    "biological_replicate", "technical_replicate")}})

    write_csv(OUT / "_manifest.csv", manifest)

    rng.shuffle(trace)
    rows = []
    for n, e in enumerate(trace, start=1):
        bid = f"B{n:02d}"
        write_tif(OUT / "trace" / f"{bid}.tif", e.pop("raw"))
        write_tif(OUT / "_labels" / f"{bid}.tif", e.pop("labels"))
        rows.append({"id": bid, **e})
    write_csv(OUT / "_manifest_trace.csv", rows)
    for d in ("points", "outlines"):
        (OUT / d).mkdir(exist_ok=True)
    print(len(manifest), "Part A images,",
          sum(int(r["objects_in_image"]) for r in manifest), "objects;",
          len(rows), "Part B slices")


if __name__ == "__main__":
    main()

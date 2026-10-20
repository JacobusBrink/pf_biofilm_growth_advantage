# Final segmentation validation (Table S4, comment on segmentation)

This answers the request for validation against manual ground truth, using
precision/recall, IoU, manual count comparison and representative overlays, across carbon
sources, time points, and low- and high-density biofilms. The design was fixed (seed
20261010) before any annotation was done.

The source is the folder of re-segmented StarDist + BiofilmQ image stacks (not in this repository), passed to `setup_final.py`. Nothing there is changed.
Every image is new to you: technical replicate 1 (round 1) and the round-2 glucose samples are excluded.

## What is in the set

| Part | Images | What you do | What it measures | Time |
|---|---|---|---|---|
| A | `images/F01.tif` … `F30.tif` (627 objects) | One click per cell (3D) | Precision, recall, F1, merged cells, manual vs automated count, bias in µ, nearest-neighbour distance | ~1.5 h |
| B | `trace/B01.tif` … `B12.tif` (~310 cells) | Trace the outline of every cell (one slice) | IoU, cell area (biovolume), centroid position (spacing), orientation (nematic order) | ~2.5 h |

**Part A** covers glucose, glycerol and succinate. With CAA it has every hour from 2 to 8 h (21 images); without CAA it has 2, 5 and 8 h (9 images).
Stacks with ≤ 50 cells are whole (low density). Denser stacks are cropped to ~30 cells in a 9-slice band around the densest layer (high density, the worst case for merging). This is the same rule as rounds 1 and 2.
Biological replicates rotate.

**Part B** uses one slice of 12 of those samples: 3 h, 5 h and 8 h with CAA, and 8 h without CAA, for each carbon.

**Decided 2026-10-10 (time):** do all of Part A, and in Part B trace only **B01, B06 and B11**
(the 8 h + CAA slices, one per carbon, ~100 cells, ~45-60 min). All architecture metrics in the
paper are measured at 8 h, so these slices test the masks they come from. The other B images can
be ignored; the script skips slices without an outline file.

Do **not** open `_labels/`, `_eligible/`, `_manifest.csv` or `_manifest_trace.csv` until you have finished.
Do Part A first, then Part B.

## Part A: clicking cells (same protocol as rounds 1 and 2)

One-time Fiji setup: go to **Analyze › Set Measurements**, tick **Stack position**, then click OK. Then set the Multi-point tool options: use a visible colour, tick **Label points**, and tick **Show on all slices** if your version has it.

For each `images/Fxx.tif`:

1. Open the image. Set **Image › Adjust › Brightness/Contrast › Auto** (display only, do not press Apply).
2. With the **Multi-point** tool, scroll through z and **click once on every cell, on its sharpest slice**.
   Touching cells get one click each. Click cells cut by the crop edge or by the top/bottom slice
   too; the script ignores cells that are mostly outside the crop. To delete a wrong click, use **Alt + click**.
3. Press **Ctrl + M**. In the Results window, choose **File › Save As** → `points/Fxx.csv`. Then clear the Results before the next image.

Some low-density images (2 h) contain only 2–5 cells. That is expected; annotate them anyway.

## Part B: tracing cell outlines

For each `trace/Bxx.tif` (a single slice):

1. Open the image and set **Auto** contrast (display only). Zoom in to 400–800%.
2. Open the ROI Manager: **Analyze › Tools › ROI Manager**.
3. With the **Freehand** (or **Polygon**) selection tool, trace the outline of a cell along the edge
   of its fluorescence, then press **T** to add it to the ROI Manager. Repeat for **every cell
   visible in focus on this slice**, including cells cut by the image edge (the script ignores
   any outline that touches the edge). Skip faint, out-of-focus blurs from cells above or below.
   Touching cells get separate outlines.
4. In the ROI Manager, choose **More › Save…** → `outlines/Bxx.zip`. Then **Delete** all ROIs before the next image.

Trace by what you see in the raw image only.

## Optional: a second annotator

If a colleague clicks some of the Part A images (blind, same protocol) into `points_rater2/Fxx.csv`,
the script reports both annotators side by side (`results/rater_comparison.csv`). Even 4–6 images
(2 low, 2–4 high density) show how much of any disagreement is annotator variability.

## Scoring and figures

From the repository root:

```
python validation/final/score_final.py
python validation/final/figures_final.py
```

You can score partway through; images without annotation files are skipped.
`python validation/final/score_final.py --selftest` checks the scoring against synthetic perfect annotations.

Outputs (`validation/final/results/`):

- `partA_per_image.csv`, `partA_summary.csv`: precision, recall, F1, merged and missed cells, and count error.
  These are given overall, by density, with and without CAA, per carbon × density, and per timepoint.
  Also the mean nearest-neighbour distance from the clicks vs from the segmented centroids.
- `partA_mu_bias.csv`: the slope of ln(automated / manual count) against time over 3–7 h, with CAA, per carbon and pooled.
  This is the error, in h⁻¹, that segmentation adds to the specific growth rate. A slope near 0 means
  any count error is constant over time and cancels out of µ.
- `partB_per_image.csv`, `partB_per_cell.csv`, `partB_summary.csv`: objects matched one-to-one by IoU (Hungarian assignment), with F1 at IoU 0.5–0.9 and mean IoU of matched cells.
  Also the automated/manual area ratio (biovolume bias), centroid offset (µm), orientation error of elongated cells (axis ratio ≥ 1.5), and the 2D nematic order parameter from manual vs automated cells.
- Figures in `output/segmentation_validation/`: `Figure_S_seg_overlays.png` (representative overlays, coloured by outcome once scored) and `Figure_S_seg_metrics.png`.

## Scoring rules

- **Part A** uses the round-1 matching rules, imported from `validation/score_validation.py` (unchanged):
  - A click inside an object, or within 3 px in xy and 1 slice in z, matches that object.
  - An object hit by 2+ clicks is *merged*: 1 true positive, and the extra clicks are false negatives.
  - A click matching no object is *missed*.
  - An object with no click is a *false positive*. This can be a split fragment or a cell that was not clicked; the overlays show which.
  - In crops, objects with < 50% of their volume inside the crop are ignored, together with their clicks.
- **Part B**: outlines and label cross-sections that touch the image border are ignored.
  Automated cross-sections smaller than 25% of the median traced cell area (cell tips grazing the slice) count only if they match a traced cell.
  The voxel size used is 0.0916 µm (xy) × 0.31 µm (z).

## Interpretation, and limits to state in the paper

- Clicks validate detection and counts in 3D. Outlines validate shape in 2D on one slice, so IoU and orientation are 2D measures of the 3D masks.
- The ground truth comes from a single annotator, unless the second-annotator option is used.
- Dense crops sample the densest layer. Accuracy elsewhere in the biofilm is likely as good or better.

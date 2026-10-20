# Segmentation validation: annotation protocol

You click once on every cell in 18 small image stacks
(about 260–350 clicks, 30–45 min). A script then compares your clicks with the
StarDist segmentation.

The images are copies. Nothing in `E:\` has been changed.

## Folder contents

| Path | What it is |
|---|---|
| `images/V01.tif` … `V18.tif` | Raw image stacks to annotate. Names are shuffled, so you do not know the condition while clicking |
| `points/` | Where you save your clicks (one CSV per image) |
| `_labels/` | StarDist segmentations. **Do not open these while annotating** |
| `_manifest.csv` | Which sample each V-number is. **Do not open until you have finished** |
| `setup_validation.py` | Made this folder (already run) |
| `score_validation.py` | Scores your clicks against the segmentation |

## One-time Fiji setup

1. **Analyze › Set Measurements…** and tick **Stack position**. The other boxes do not matter. Click OK.
   Without this, the Slice number is not saved and the script will stop with an error.
2. Double-click the **Multi-point** tool icon. In the options, set a clearly visible colour and size, and tick **Label points**.
   If your version has **Show on all slices**, tick it too. It helps you see which cells you have already clicked in neighbouring slices.

## For each image (V01 to V18)

1. **File › Open** `images/Vxx.tif`.
2. **Image › Adjust › Brightness/Contrast › Auto**, for display only. Do not press Apply.
3. Select the **Multi-point** tool and scroll through the z-slices.
4. **Click once on every cell**, on the slice where it looks sharpest.
   - One click per cell, even if the cell is visible across several slices.
   - Two touching cells count as two cells, so give each its own click.
   - Dense crops (V03, V04, V09, V11, V15, V16) are thin, only 9 slices, taken from the middle of a biofilm. Click every cell you can clearly see, including cells cut off by the image edges or by the top or bottom slice. The script automatically ignores cells that are mostly outside the crop, so you do not need to judge this yourself.
   - A wrong click can be deleted with **Alt + click** on the point.
5. Press **Ctrl + M** (Analyze › Measure). The Results window shows one row per click, with X, Y and Slice.
6. In the Results window: **File › Save As…** and save as `points/Vxx.csv`, using the same number as the image.
7. Close the image. In the Results window, choose **Edit › Clear Results** (or close the window) before the next image.

Annotate by what you see in the raw image. Do not look at the segmentation first.

## When finished

Tell Claude, or run from the repository root:

```
python validation/score_validation.py
```

Results are written to `results/per_image.csv` and `results/summary.csv`: precision, recall, F1, merged objects, and the percentage difference between manual and automated counts, for low-density stacks (2–4 h) and dense crops (6–8 h).

## What was selected

- **Low density:** whole stacks at 2 h and 4 h for glucose, glycerol and succinate, with and without CAA (12 stacks).
- **High density:** crops of about 30 cells (a 9-slice band around the densest layer) from 6 h and 8 h biofilms of glucose, glycerol and succinate with CAA (6 crops).
- Biological replicates rotate across samples.
- The segmentation used for each sample is the one whose object count equals the BiofilmQ cell count used in the paper.

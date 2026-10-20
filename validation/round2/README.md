# Round 2: preliminary validation of the new segmentation

A quick check of the re-segmented glucose + CAA stacks (`D:\validation_stacks`, 2–7 h) before
re-segmenting the other stacks. There are 10 blinded images, about 470 cells in total. Allow
about an hour.

Nothing in `D:\validation_stacks` or `E:\` has been changed. Round 1 (`validation/`) is untouched.

## Folder contents

| Path | What it is |
|---|---|
| `images/W01.tif` … `W10.tif` | Raw stacks to annotate. The names are shuffled |
| `points/` | Save your clicks here (`W01.csv` etc.) |
| `_labels/` | New StarDist segmentations. **Do not open these while annotating** |
| `_manifest.csv` | Which sample each W-number is. **Do not open it until you have finished** |
| `rescore/` | Round-1 stacks V07 and V09 with the new labels, scored with your round-1 clicks. You do not need to annotate these |
| `setup_round2.py` / `score_round2.py` | Made this folder / scores it |

## What is in it

- Whole stacks at 2 h, 3 h and 4 h (one each).
- Crops of about 30 cells (a 9-slice band around the densest layer, the same rule as round 1): one at 4 h and two each at 5 h, 6 h and 7 h, from different biological replicates.
- Technical replicate 1 was excluded, because round 1 used tr1 stacks.

## Annotating

Use the same Fiji protocol as round 1 (`validation/README.md`): open `images/Wxx.tif` and set
Auto contrast for display only. Click once per cell on its sharpest slice, and give touching cells
one click each. Click cells cut by the crop edge too; the script ignores cells that are mostly
outside the crop. Then press **Ctrl + M** and save the Results as `points/Wxx.csv`. Clear the
Results before the next image.

The crops are denser than in round 1 (44–78 objects each, because the new segmentation finds
more cells). Take your time on the 6 h and 7 h crops.

## Scoring

From the repository root:

```
python validation/round2/score_round2.py
```

The results go to `results/summary.csv`, which gives precision, recall, merged objects and count
error per timepoint, for 2–4 h, for 5–7 h and overall. The file `results/rescore_per_image.csv`
compares the old and new segmentation on V07 and V09.

## What would count as good enough

- Round 1 (old segmentation), dense crops: precision 0.96, recall 0.62, counts 23–44% too low.
- The goal is higher recall **with precision still high**. If precision falls, cells are being
  split, so counts rise for the wrong reason.
- First sign from V09 (6 h, round-1 clicks): recall 0.65 → 0.94, but precision 0.93 → 0.84 and
  the count error −30% → +12%. So watch for splitting.

"""Table S4: segmentation validation summary (from results/partA_summary.csv).

Run from the repository root:  python validation/final/table_s4.py
Output: output/segmentation_validation/Table_S4_segmentation_validation.csv
"""

import csv
from pathlib import Path

RES = Path("validation/final/results")
OUT = Path("output/segmentation_validation/Table_S4_segmentation_validation.csv")
GROUPS = [("all", "All images"),
          ("low density (whole stacks)", "Low density (whole stacks)"),
          ("high density (crops)", "High density (crops)"),
          ("glucose", "Glucose"), ("glycerol", "Glycerol"), ("succinate", "Succinate"),
          ("+CAA", "With casamino acids"), ("-CAA", "Without casamino acids")]

with (RES / "partA_summary.csv").open(newline="") as fh:
    summary = {r["group"]: r for r in csv.DictReader(fh)}
rows = [{"Group": label, "Images": summary[g]["images"],
         "Manual cells": summary[g]["manual_cells"],
         "Automated objects": summary[g]["automated_objects"],
         "Precision": f"{float(summary[g]['precision']):.2f}",
         "Recall": f"{float(summary[g]['recall']):.2f}",
         "F1": f"{float(summary[g]['f1']):.2f}",
         "Count error (%)": f"{float(summary[g]['count_error_pct_pooled']):+.1f}"}
        for g, label in GROUPS]
OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="", encoding="utf-8-sig") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
for r in rows:
    print(r)

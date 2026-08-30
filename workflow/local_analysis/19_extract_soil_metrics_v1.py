#!/usr/bin/env python3
import csv
import hashlib
import sys
from pathlib import Path

import openpyxl

PROJECT = Path(__file__).resolve().parents[2]
SOURCE = PROJECT / "羊肚菌连作土壤土壤8项检测-C2604129吕锐玲9个土壤数据-汇总(1).xlsx"
OUTDIR = PROJECT / "analysis" / "downstream" / "ITS_v5_soil_integration" / "input"
DATA_OUT = OUTDIR / "soil_metrics_normalized_v1.tsv"
DICT_OUT = OUTDIR / "soil_metric_dictionary_v1.tsv"
PROV_OUT = OUTDIR / "soil_extraction_provenance_v1.tsv"

for target in (DATA_OUT, DICT_OUT, PROV_OUT):
    if target.exists():
        raise SystemExit(f"Refusing to overwrite existing output: {target}")
if not SOURCE.is_file():
    raise SystemExit(f"Missing source workbook: {SOURCE}")

wb = openpyxl.load_workbook(SOURCE, data_only=True, read_only=True)
if wb.sheetnames != ["Sheet1"]:
    raise SystemExit(f"Unexpected worksheets: {wb.sheetnames}")
ws = wb["Sheet1"]
if ws.max_row != 22 or ws.max_column != 10:
    raise SystemExit(f"Unexpected worksheet dimensions: {ws.max_row}x{ws.max_column}")

metrics = [
    ("C", "pH", "pH", "confirmed_from_readable_header"),
    ("D", "metric_D_total_category_1", "g/kg", "formal_analyte_name_pending_confirmation"),
    ("E", "metric_E_total_category_2", "g/kg", "formal_analyte_name_pending_confirmation"),
    ("F", "metric_F_total_category_3", "g/kg", "formal_analyte_name_pending_confirmation"),
    ("G", "metric_G_available_category_1", "g/kg", "formal_analyte_name_pending_confirmation"),
    ("H", "metric_H_available_category_2", "mg/kg", "formal_analyte_name_pending_confirmation"),
    ("I", "metric_I_available_category_3", "mg/kg", "formal_analyte_name_pending_confirmation"),
    ("J", "organic_matter", "g/kg", "confirmed_from_header_pattern_and_project_plan"),
]

blocks = [(3, 11, "after_cultivation"), (14, 22, "before_cultivation")]
rows = []
for start, end, phase in blocks:
    for r in range(start, end + 1):
        lab_id = ws[f"A{r}"].value
        sample_id = ws[f"B{r}"].value
        if not isinstance(lab_id, str) or not isinstance(sample_id, str):
            raise SystemExit(f"Missing identifiers at row {r}")
        expected_prefix = "A" if phase == "after_cultivation" else "M"
        if not sample_id.startswith(expected_prefix):
            raise SystemExit(f"Unexpected phase/sample mapping at row {r}: {sample_id}")
        record = {"sample_id": sample_id, "phase": phase, "lab_sample_id": lab_id}
        for col, name, _, _ in metrics:
            value = ws[f"{col}{r}"].value
            if not isinstance(value, (int, float)):
                raise SystemExit(f"Non-numeric soil value at {col}{r}: {value!r}")
            record[name] = value
        rows.append(record)

sample_ids = {x["sample_id"] for x in rows}
expected = {f"{phase}{year}-{rep}" for phase in ("A", "M")
            for year in (2, 3, 4) for rep in (1, 2, 3)}
if sample_ids != expected or len(rows) != 18:
    raise SystemExit("Soil sample set does not match the expected 18 paired M/A samples")

OUTDIR.mkdir(parents=True, exist_ok=False)
fields = ["sample_id", "phase", "lab_sample_id"] + [m[1] for m in metrics]
with DATA_OUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

with DICT_OUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["source_column", "derived_name", "unit", "name_status"])
    writer.writerows(metrics)

sha256 = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
with PROV_OUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["field", "value"])
    writer.writerow(["source_workbook", SOURCE.name])
    writer.writerow(["source_sha256", sha256])
    writer.writerow(["worksheet", "Sheet1"])
    writer.writerow(["after_rows", "3-11"])
    writer.writerow(["before_rows", "14-22"])
    writer.writerow(["extracted_samples", len(rows)])
    writer.writerow(["formula_cells", 0])

print(f"SOIL_EXTRACTION_OK samples={len(rows)} metrics={len(metrics)} sha256={sha256}")


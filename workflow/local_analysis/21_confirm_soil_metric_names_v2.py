#!/usr/bin/env python3
import csv
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
OLD = PROJECT / "analysis" / "downstream" / "ITS_v5_soil_integration" / "input"
NEW = PROJECT / "analysis" / "downstream" / "ITS_v5_soil_integration_v2" / "input"
SOURCE = OLD / "soil_metrics_normalized_v1.tsv"
OLD_PROV = OLD / "soil_extraction_provenance_v1.tsv"
DATA_OUT = NEW / "soil_metrics_normalized_v2.tsv"
DICT_OUT = NEW / "soil_metric_dictionary_v2.tsv"
PROV_OUT = NEW / "soil_name_confirmation_provenance_v2.tsv"

for target in (DATA_OUT, DICT_OUT, PROV_OUT):
    if target.exists():
        raise SystemExit(f"Refusing to overwrite: {target}")
if not SOURCE.is_file() or not OLD_PROV.is_file():
    raise SystemExit("Missing v1 soil extraction inputs")

mapping = [
    ("pH", "pH", "pH", "confirmed_from_workbook"),
    ("metric_D_total_category_1", "total_nitrogen", "g/kg", "confirmed_by_user_2026-08-21"),
    ("metric_E_total_category_2", "total_phosphorus", "g/kg", "confirmed_by_user_2026-08-21"),
    ("metric_F_total_category_3", "total_potassium", "g/kg", "confirmed_by_user_2026-08-21"),
    ("metric_G_available_category_1", "available_nitrogen", "g/kg", "confirmed_by_user_2026-08-21"),
    ("metric_H_available_category_2", "available_phosphorus", "mg/kg", "confirmed_by_user_2026-08-21"),
    ("metric_I_available_category_3", "available_potassium", "mg/kg", "confirmed_by_user_2026-08-21"),
    ("organic_matter", "organic_matter", "g/kg", "confirmed_from_workbook_and_project_plan"),
]

with SOURCE.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if len(rows) != 18:
    raise SystemExit(f"Expected 18 rows, found {len(rows)}")

NEW.mkdir(parents=True, exist_ok=False)
fields = ["sample_id", "phase", "lab_sample_id"] + [new for _, new, _, _ in mapping]
with DATA_OUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        out = {k: row[k] for k in ("sample_id", "phase", "lab_sample_id")}
        for old, new, _, _ in mapping:
            out[new] = row[old]
        writer.writerow(out)

with DICT_OUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["v1_name", "confirmed_name", "unit", "name_status"])
    writer.writerows(mapping)

with PROV_OUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["field", "value"])
    writer.writerow(["source_numeric_table", str(SOURCE.relative_to(PROJECT))])
    writer.writerow(["source_extraction_provenance", str(OLD_PROV.relative_to(PROJECT))])
    writer.writerow(["name_confirmation_source", "user_message"])
    writer.writerow(["name_confirmation_date", "2026-08-21"])
    writer.writerow(["yield_data_status", "not_currently_available"])
    writer.writerow(["disease_numeric_data_status", "not_currently_available"])

print(f"SOIL_NAME_CONFIRMATION_V2_OK samples={len(rows)} metrics={len(mapping)}")


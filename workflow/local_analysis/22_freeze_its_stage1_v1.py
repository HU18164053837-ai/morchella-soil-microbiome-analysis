#!/usr/bin/env python3
import csv
import hashlib
import shutil
import zipfile
from datetime import datetime, timezone
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
FREEZE_ROOT = PROJECT / "analysis" / "frozen"
NAME = "ITS_STAGE1_FREEZE_2026-08-21_v1"
DEST = FREEZE_ROOT / NAME
ZIP_PATH = FREEZE_ROOT / f"{NAME}.zip"
ZIP_SHA = FREEZE_ROOT / f"{NAME}.zip.sha256"

for target in (DEST, ZIP_PATH, ZIP_SHA):
    if target.exists():
        raise SystemExit(f"Refusing to overwrite frozen artifact: {target}")

source_dirs = [
    Path("analysis/hpc_results/ITS_dada2_v1_corrected_v1"),
    Path("analysis/hpc_results/ITS_UNITE_v2_job119642037"),
    Path("analysis/downstream/ITS_v1"),
    Path("analysis/downstream/ITS_v2_taxonomy"),
    Path("analysis/downstream/ITS_v3_design_aware"),
    Path("analysis/downstream/ITS_v4_genus_screening"),
    Path("analysis/downstream/ITS_v5_soil_integration_v2"),
]
single_files = [
    Path("analysis/metadata/sample_metadata.tsv"),
    Path("analysis/reports/ITS_STAGE1_TECHNICAL_SUMMARY_ZH_v1.md"),
]

script_dir = PROJECT / "analysis" / "scripts"
for path in sorted(script_dir.iterdir()):
    if path.is_file() and path.name[:2].isdigit() and 11 <= int(path.name[:2]) <= 22:
        single_files.append(path.relative_to(PROJECT))

slurm_dir = PROJECT / "analysis" / "slurm"
if slurm_dir.is_dir():
    for path in sorted(slurm_dir.iterdir()):
        if path.is_file() and "its" in path.name.lower():
            single_files.append(path.relative_to(PROJECT))

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

selected = []
for rel_dir in source_dirs:
    src_dir = PROJECT / rel_dir
    if not src_dir.is_dir():
        raise SystemExit(f"Missing required source directory: {src_dir}")
    for src in sorted(src_dir.rglob("*")):
        if src.is_file():
            selected.append((src, Path("evidence") / src.relative_to(PROJECT), "analysis_evidence"))
for rel in single_files:
    src = PROJECT / rel
    if not src.is_file():
        raise SystemExit(f"Missing required source file: {src}")
    category = "technical_report" if "reports" in rel.parts else "metadata_or_code"
    selected.append((src, Path("evidence") / rel, category))

dedup = {}
for src, rel, category in selected:
    dedup[str(rel).lower()] = (src, rel, category)
selected = sorted(dedup.values(), key=lambda x: str(x[1]).lower())

DEST.mkdir(parents=True, exist_ok=False)
manifest_rows = []
for src, rel, category in selected:
    dst = DEST / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    source_hash = sha256(src)
    frozen_hash = sha256(dst)
    if source_hash != frozen_hash:
        raise SystemExit(f"Hash mismatch after copy: {src}")
    manifest_rows.append({
        "frozen_relative_path": rel.as_posix(),
        "source_relative_path": src.relative_to(PROJECT).as_posix(),
        "category": category,
        "bytes": dst.stat().st_size,
        "sha256": frozen_hash,
        "copy_verified": "yes",
    })

manifest = DEST / "FREEZE_MANIFEST_SHA256.tsv"
with manifest.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]),
                            delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(manifest_rows)

metadata = DEST / "FREEZE_METADATA.tsv"
with metadata.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["field", "value"])
    writer.writerow(["freeze_name", NAME])
    writer.writerow(["freeze_created_utc", datetime.now(timezone.utc).isoformat()])
    writer.writerow(["project", "Morchella soil microbiome"])
    writer.writerow(["marker", "ITS"])
    writer.writerow(["status", "stage1_analysis_complete"])
    writer.writerow(["frozen_file_count", len(manifest_rows)])
    writer.writerow(["frozen_bytes", sum(int(r["bytes"]) for r in manifest_rows)])
    writer.writerow(["sample_count", 24])
    writer.writerow(["asv_count_nochim", 8056])
    writer.writerow(["nonchim_reads", 942492])
    writer.writerow(["yield_data", "not_available"])
    writer.writerow(["disease_severity_numeric_data", "not_available"])

readme = DEST / "FROZEN_DO_NOT_EDIT.md"
readme.write_text(
    "# ITS Stage 1 frozen release\n\n"
    "This directory is a content-hashed snapshot of the completed ITS analysis as of 2026-08-21.\n\n"
    "- Do not edit files in place.\n"
    "- New analyses must create a new versioned directory.\n"
    "- Verify files against `FREEZE_MANIFEST_SHA256.tsv`.\n"
    "- The reader-facing technical summary is under `evidence/analysis/reports/`.\n"
    "- Raw FASTQ files and failed-attempt directories are intentionally excluded.\n",
    encoding="utf-8")

exclusions = DEST / "EXCLUSIONS_AND_OPEN_ITEMS.md"
exclusions.write_text(
    "# Exclusions and open items\n\n"
    "## Intentionally excluded\n\n"
    "- Raw and primer-trimmed FASTQ files: preserved in the project but not duplicated in this release.\n"
    "- Directories explicitly named `failed_attempt` or `preliminary_thresholds`.\n"
    "- The ambiguous-name soil integration v1; formal confirmed-name v2 is included.\n"
    "- 16S results, because this is an ITS-only stage freeze.\n\n"
    "## Open items\n\n"
    "- Yield values are not currently available.\n"
    "- Numeric disease incidence or severity is not currently available.\n"
    "- Cross-domain 16S-ITS integration awaits completion of the 16S workflow.\n",
    encoding="utf-8")

for row in manifest_rows:
    path = DEST / row["frozen_relative_path"]
    if path.stat().st_size != int(row["bytes"]) or sha256(path) != row["sha256"]:
        raise SystemExit(f"Final frozen verification failed: {path}")

with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED,
                     compresslevel=6, allowZip64=True) as archive:
    for path in sorted(DEST.rglob("*")):
        if path.is_file():
            archive.write(path, arcname=(Path(NAME) / path.relative_to(DEST)).as_posix())

zip_hash = sha256(ZIP_PATH)
ZIP_SHA.write_text(f"{zip_hash}  {ZIP_PATH.name}\n", encoding="ascii")
print(f"ITS_FREEZE_OK files={len(manifest_rows)} bytes={sum(int(r['bytes']) for r in manifest_rows)}")
print(f"FREEZE_DIR={DEST}")
print(f"ZIP={ZIP_PATH}")
print(f"ZIP_SHA256={zip_hash}")


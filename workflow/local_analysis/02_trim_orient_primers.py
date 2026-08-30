#!/usr/bin/env python3
"""Orient PacBio CCS reads and remove both full-length amplicon primers."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path


ASCII_ROOT = Path(
    r"C:\path\to\workspace"
    r"\01a01545-8f20-73e0-bf2c-0d1bc52467df"
)
PROJECT = ASCII_ROOT / "morchella_project"
INPUT_ROOT = ASCII_ROOT / "morchella_inputs"
OUT = PROJECT / "analysis"

PRIMERS = {
    "16S": ("AGRGTTYGATYMTGGCTCAG", "AAGTCGTAACAAGGTARCY"),
    "ITS": ("TACACACCGCCCGTCG", "GCATATHANTAAGSGSAGG"),
}

DATASETS = {
    "M_16S": ("16S", [f"M{year}-{rep}" for year in (2, 3, 4) for rep in (1, 2, 3)]),
    "A_16S": ("16S", [f"A{year}-{rep}" for year in (2, 3, 4) for rep in (1, 2, 3)]),
    "DH_16S": ("16S", [f"{group}-{rep}" for group in ("D", "H") for rep in (1, 2, 3)]),
    "M_ITS": ("ITS", [f"M{year}-{rep}" for year in (2, 3, 4) for rep in (1, 2, 3)]),
    "A_ITS": ("ITS", [f"A{year}-{rep}" for year in (2, 3, 4) for rep in (1, 2, 3)]),
    "DH_ITS": ("ITS", [f"{group}-{rep}" for group in ("D", "H") for rep in (1, 2, 3)]),
}


def main() -> None:
    trimmed = OUT / "trimmed"
    reports = OUT / "qc" / "cutadapt_json"
    logs = OUT / "qc" / "cutadapt_logs"
    for directory in (trimmed, reports, logs):
        directory.mkdir(parents=True, exist_ok=True)

    total = sum(len(samples) for _, samples in DATASETS.values())
    completed = 0
    summaries = []
    for dataset, (marker, samples) in DATASETS.items():
        front, back = PRIMERS[marker]
        linked_adapter = f"^{front}...{back}$"
        for sample in samples:
            completed += 1
            input_path = INPUT_ROOT / dataset / sample / f"{sample}_withprimer.fastq.gz"
            output_path = trimmed / marker / f"{sample}_{marker}_trimmed.fastq.gz"
            json_path = reports / f"{sample}_{marker}.json"
            log_path = logs / f"{sample}_{marker}.log"
            output_path.parent.mkdir(parents=True, exist_ok=True)
            if not input_path.exists():
                raise FileNotFoundError(input_path)
            command = [
                sys.executable, "-m", "cutadapt",
                "-g", linked_adapter,
                "--revcomp",
                "--discard-untrimmed",
                "--error-rate", "0.1",
                "--overlap", str(min(len(front), len(back))),
                "--cores", "1",
                "--json", str(json_path),
                "--output", str(output_path),
                str(input_path),
            ]
            print(f"[{completed:02d}/{total}] {sample} {marker}", flush=True)
            result = subprocess.run(command, capture_output=True, text=True)
            log_path.write_text(result.stdout + "\n" + result.stderr, encoding="utf-8")
            if result.returncode != 0:
                raise RuntimeError(f"cutadapt failed for {sample} {marker}; see {log_path}")
            report = json.loads(json_path.read_text(encoding="utf-8"))
            input_reads = report["read_counts"]["input"]
            output_reads = report["read_counts"]["output"]
            reverse_complemented = report["read_counts"].get("reverse_complemented", 0)
            summaries.append({
                "sample_marker_id": f"{sample}_{marker}",
                "biological_sample_id": sample,
                "marker": marker,
                "input_reads": input_reads,
                "trimmed_reads": output_reads,
                "retention_rate": output_reads / input_reads,
                "reverse_complemented_reads": reverse_complemented,
                "reverse_complemented_rate": reverse_complemented / input_reads,
                "trimmed_fastq": str(output_path),
            })

    import csv
    summary_path = OUT / "qc" / "primer_trimming_summary.tsv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)
    print(f"Completed primer orientation/trimming for {total} FASTQ files.")


if __name__ == "__main__":
    main()


#!/usr/bin/env python3
"""Read-only inventory and raw PacBio CCS FASTQ quality audit."""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
import re
import statistics
from collections import Counter
from pathlib import Path


PROJECT = Path(r"D:\path\to\morchella_microbiome")
OUT = PROJECT / "analysis"

PRIMERS = {
    "16S": ("AGRGTTYGATYMTGGCTCAG", "RGYTACCTTGTTACGACTT"),
    "ITS": ("TACACACCGCCCGTCG", "CCTSCSCTTANTDATATGC"),
}

IUPAC = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "M": {"A", "C"},
    "K": {"G", "T"}, "S": {"G", "C"}, "W": {"A", "T"},
    "H": {"A", "C", "T"}, "B": {"C", "G", "T"},
    "V": {"A", "C", "G"}, "D": {"A", "G", "T"},
    "N": {"A", "C", "G", "T"},
}
COMP = {
    "A": "T", "T": "A", "C": "G", "G": "C", "R": "Y", "Y": "R",
    "M": "K", "K": "M", "S": "S", "W": "W", "H": "D", "D": "H",
    "B": "V", "V": "B", "N": "N",
}


def reverse_complement(sequence: str) -> str:
    return "".join(COMP.get(base, "N") for base in reversed(sequence.upper()))


def primer_matches(sequence: str, primer: str) -> bool:
    if len(sequence) != len(primer):
        return False
    return all(base in IUPAC[code] for base, code in zip(sequence.upper(), primer))


def exact_orientation(sequence: str, marker: str) -> str:
    forward, reverse = PRIMERS[marker]
    reverse_rc = reverse_complement(reverse)
    forward_rc = reverse_complement(forward)
    if primer_matches(sequence[: len(forward)], forward) and primer_matches(
        sequence[-len(reverse_rc):], reverse_rc
    ):
        return "forward"
    if primer_matches(sequence[: len(reverse)], reverse) and primer_matches(
        sequence[-len(forward_rc):], forward_rc
    ):
        return "reverse"
    return "unmatched"


def percentile(values: list[int], p: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    k = (len(ordered) - 1) * p
    lo, hi = math.floor(k), math.ceil(k)
    if lo == hi:
        return float(ordered[lo])
    return ordered[lo] * (hi - k) + ordered[hi] * (k - lo)


def marker_from_path(path: Path) -> str:
    text = str(path).lower()
    return "16S" if "16s" in text else "ITS"


def sample_metadata(sample: str) -> dict[str, str]:
    prefix = sample[0].upper()
    if prefix in {"M", "A"}:
        history, replicate = sample[1:].split("-")
        return {
            "cohort": "cultivation_paired",
            "condition": "before_cultivation" if prefix == "M" else "after_cultivation",
            "continuous_cropping_years": history,
            "greenhouse_id": f"GH{history}",
            "spatial_point": replicate,
            "pair_id": f"{history}-{replicate}",
            "independent_unit": f"GH{history}",
            "subsample_role": "within_greenhouse_spatial_subsample",
        }
    disease = prefix == "D"
    replicate = sample.split("-")[1]
    return {
        "cohort": "white_mold_case_example",
        "condition": "diseased" if disease else "healthy",
        "continuous_cropping_years": "",
        "greenhouse_id": "GH_white_mold",
        "spatial_point": "disease_point" if disease else "healthy_point",
        "pair_id": "",
        "independent_unit": "disease_point" if disease else "healthy_point",
        "subsample_role": f"subsample_{replicate}_within_one_point",
    }


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_fastq(path: Path, marker: str) -> dict[str, object]:
    lengths: list[int] = []
    quality_sum = 0
    quality_bases = 0
    orientations: Counter[str] = Counter()
    reads = 0
    bases = 0
    min_q = 99
    max_q = 0
    with gzip.open(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip().upper()
            plus = handle.readline()
            quality = handle.readline().strip()
            if not plus or len(sequence) != len(quality):
                raise ValueError(f"Malformed FASTQ record in {path}")
            qvals = [ord(char) - 33 for char in quality]
            reads += 1
            bases += len(sequence)
            lengths.append(len(sequence))
            quality_sum += sum(qvals)
            quality_bases += len(qvals)
            min_q = min(min_q, min(qvals))
            max_q = max(max_q, max(qvals))
            orientations[exact_orientation(sequence, marker)] += 1
    return {
        "raw_reads": reads,
        "raw_bases": bases,
        "mean_length": round(statistics.mean(lengths), 2),
        "median_length": round(statistics.median(lengths), 2),
        "length_p01": round(percentile(lengths, 0.01), 2),
        "length_p05": round(percentile(lengths, 0.05), 2),
        "length_p95": round(percentile(lengths, 0.95), 2),
        "length_p99": round(percentile(lengths, 0.99), 2),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_phred": round(quality_sum / quality_bases, 3),
        "min_phred": min_q,
        "max_phred": max_q,
        "exact_forward_reads": orientations["forward"],
        "exact_reverse_reads": orientations["reverse"],
        "exact_both_primer_reads": orientations["forward"] + orientations["reverse"],
        "exact_both_primer_rate": round(
            (orientations["forward"] + orientations["reverse"]) / reads, 6
        ),
        "unmatched_reads": orientations["unmatched"],
    }


def expected_md5(path: Path) -> str:
    manifest = next(path.parents[1].glob("Rawdata_MD5.txt"), None)
    if manifest is None:
        for parent in path.parents:
            candidate = parent / "Rawdata_MD5.txt"
            if candidate.exists():
                manifest = candidate
                break
    if manifest is None:
        return ""
    for line in manifest.read_text(encoding="utf-8").splitlines():
        parts = line.split(maxsplit=1)
        if len(parts) == 2 and Path(parts[1]).name == path.name:
            return parts[0].lower()
    return ""


def main() -> None:
    (OUT / "metadata").mkdir(parents=True, exist_ok=True)
    (OUT / "qc").mkdir(parents=True, exist_ok=True)
    fastqs = sorted(PROJECT.glob("**/*_withprimer.fastq.gz"))
    if len(fastqs) != 48:
        raise RuntimeError(f"Expected 48 FASTQ files, found {len(fastqs)}")

    metadata_rows = []
    qc_rows = []
    for index, path in enumerate(fastqs, start=1):
        sample = path.name.replace("_withprimer.fastq.gz", "")
        marker = marker_from_path(path)
        meta = sample_metadata(sample)
        forward, reverse = PRIMERS[marker]
        actual_md5 = md5sum(path)
        expected = expected_md5(path)
        metadata_rows.append({
            "sample_marker_id": f"{sample}_{marker}",
            "biological_sample_id": sample,
            "marker": marker,
            **meta,
            "forward_primer": forward,
            "reverse_primer": reverse,
            "raw_fastq": str(path),
            "expected_md5": expected,
            "actual_md5": actual_md5,
            "md5_pass": str(bool(expected) and expected == actual_md5).upper(),
        })
        print(f"[{index:02d}/{len(fastqs)}] auditing {sample} {marker}", flush=True)
        qc_rows.append({
            "sample_marker_id": f"{sample}_{marker}",
            "biological_sample_id": sample,
            "marker": marker,
            **audit_fastq(path, marker),
        })

    metadata_fields = list(metadata_rows[0])
    with (OUT / "metadata" / "sample_metadata.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=metadata_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(metadata_rows)

    qc_fields = list(qc_rows[0])
    with (OUT / "qc" / "raw_fastq_qc.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=qc_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(qc_rows)

    if not all(row["md5_pass"] == "TRUE" for row in metadata_rows):
        raise RuntimeError("At least one FASTQ failed MD5 verification")
    print("Completed: 48 FASTQ files; all MD5 checks passed.")


if __name__ == "__main__":
    main()


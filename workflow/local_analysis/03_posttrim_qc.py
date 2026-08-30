#!/usr/bin/env python3
"""Summarize primer-trimmed PacBio CCS reads and propose length guardrails."""

from __future__ import annotations

import csv
import gzip
import math
import statistics
from collections import defaultdict
from pathlib import Path


PROJECT = Path(r"D:\path\to\morchella_microbiome")
ANALYSIS = PROJECT / "analysis"


def percentile(values: list[int], p: float) -> float:
    ordered = sorted(values)
    k = (len(ordered) - 1) * p
    lo, hi = math.floor(k), math.ceil(k)
    if lo == hi:
        return float(ordered[lo])
    return ordered[lo] * (hi - k) + ordered[hi] * (k - lo)


def fastq_stats(path: Path) -> dict[str, float | int]:
    lengths = []
    qsum = 0
    bases = 0
    with gzip.open(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip()
            handle.readline()
            quality = handle.readline().strip()
            if len(sequence) != len(quality):
                raise ValueError(f"Malformed FASTQ: {path}")
            lengths.append(len(sequence))
            qsum += sum(ord(char) - 33 for char in quality)
            bases += len(quality)
    return {
        "reads": len(lengths),
        "mean_length": round(statistics.mean(lengths), 2),
        "median_length": round(statistics.median(lengths), 2),
        "p001": round(percentile(lengths, 0.001), 2),
        "p01": round(percentile(lengths, 0.01), 2),
        "p05": round(percentile(lengths, 0.05), 2),
        "p95": round(percentile(lengths, 0.95), 2),
        "p99": round(percentile(lengths, 0.99), 2),
        "p999": round(percentile(lengths, 0.999), 2),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_phred": round(qsum / bases, 3),
    }


def main() -> None:
    trimming = {}
    with (ANALYSIS / "qc" / "primer_trimming_summary.tsv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            trimming[row["sample_marker_id"]] = row

    rows = []
    pooled_lengths: dict[str, list[int]] = defaultdict(list)
    for path in sorted((ANALYSIS / "trimmed").glob("**/*_trimmed.fastq.gz")):
        marker = path.parent.name
        sample_marker = path.name.replace("_trimmed.fastq.gz", "")
        stats = fastq_stats(path)
        input_reads = int(trimming[sample_marker]["input_reads"])
        rows.append({
            "sample_marker_id": sample_marker,
            "marker": marker,
            "input_reads": input_reads,
            **stats,
            "retention_rate": round(int(stats["reads"]) / input_reads, 6),
        })
        with gzip.open(path, "rt") as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline().strip()
                handle.readline(); handle.readline()
                pooled_lengths[marker].append(len(sequence))

    fields = list(rows[0])
    with (ANALYSIS / "qc" / "posttrim_fastq_qc.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)

    lines = [
        "# 第一阶段PacBio全长扩增子QC报告",
        "",
        "## 数据完整性",
        "",
        "- 48个原始FASTQ全部通过MD5校验。",
        "- 24个生物样品均同时具有16S和ITS数据。",
        "- 所有reads使用实际实验引物进行双端识别、反向互补定向和引物去除。",
        "- 处理参数：cutadapt linked adapter、IUPAC匹配、最大错误率10%、必须同时识别两端引物。",
        "",
        "## 去引物总体结果",
        "",
        "| Marker | 样本数 | 输入reads | 保留reads | 加权保留率 | 合并P1–P99读长 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for marker in ("16S", "ITS"):
        subset = [r for r in rows if r["marker"] == marker]
        input_n = sum(int(r["input_reads"]) for r in subset)
        kept_n = sum(int(r["reads"]) for r in subset)
        lens = pooled_lengths[marker]
        lines.append(
            f"| {marker} | {len(subset)} | {input_n:,} | {kept_n:,} | "
            f"{100*kept_n/input_n:.2f}% | {percentile(lens,0.01):.0f}–{percentile(lens,0.99):.0f} bp |"
        )
    lines.extend([
        "",
        "## ASV去噪前建议",
        "",
        "- 16S和ITS必须分开建立误差模型和ASV表。",
        "- 使用PacBio CCS单端模式，不采用Illumina PE250双端拼接流程。",
        "- DADA2阶段采用maxEE过滤，并根据本报告的全体读长分布设置宽松长度护栏。",
        "- 白霉病D/H为每组一个点位的三个子样本，只作描述性支持，不作n=3 vs n=3总体疾病推断。",
        "- M/A虽为同点配对，但采样和测序批次与种植前后重合，需保留批次局限说明。",
        "",
    ])
    (ANALYSIS / "qc" / "第一阶段_QC报告.md").write_text("\n".join(lines), encoding="utf-8")
    print("Wrote posttrim_fastq_qc.tsv and 第一阶段_QC报告.md")


if __name__ == "__main__":
    main()


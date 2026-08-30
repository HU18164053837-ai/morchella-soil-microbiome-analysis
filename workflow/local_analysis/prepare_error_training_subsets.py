#!/usr/bin/env python3
"""Create small, balanced FASTQ subsets used only to learn PacBio error rates."""

from pathlib import Path
import gzip


ROOT = Path(
    r"C:\path\to\workspace"
    r"\01a01545-8f20-73e0-bf2c-0d1bc52467df"
)
SAMPLES = {
    "M": ["M2-1", "M3-1", "M4-1"],
    "A": ["A2-1", "A3-1", "A4-1"],
    "DH": ["D-1", "H-1", "H-2"],
}
READS_PER_SAMPLE = 2500


def subset_fastq(source: Path, target: Path, n: int) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source, "rt") as src, gzip.open(target, "wt", compresslevel=6) as dst:
        for _ in range(n):
            record = [src.readline() for _ in range(4)]
            if not record[0]:
                break
            if not record[3]:
                raise ValueError(f"Incomplete FASTQ record in {source}")
            dst.writelines(record)


def main() -> None:
    for marker in ("16S", "ITS"):
        filtered = ROOT / "dada2_runtime" / marker / "filtered"
        for batch, sample_ids in SAMPLES.items():
            for sample in sample_ids:
                source = filtered / f"{sample}_{marker}_filtered.fastq.gz"
                target = (
                    ROOT / "dada2_runtime" / marker / "error_training" / batch
                    / f"{sample}_{marker}_errortrain.fastq.gz"
                )
                subset_fastq(source, target, READS_PER_SAMPLE)
                print(marker, batch, sample, target)


if __name__ == "__main__":
    main()


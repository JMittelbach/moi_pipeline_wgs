#!/usr/bin/env python3
"""Generate a deterministic, non-biological paired-end pipeline fixture."""

from __future__ import annotations

import csv
import gzip
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parent
FASTQ_DIR = ROOT / "fastq"
REFERENCE_DIR = ROOT / "reference"
CONTIG = "synthetic_chr"
REFERENCE_LENGTH = 7500
READ_LENGTH = 150
FRAGMENT_LENGTH = 250
DEPTH_PER_TARGET = 60
TARGET_POSITIONS = tuple(range(400, 400 + 24 * 280, 280))
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def wrapped_fasta(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def choose_alt(reference_base: str, index: int) -> str:
    alternatives = [base for base in "ACGT" if base != reference_base]
    return alternatives[index % len(alternatives)]


def alt_count(sample: str, target_index: int) -> int:
    if sample == "synthetic_clonal":
        return DEPTH_PER_TARGET if target_index % 2 else 0
    # Four repeatable within-sample allele fractions: 0.20, 0.35, 0.65, 0.80.
    return (12, 21, 39, 48)[target_index % 4]


def write_fastq_pair(sample: str, reference: str, targets: list[tuple[int, str, str]]) -> None:
    r1_path = FASTQ_DIR / f"{sample}_R1.fastq.gz"
    r2_path = FASTQ_DIR / f"{sample}_R2.fastq.gz"
    quality = "I" * READ_LENGTH
    with gzip.open(r1_path, "wt", encoding="ascii", newline="\n") as r1_handle, gzip.open(
        r2_path, "wt", encoding="ascii", newline="\n"
    ) as r2_handle:
        for target_index, (position, _ref, alt) in enumerate(targets):
            alternate_fragments = alt_count(sample, target_index)
            for fragment_index in range(DEPTH_PER_TARGET):
                # Unique endpoints prevent the synthetic molecules from being
                # removed as PCR duplicates. The target always lies in read 1.
                offset = 40 + fragment_index
                fragment_start = position - 1 - offset
                fragment = reference[fragment_start : fragment_start + FRAGMENT_LENGTH]
                read1 = fragment[:READ_LENGTH]
                read2 = reverse_complement(fragment[-READ_LENGTH:])
                if fragment_index < alternate_fragments:
                    read1 = read1[:offset] + alt + read1[offset + 1 :]
                name = f"{sample}:target{target_index + 1:02d}:molecule{fragment_index + 1:02d}"
                r1_handle.write(f"@{name}/1\n{read1}\n+\n{quality}\n")
                r2_handle.write(f"@{name}/2\n{read2}\n+\n{quality}\n")


def main() -> None:
    FASTQ_DIR.mkdir(parents=True, exist_ok=True)
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)

    rng = random.Random(20260828)
    reference = "".join(rng.choice("ACGT") for _ in range(REFERENCE_LENGTH))
    targets = [
        (position, reference[position - 1], choose_alt(reference[position - 1], index))
        for index, position in enumerate(TARGET_POSITIONS)
    ]

    fasta_path = REFERENCE_DIR / "synthetic_reference.fna"
    fasta_path.write_text(f">{CONTIG}\n{wrapped_fasta(reference)}\n", encoding="ascii")

    target_path = REFERENCE_DIR / "synthetic.targets.tsv"
    with target_path.open("w", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for position, ref, alt in targets:
            writer.writerow((CONTIG, position, f"{ref},{alt}"))

    panel_path = REFERENCE_DIR / "synthetic_population_panel.tsv.gz"
    with gzip.open(panel_path, "wt", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("chrom", "pos", "ref", "alt", "Global"))
        frequencies = (0.10, 0.20, 0.30, 0.40, 0.50)
        for index, (position, ref, alt) in enumerate(targets):
            writer.writerow((CONTIG, position, ref, alt, frequencies[index % len(frequencies)]))

    for sample in ("synthetic_clonal", "synthetic_mixed"):
        write_fastq_pair(sample, reference, targets)

    with (ROOT / "samples.tsv").open("w", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("sample_id", "R1", "R2"))
        for sample in ("synthetic_clonal", "synthetic_mixed"):
            writer.writerow((sample, f"{sample}_R1.fastq.gz", f"{sample}_R2.fastq.gz"))

    print(f"Created deterministic fixture in {ROOT}")
    print(f"Targets: {len(targets)}; read pairs per sample: {len(targets) * DEPTH_PER_TARGET}")


if __name__ == "__main__":
    main()

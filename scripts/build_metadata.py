#!/usr/bin/env python3
"""Build a reviewable sample/lane/read manifest from FASTQ files in data/."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

from metadata_utils import FASTQ_SUFFIXES, is_index_read, parse_fastq_name


def is_fastq(path: Path) -> bool:
    lowered = path.name.lower()
    return any(lowered.endswith(suffix) for suffix in FASTQ_SUFFIXES)


def discover_fastqs(data_dir: Path) -> list[Path]:
    files: list[Path] = []
    skip_dirs = {"merged", ".metadata_build"}
    for root, dirs, names in os.walk(data_dir):
        dirs[:] = sorted(name for name in dirs if name not in skip_dirs and not name.startswith("."))
        for name in sorted(names):
            path = Path(root) / name
            if path.is_file() and is_fastq(path):
                files.append(path)
    return sorted(files)


def build(data_dir: Path, metadata_path: Path) -> int:
    data_dir = data_dir.resolve()
    metadata_path = metadata_path.resolve()
    if not data_dir.is_dir():
        print(f"[ERROR] data directory does not exist: {data_dir}", file=sys.stderr)
        print("        Create it, copy the FASTQ/FASTQ.GZ files into it, and rerun.", file=sys.stderr)
        return 1

    fastqs = discover_fastqs(data_dir)
    if not fastqs:
        print(f"[ERROR] no FASTQ or FASTQ.GZ files found below {data_dir}", file=sys.stderr)
        print("        Expected names such as sample_S1_L001_R1_001.fastq.gz and R2.", file=sys.stderr)
        return 1

    rows: list[tuple[str, str, str, str]] = []
    ignored_index_reads: list[Path] = []
    errors: list[str] = []
    for path in fastqs:
        if is_index_read(path):
            ignored_index_reads.append(path)
            continue
        try:
            sample, lane, read = parse_fastq_name(path)
        except ValueError as exc:
            errors.append(str(exc))
            continue
        rel = path.relative_to(data_dir).as_posix()
        rows.append((rel, lane, sample, read))

    if errors:
        print("[ERROR] could not classify all FASTQ filenames:", file=sys.stderr)
        for error in errors:
            print(f"        - {error}", file=sys.stderr)
        print("        Rename files to a supported Illumina pattern or edit the manifest manually.", file=sys.stderr)
        return 1
    if not rows:
        print("[ERROR] no paired R1/R2 FASTQs found (only index reads were present)", file=sys.stderr)
        return 1

    lane_reads: dict[tuple[str, str], set[str]] = defaultdict(set)
    for _file, lane, sample, read in rows:
        lane_reads[(sample, lane)].add(read)
    incomplete = [(sample, lane, sorted(reads)) for (sample, lane), reads in sorted(lane_reads.items()) if reads != {"R1", "R2"}]
    if incomplete:
        print("[ERROR] every sample/lane must have both R1 and R2:", file=sys.stderr)
        for sample, lane, reads in incomplete:
            print(f"        - {sample} {lane}: found {','.join(reads)}", file=sys.stderr)
        print("        Add the missing mate or correct data/metadata.txt manually.", file=sys.stderr)
        return 1

    rows.sort(key=lambda row: (row[2], row[1], row[3], row[0]))
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = metadata_path.with_suffix(metadata_path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("file", "lane", "sample", "read"))
        writer.writerows(rows)
    os.replace(temporary, metadata_path)

    samples = sorted({row[2] for row in rows})
    lanes = sorted({(row[2], row[1]) for row in rows})
    print(f"Metadata written: {metadata_path}")
    print(f"Detected {len(rows)} FASTQ files, {len(samples)} samples, and {len(lanes)} sample/lane groups.")
    if ignored_index_reads:
        print(f"Ignored {len(ignored_index_reads)} Illumina index-read file(s) (I1/I2).")
    for sample in samples:
        sample_lanes = sorted(lane for current_sample, lane in lanes if current_sample == sample)
        print(f"  {sample}: {', '.join(sample_lanes)}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_dir", nargs="?", type=Path, default=Path("data"))
    parser.add_argument("metadata", nargs="?", type=Path, default=Path("data/metadata.txt"))
    args = parser.parse_args()
    return build(args.data_dir, args.metadata)


if __name__ == "__main__":
    raise SystemExit(main())

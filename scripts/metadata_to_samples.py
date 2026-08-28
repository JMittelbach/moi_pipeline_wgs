#!/usr/bin/env python3
"""Validate a reviewed metadata manifest and make the pipeline sample sheet."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

from metadata_utils import parse_fastq_name


def read_metadata(metadata_path: Path, data_dir: Path) -> dict[tuple[str, str, str], list[Path]]:
    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise ValueError(f"metadata file is empty: {metadata_path}")

    reader = csv.DictReader(lines, delimiter="\t")
    required = {"file", "lane", "sample"}
    if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
        raise ValueError("metadata header must contain: file, lane, sample (optional: read)")

    grouped: dict[tuple[str, str, str], list[Path]] = defaultdict(list)
    path_assignments: dict[Path, set[tuple[str, str, str]]] = defaultdict(set)
    for line_number, row in enumerate(reader, start=2):
        filename = (row.get("file") or "").strip()
        lane = (row.get("lane") or "").strip().upper()
        sample = (row.get("sample") or "").strip()
        read = (row.get("read") or "").strip().upper()
        if not filename and not lane and not sample and not read:
            continue
        if not filename or not lane or not sample:
            raise ValueError(f"metadata line {line_number} has an empty file, lane, or sample field")
        if read not in {"R1", "R2", "1", "2"}:
            try:
                _parsed_sample, _parsed_lane, parsed_read = parse_fastq_name(filename)
            except ValueError as exc:
                raise ValueError(f"metadata line {line_number}: {exc}") from exc
            read = parsed_read
        else:
            read = f"R{read[-1]}"
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
            raise ValueError(f"metadata line {line_number}: unsafe sample id '{sample}'")
        path = Path(filename)
        if not path.is_absolute():
            path = data_dir / path
        path = path.resolve()
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"metadata line {line_number}: FASTQ is missing or empty: {path}")
        assignment = (sample, lane, read)
        grouped[assignment].append(path)
        path_assignments[path].add(assignment)

    if not grouped:
        raise ValueError(f"metadata file contains no data rows: {metadata_path}")
    reused = {path: assignments for path, assignments in path_assignments.items() if len(assignments) > 1}
    if reused:
        details = "; ".join(f"{path}: {sorted(assignments)}" for path, assignments in sorted(reused.items()))
        raise ValueError(f"the same FASTQ is assigned to multiple sample/lane/read rows ({details})")
    return grouped


def open_fastq(path: Path, mode: str):
    return gzip.open(path, mode) if path.name.lower().endswith(".gz") else path.open(mode)


def merge_fastqs(paths: list[Path], output: Path, metadata_mtime_ns: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    newest_input_ns = max(metadata_mtime_ns, *(path.stat().st_mtime_ns for path in paths))
    if output.is_file() and output.stat().st_mtime_ns >= newest_input_ns:
        return
    temporary = output.with_suffix(output.suffix + ".tmp")
    with gzip.open(temporary, "wb", compresslevel=6) as destination:
        for path in paths:
            with open_fastq(path, "rb") as source:
                shutil.copyfileobj(source, destination, length=1024 * 1024)
    os.replace(temporary, output)


def convert(metadata_path: Path, data_dir: Path, samples_path: Path) -> int:
    metadata_path = metadata_path.resolve()
    data_dir = data_dir.resolve()
    samples_path = samples_path.resolve()
    if not metadata_path.is_file():
        print(f"[ERROR] metadata file not found: {metadata_path}", file=sys.stderr)
        print("        First run: ./run_pipeline.sh build-metadata", file=sys.stderr)
        return 1
    try:
        grouped = read_metadata(metadata_path, data_dir)
    except (OSError, ValueError) as exc:
        print(f"[ERROR] invalid metadata: {exc}", file=sys.stderr)
        print("        Correct data/metadata.txt manually, then rerun the pipeline.", file=sys.stderr)
        return 1
    metadata_mtime_ns = metadata_path.stat().st_mtime_ns

    sample_lane_reads: dict[tuple[str, str], dict[str, list[Path]]] = defaultdict(lambda: defaultdict(list))
    for (sample, lane, read), paths in grouped.items():
        sample_lane_reads[(sample, lane)][read].extend(sorted(paths))
    incomplete = [
        (sample, lane, sorted(reads))
        for (sample, lane), reads in sorted(sample_lane_reads.items())
        if set(reads) != {"R1", "R2"}
    ]
    if incomplete:
        print("[ERROR] every sample/lane in metadata must have both R1 and R2:", file=sys.stderr)
        for sample, lane, reads in incomplete:
            print(f"        - {sample} {lane}: found {','.join(reads)}", file=sys.stderr)
        return 1

    merged_dir = data_dir / "merged"
    samples: dict[str, dict[str, list[Path]]] = defaultdict(lambda: {"R1": [], "R2": []})
    for (sample, lane), reads in sorted(sample_lane_reads.items()):
        for read in ("R1", "R2"):
            samples[sample][read].extend(reads[read])

    output_rows: list[tuple[str, str, str]] = []
    for sample in sorted(samples):
        paths: dict[str, Path] = {}
        for read in ("R1", "R2"):
            source_paths = samples[sample][read]
            if len(source_paths) == 1:
                paths[read] = source_paths[0]
            else:
                merged = merged_dir / f"{sample}_{read}.fastq.gz"
                merge_fastqs(source_paths, merged, metadata_mtime_ns)
                paths[read] = merged.resolve()
        try:
            r1 = paths["R1"].relative_to(data_dir).as_posix()
            r2 = paths["R2"].relative_to(data_dir).as_posix()
        except ValueError:
            # Absolute paths outside data/ are valid manual overrides.
            r1, r2 = paths["R1"].as_posix(), paths["R2"].as_posix()
        output_rows.append((sample, r1, r2))

    samples_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = samples_path.with_suffix(samples_path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("sample_id", "R1", "R2"))
        writer.writerows(output_rows)
    os.replace(temporary, samples_path)
    print(f"Validated metadata and wrote {samples_path} for {len(output_rows)} samples.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", nargs="?", type=Path, default=Path("data/metadata.txt"))
    parser.add_argument("data_dir", nargs="?", type=Path, default=Path("data"))
    parser.add_argument("samples", nargs="?", type=Path, default=Path("data/samples.tsv"))
    args = parser.parse_args()
    return convert(args.metadata, args.data_dir, args.samples)


if __name__ == "__main__":
    raise SystemExit(main())

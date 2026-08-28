#!/usr/bin/env python3
"""Shared filename and metadata helpers for the root metadata workflow."""

from __future__ import annotations

import re
from pathlib import Path


FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
_READ_AT_END = re.compile(
    r"(?:^|[._-])(?:R)?(?P<read>[12])(?:[._-](?:[0-9]{3}|[0-9]+))?$",
    re.IGNORECASE,
)
_INDEX_READ_AT_END = re.compile(r"(?:^|[._-])I[12](?:[._-](?:[0-9]{3}|[0-9]+))?$", re.IGNORECASE)
_LANE = re.compile(r"(?:^|[._-])(?:L|lane)(?P<lane>[0-9]{1,4})(?=$|[._-])", re.IGNORECASE)
_SAMPLE_NUMBER = re.compile(r"(?:^|[._-])S[0-9]+$", re.IGNORECASE)


def strip_fastq_suffix(filename: str) -> str:
    """Return a FASTQ basename without its supported extension."""

    lowered = filename.lower()
    for suffix in FASTQ_SUFFIXES:
        if lowered.endswith(suffix):
            return filename[: -len(suffix)]
    raise ValueError(f"not a FASTQ/FASTQ.GZ filename: {filename}")


def is_index_read(path: str | Path) -> bool:
    """Return true for Illumina index reads (I1/I2), which are not paired reads."""

    try:
        stem = strip_fastq_suffix(Path(path).name)
    except ValueError:
        return False
    return _INDEX_READ_AT_END.search(stem) is not None


def parse_fastq_name(path: str | Path) -> tuple[str, str, str]:
    """Infer ``(sample, lane, read)`` from a conventional Illumina filename.

    Supported examples include ``sample_S1_L001_R1_001.fastq.gz``,
    ``sample_L002_R2.fastq.gz``, ``sample_R1.fastq.gz`` and
    ``sample.2.fq.gz``.  A missing lane is represented as ``L001``.
    """

    name = Path(path).name
    stem = strip_fastq_suffix(name)
    read_match = _READ_AT_END.search(stem)
    if read_match is None:
        raise ValueError(
            f"cannot identify R1/R2 in '{name}'; expected a suffix such as "
            "_R1_001.fastq.gz, _R2.fastq.gz, or .1.fq.gz"
        )

    read = f"R{read_match.group('read')}"
    prefix = stem[: read_match.start()].rstrip("._-")

    lane_match = _LANE.search(prefix)
    has_explicit_lane = lane_match is not None
    if lane_match is None:
        lane = "L001"
    else:
        lane = f"L{int(lane_match.group('lane')):03d}"
        prefix = (prefix[: lane_match.start()] + prefix[lane_match.end() :]).strip("._-")

    # Illumina sample sheets commonly add _S1/_S2 between the sample and an
    # explicit lane. Without a lane token, retain that suffix because it may
    # be part of the user-provided sample name.
    if has_explicit_lane:
        prefix = _SAMPLE_NUMBER.sub("", prefix)
    prefix = prefix.strip("._-")
    if not prefix:
        raise ValueError(f"cannot infer a sample name from '{name}'")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", prefix):
        raise ValueError(
            f"inferred unsafe sample name '{prefix}' from '{name}'; "
            "rename the file or correct data/metadata.txt manually"
        )
    return prefix, lane, read

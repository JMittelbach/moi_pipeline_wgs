#!/usr/bin/env python3
"""Validate the compact result table and create small dependency-free SVG plots."""

from __future__ import annotations

import argparse
import csv
from html import escape
import math
from pathlib import Path
import re
import shlex
import sys
from typing import NoReturn


MAIN_FIELDS = [
    "sample_id",
    "moi",
    "moi_status",
    "bic",
    "bic_delta",
    "bic_weight",
    "confidence",
    "fws",
]
LONG_FIELDS = [
    "sample_id",
    "metric",
    "status",
    "reason",
    "value",
    "callable_sites",
    "bins_used",
    "model_k",
    "bic",
    "bic_delta",
    "bic_weight",
    "confidence",
    "pi_hat",
    "mu_hat",
]
METRICS = {"fws_moimix_compatible", "fws_direct", "binommix_moi"}
SAFE_SAMPLE = re.compile(r"^[A-Za-z0-9_.-]+$")


def fail(message: str, fix: str) -> NoReturn:
    print(f"[ERROR] {message}\n        Fix: {fix}", file=sys.stderr)
    raise SystemExit(1)


def parse_config(path: Path) -> dict[str, str]:
    if not path.exists():
        fail(f"configuration file not found: {path}", "edit config/pipeline.env or pass its path")
    values: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        try:
            parsed = shlex.split(value.strip())
        except ValueError:
            fail(f"invalid setting in {path}: {raw_line}", "use KEY=value with quoted paths when needed")
        values[key.strip()] = parsed[0] if parsed else ""
    return values


def configured_path(value: str, base: Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else base / path


def finite_number(value: str | None, field: str, sample: str, *, allow_empty: bool = True) -> float | None:
    if value is None:
        fail(f"main table has a missing {field} field for sample {sample}", "rerun steps 07-09 and inspect the metrics files")
    if value == "":
        if allow_empty:
            return None
        fail(f"main table has an empty {field} for sample {sample}", "rerun steps 07-09 and inspect OUTPUT_DIR/logs")
    try:
        number = float(value)
    except ValueError:
        fail(f"main table has a non-numeric {field} for sample {sample}: {value!r}", "rerun steps 07-09 and inspect the metrics files")
    if not math.isfinite(number):
        fail(f"main table has a non-finite {field} for sample {sample}", "rerun BinomMix and inspect its log")
    return number


def read_table(path: Path, expected_fields: list[str]) -> list[dict[str, str]]:
    if not path.exists() or not path.is_file():
        fail(f"required table is missing: {path}", "complete step 08 before running the table check")
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != expected_fields:
                fail(
                    f"{path.name} has unexpected columns: {reader.fieldnames}",
                    f"rerun step 08 to regenerate {path.name} from the current metrics files",
                )
            rows = list(reader)
    except OSError as error:
        fail(f"cannot read {path}: {error}", "check OUTPUT_DIR permissions")
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        fail(f"{path.name} contains a row with fewer fields than its header", "rerun step 08 to regenerate the table")
    if not rows:
        fail(f"{path.name} is empty", "complete step 07 and rerun step 08")
    return rows


def validate_long_table(rows: list[dict[str, str]], main_samples: set[str]) -> None:
    seen: set[tuple[str, str]] = set()
    samples: set[str] = set()
    for row in rows:
        sample = row["sample_id"]
        metric = row["metric"]
        if not SAFE_SAMPLE.fullmatch(sample):
            fail(f"long table contains unsafe sample_id {sample!r}", "use only letters, numbers, '.', '_' or '-' in sample IDs")
        if metric not in METRICS:
            fail(f"long table contains unknown metric {metric!r} for {sample}", "rerun step 07 with the repository scripts")
        key = (sample, metric)
        if key in seen:
            fail(f"long table contains a duplicate {sample}/{metric} row", "remove duplicate metrics or rerun step 08")
        seen.add(key)
        samples.add(sample)
        if metric in {"fws_moimix_compatible", "fws_direct"}:
            value = finite_number(row["value"], "Fws value", sample)
            if value is not None and not 0.0 <= value <= 1.0:
                print(f"[WARN] Fws value for {sample} is outside the usual [0,1] range: {value:g}", file=sys.stderr)
        elif row["status"] == "estimated":
            finite_number(row["value"], "MOI", sample, allow_empty=False)
            finite_number(row["bic"], "BIC", sample, allow_empty=False)
            weight = finite_number(row["bic_weight"], "BIC weight", sample, allow_empty=False)
            if not 0.0 <= weight <= 1.0:
                fail(f"long table BIC weight for {sample} is outside [0,1]: {weight}", "rerun the BinomMix model step")
    if samples != main_samples:
        fail(
            f"main table samples differ from long table (main-only={sorted(main_samples - samples)}, long-only={sorted(samples - main_samples)})",
            "rerun step 08 so both tables are generated from the same metrics directory",
        )
    for sample in sorted(samples):
        present = {metric for sid, metric in seen if sid == sample}
        if present != METRICS:
            fail(f"long table is incomplete for {sample}: found {sorted(present)}", "rerun step 07 for every sample")


def validate_main_table(rows: list[dict[str, str]]) -> set[str]:
    samples: set[str] = set()
    for row in rows:
        sample = row["sample_id"]
        if not SAFE_SAMPLE.fullmatch(sample):
            fail(f"main table contains unsafe sample_id {sample!r}", "use only letters, numbers, '.', '_' or '-' in sample IDs")
        if sample in samples:
            fail(f"main table contains duplicate sample {sample}", "rerun step 08 and check duplicate metrics files")
        samples.add(sample)
        status = row["moi_status"]
        confidence = row["confidence"]
        if status == "estimated":
            try:
                moi = int(row["moi"])
            except ValueError:
                fail(f"main table MOI is not an integer for {sample}: {row['moi']!r}", "rerun step 07 and inspect BinomMix output")
            if moi < 1:
                fail(f"main table MOI is below 1 for {sample}: {moi}", "rerun the BinomMix model step")
            finite_number(row["bic"], "BIC", sample, allow_empty=False)
            delta = finite_number(row["bic_delta"], "BIC delta", sample)
            if delta is not None and delta < 0:
                fail(f"main table BIC delta is negative for {sample}: {delta}", "rerun the BinomMix model step")
            weight = finite_number(row["bic_weight"], "BIC weight", sample, allow_empty=False)
            if not 0.0 <= weight <= 1.0:
                fail(f"main table BIC weight is outside [0,1] for {sample}: {weight}", "rerun the BinomMix model step")
            expected = "high" if weight >= 0.90 else "moderate" if weight >= 0.60 else "low"
            if confidence != expected:
                fail(f"main table confidence disagrees with BIC weight for {sample}: {confidence!r} vs {expected!r}", "rerun step 08 from the current metrics files")
        elif status == "abstain_model_failure":
            if confidence not in {"", "not_available"}:
                fail(f"abstained sample {sample} has unexpected confidence {confidence!r}", "rerun step 08 from the current metrics files")
            if row["moi"] or row["bic"] or row["bic_weight"]:
                fail(f"abstained sample {sample} has numeric MOI model fields", "inspect the corresponding metrics file")
        else:
            fail(f"main table has unknown MOI status for {sample}: {status!r}", "rerun step 07 with the repository scripts")
        fws_value = finite_number(row["fws"], "Fws value", sample)
        if fws_value is not None and not 0.0 <= fws_value <= 1.0:
            print(f"[WARN] Fws value for {sample} is outside the usual [0,1] range: {fws_value:g}", file=sys.stderr)
    return samples


def svg_document(title: str, body: str, width: int, height: int) -> str:
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n'
        f"  <title>{escape(title)}</title>\n{body}\n</svg>\n"
    )


def write_bar_plot(path: Path, title: str, rows: list[dict[str, str]], field: str, ylabel: str, colors: dict[str, str] | None = None, fixed_max: float | None = None) -> None:
    ordered = sorted(rows, key=lambda row: row["sample_id"])
    values: list[float | None] = []
    for row in ordered:
        try:
            value = float(row[field]) if row[field] else None
        except ValueError:
            value = None
        values.append(value if value is not None and math.isfinite(value) else None)
    n = max(1, len(ordered))
    width = max(640, min(1600, 110 * n + 100))
    height = 340
    left, right, top, bottom = 60, 20, 52, 110
    plot_width = width - left - right
    plot_height = height - top - bottom
    finite = [value for value in values if value is not None]
    lower = min(0.0, min(finite)) if finite else 0.0
    maximum = fixed_max if fixed_max is not None else (max(finite) if finite else 1.0)
    maximum = max(maximum, max(finite) if finite else maximum, 1.0)
    if maximum <= lower:
        maximum = lower + 1.0
    span = maximum - lower
    zero_y = top + plot_height - ((0.0 - lower) / span) * plot_height
    bar_width = max(8.0, min(54.0, plot_width / n * 0.62))
    body = [f'  <text x="{left}" y="28" font-family="sans-serif" font-size="18" font-weight="bold">{escape(title)}</text>']
    body.append(f'  <text x="16" y="{top + plot_height / 2:.1f}" transform="rotate(-90 16 {top + plot_height / 2:.1f})" text-anchor="middle" font-family="sans-serif" font-size="12">{escape(ylabel)}</text>')
    body.append(f'  <line x1="{left}" y1="{zero_y:.1f}" x2="{width - right}" y2="{zero_y:.1f}" stroke="#333"/>')
    for tick in (lower, lower + span / 2.0, maximum):
        y = top + plot_height - ((tick - lower) / span) * plot_height
        body.append(f'  <line x1="{left}" y1="{y:.1f}" x2="{width - right}" y2="{y:.1f}" stroke="#ddd"/>')
        body.append(f'  <text x="{left - 8}" y="{y + 4:.1f}" text-anchor="end" font-family="sans-serif" font-size="11">{tick:.2g}</text>')
    if colors:
        legend_x = left
        for label, color in colors.items():
            body.append(f'  <rect x="{legend_x}" y="36" width="10" height="10" fill="{color}"/>')
            body.append(f'  <text x="{legend_x + 14}" y="45" font-family="sans-serif" font-size="10">{escape(label)}</text>')
            legend_x += 82
    for index, (row, value) in enumerate(zip(ordered, values)):
        center = left + (index + 0.5) * plot_width / n
        if value is not None:
            value_y = top + plot_height - ((value - lower) / span) * plot_height
            y = min(zero_y, value_y)
            bar_height = abs(value_y - zero_y)
            key = row.get("confidence", "")
            color = colors.get(key, "#4c78a8") if colors else "#4c78a8"
            body.append(f'  <rect x="{center - bar_width / 2:.1f}" y="{y:.1f}" width="{bar_width:.1f}" height="{bar_height:.1f}" fill="{color}"/>')
            body.append(f'  <text x="{center:.1f}" y="{max(top + 12, y - 5):.1f}" text-anchor="middle" font-family="sans-serif" font-size="11">{escape(row[field])}</text>')
        else:
            body.append(f'  <line x1="{center - bar_width / 2:.1f}" y1="{zero_y - 2:.1f}" x2="{center + bar_width / 2:.1f}" y2="{zero_y - 2:.1f}" stroke="#999" stroke-width="3"/>')
            body.append(f'  <text x="{center:.1f}" y="{zero_y - 10:.1f}" text-anchor="middle" font-family="sans-serif" font-size="10">NA</text>')
        body.append(f'  <text x="{center:.1f}" y="{height - bottom + 24}" transform="rotate(-45 {center:.1f} {height - bottom + 24})" text-anchor="end" font-family="sans-serif" font-size="11">{escape(row["sample_id"])}</text>')
    path.write_text(svg_document(title, "\n".join(body), width, height), encoding="utf-8")


def make_plots(rows: list[dict[str, str]], plot_dir: Path) -> list[Path]:
    plot_dir.mkdir(parents=True, exist_ok=True)
    colors = {"high": "#2ca02c", "moderate": "#f1a340", "low": "#d73027", "not_available": "#999999"}
    outputs = [
        plot_dir / "moi_per_sample.svg",
        plot_dir / "fws_per_sample.svg",
        plot_dir / "bic_support_per_sample.svg",
    ]
    write_bar_plot(outputs[0], "MOI per sample", rows, "moi", "MOI", colors)
    write_bar_plot(outputs[1], "Fws per sample", rows, "fws", "Fws", None, fixed_max=1.0)
    write_bar_plot(outputs[2], "BinomMix BIC support per sample", rows, "bic_weight", "relative BIC support", colors, fixed_max=1.0)
    return outputs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", nargs="?", type=Path, help="pipeline.env; defaults to config/pipeline.env")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    config = (args.config or root / "config" / "pipeline.env").resolve()
    settings = parse_config(config)
    output = configured_path(settings.get("OUTPUT_DIR", "results"), root)
    plots_value = settings.get("PLOTS_DIR", "plots") or "plots"
    plot_dir = configured_path(plots_value, output if not Path(plots_value).is_absolute() else root)
    main_path = output / "moi_per_sample.tsv"
    long_path = output / "moi_fws_summary.tsv"
    print("[9/9] checking the main table and writing small plots", flush=True)
    main_rows = read_table(main_path, MAIN_FIELDS)
    samples = validate_main_table(main_rows)
    long_rows = read_table(long_path, LONG_FIELDS)
    validate_long_table(long_rows, samples)
    outputs = make_plots(main_rows, plot_dir)
    print(f"[9/9] OK: main table validated for {len(samples)} samples")
    print(f"        main table: {main_path}")
    print(f"        plots: {plot_dir} ({', '.join(path.name for path in outputs)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

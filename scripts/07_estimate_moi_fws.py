#!/usr/bin/env python3
"""Calculate Fws and run the BinomMix MOI estimator for each sample."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
from typing import NoReturn


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.name.endswith(".gz") else path.open(encoding="utf-8")


def fail(message: str, fix: str) -> NoReturn:
    print(f"[ERROR] {message}\n        Fix: {fix}", file=sys.stderr)
    raise SystemExit(1)


def read_targets(path: Path) -> list[tuple[str, int, str, str]]:
    targets: list[tuple[str, int, str, str]] = []
    seen: set[tuple[str, int, str, str]] = set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if fields[0].lower() in {"chrom", "chr"}:
                continue
            if len(fields) == 3:
                chrom, position, allele = fields
                try:
                    ref, alt = allele.split(",", 1)
                except ValueError:
                    fail(f"{path}:{line_number}: target allele must be REF,ALT", "use three tab-separated fields: chrom, pos, REF,ALT")
            elif len(fields) >= 4:
                chrom, position, ref, alt = fields[:4]
            else:
                fail(f"{path}:{line_number}: target row has too few columns", "use chrom, pos, and REF,ALT")
            try:
                key = (chrom, int(position), ref.upper(), alt.upper())
            except ValueError:
                fail(f"{path}:{line_number}: target position is not an integer", "check the target file")
            if key in seen:
                fail(f"{path}:{line_number}: duplicate target {key}", "remove duplicate sites")
            seen.add(key)
            targets.append(key)
    if not targets:
        fail(f"target file is empty: {path}", "provide the frozen Plasmodium target set")
    return targets


def read_panel(path: Path, population: str) -> dict[tuple[str, int, str, str], float]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chrom", "pos", "ref", "alt", population}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            fail(f"population panel lacks columns {sorted(required)}: {path}", f"provide chrom, pos, ref, alt, and {population} columns")
        panel: dict[tuple[str, int, str, str], float] = {}
        for line_number, row in enumerate(reader, 2):
            try:
                key = (row["chrom"], int(row["pos"]), row["ref"].upper(), row["alt"].upper())
                value = float(row[population])
            except (KeyError, TypeError, ValueError):
                fail(f"{path}:{line_number}: invalid panel row", "check positions and Global frequencies")
            if not 0.0 <= value <= 1.0:
                fail(f"{path}:{line_number}: Global frequency outside [0,1]", "correct the panel")
            if key in panel:
                fail(f"{path}:{line_number}: duplicate panel site {key}", "remove duplicate sites")
            panel[key] = value
    return panel


def read_counts(path: Path) -> dict[tuple[str, int, str, str], tuple[int, int, int]]:
    counts: dict[tuple[str, int, str, str], tuple[int, int, int]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n\r").split("\t")
            if fields[0].lower() == "chrom":
                continue
            if len(fields) != 6:
                fail(f"{path}:{line_number}: expected six count columns", "rerun step 06 to create chrom/pos/ref/alt/AD/DP")
            chrom, position, ref, alt, ad, depth = fields
            try:
                ad_ref, ad_alt = (0, 0) if ad in {"", ".", ".,."} else (int(value) for value in ad.split(","))
                dp = 0 if depth in {"", "."} else int(depth)
                key = (chrom, int(position), ref.upper(), alt.upper())
            except (ValueError, TypeError):
                fail(f"{path}:{line_number}: invalid AD/DP or position", "inspect the bcftools log")
            if min(ad_ref, ad_alt, dp) < 0 or ad_ref + ad_alt > dp:
                fail(f"{path}:{line_number}: AD is incompatible with DP", "check mapping/counting parameters")
            if key in counts:
                fail(f"{path}:{line_number}: duplicate count site {key}", "rerun step 06 with one indexed target set")
            counts[key] = (ad_ref, ad_alt, dp)
    if not counts:
        fail(f"count table is empty: {path}", "check BAM coverage and the target/reference match")
    return counts


def fws(counts: dict[tuple[str, int, str, str], tuple[int, int, int]], panel: dict[tuple[str, int, str, str], float], targets: list[tuple[str, int, str, str]], min_depth: int, max_unmodelled: float) -> tuple[str, str, int, int]:
    sums = [[0.0, 0.0, 0] for _ in range(10)]
    population_sum = 0.0
    sample_sum = 0.0
    callable_sites = 0
    bounds = [i / 20 for i in range(1, 11)]
    for key in targets:
        if key not in counts or key not in panel:
            continue
        ref_count, alt_count, depth = counts[key]
        modelled = ref_count + alt_count
        if modelled < min_depth or (depth and (depth - modelled) / depth > max_unmodelled):
            continue
        population_frequency = panel[key]
        sample_frequency = alt_count / modelled
        hp = 2.0 * population_frequency * (1.0 - population_frequency)
        hs = 2.0 * sample_frequency * (1.0 - sample_frequency)
        index = min(bisect.bisect_left(bounds, min(population_frequency, 1.0 - population_frequency)), 9)
        sums[index][0] += hp
        sums[index][1] += hs
        sums[index][2] += 1
        population_sum += hp
        sample_sum += hs
        callable_sites += 1
    numerator = denominator = 0.0
    bins_used = 0
    for hp_sum, hs_sum, n in sums:
        if n:
            hp = hp_sum / n
            hs = hs_sum / n
            numerator += hp * hs
            denominator += hp * hp
            bins_used += hp > 0
    if denominator == 0.0:
        fws_value, status, reason = "", "abstain", "no_positive_population_heterozygosity"
    else:
        fws_value, status, reason = f"{1.0 - numerator / denominator:.12g}", "estimated", ""
    direct = "" if population_sum == 0 else f"{1.0 - sample_sum / population_sum:.12g}"
    return fws_value, direct, callable_sites, bins_used


def run_sample(args: argparse.Namespace, sample: str, counts_path: Path, output: Path) -> None:
    targets = read_targets(args.targets)
    panel = read_panel(args.panel, args.population)
    counts = read_counts(counts_path)
    missing_panel = len(set(targets) - set(panel))
    extra_panel = len(set(panel) - set(targets))
    if missing_panel or extra_panel:
        fail(f"population panel and target universe differ (missing={missing_panel}, extra={extra_panel})", "use the panel generated for exactly this TARGETS file")
    missing = len(set(targets) - set(counts))
    extra = len(set(counts) - set(targets))
    if missing or extra:
        fail(f"{sample}: count table and target universe differ (missing={missing}, extra={extra})", "use exactly the same indexed TARGETS in steps 06 and 07")
    fws_value, direct, callable_sites, bins_used = fws(counts, panel, targets, args.fws_min_depth, args.max_unmodelled)
    with tempfile.NamedTemporaryFile(prefix="moi_", suffix=".tsv", delete=False) as handle:
        moi_output = Path(handle.name)
    try:
        command = ["Rscript", str(args.binommix), str(counts_path), sample, args.k_values, str(args.coverage), str(args.niter), str(moi_output)]
        result = subprocess.run(command, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stdout, end="", file=sys.stderr)
            print(result.stderr, end="", file=sys.stderr)
            fail(f"BinomMix failed for {sample}", "install moimix/flexmix with ./setup.sh and inspect the R error above")
        with moi_output.open(newline="", encoding="utf-8") as handle:
            try:
                moi = next(csv.DictReader(handle, delimiter="\t"))
            except StopIteration:
                fail(f"BinomMix returned no result for {sample}", "inspect the R log and verify that counts pass the coverage threshold")
    finally:
        moi_output.unlink(missing_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sample_id", "metric", "status", "reason", "value", "callable_sites", "bins_used", "model_k", "bic", "bic_delta", "bic_weight", "confidence", "pi_hat", "mu_hat"]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        bic_weight = moi.get("bic_weight", "")
        try:
            weight_value = float(bic_weight)
            confidence = "high" if weight_value >= 0.90 else "moderate" if weight_value >= 0.60 else "low"
        except (TypeError, ValueError):
            confidence = "not_available"
        writer.writerow({"sample_id": sample, "metric": "fws_moimix_compatible", "status": "estimated" if fws_value else "abstain", "reason": "" if fws_value else "no_positive_population_heterozygosity", "value": fws_value, "callable_sites": callable_sites, "bins_used": bins_used})
        writer.writerow({"sample_id": sample, "metric": "fws_direct", "status": "estimated" if direct else "abstain", "reason": "" if direct else "no_positive_population_heterozygosity", "value": direct, "callable_sites": callable_sites, "bins_used": bins_used})
        writer.writerow({"sample_id": sample, "metric": "binommix_moi", "status": moi.get("status", ""), "reason": moi.get("reason", ""), "value": moi.get("model_k", ""), "callable_sites": moi.get("callable_sites", ""), "model_k": moi.get("model_k", ""), "bic": moi.get("bic", ""), "bic_delta": moi.get("bic_delta", ""), "bic_weight": bic_weight, "confidence": confidence, "pi_hat": moi.get("pi_hat", ""), "mu_hat": moi.get("mu_hat", "")})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", nargs="?", type=Path, help="pipeline.env; all other options override it")
    parser.add_argument("--counts-dir", type=Path)
    parser.add_argument("--targets", type=Path)
    parser.add_argument("--panel", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--binommix", type=Path)
    parser.add_argument("--population")
    parser.add_argument("--fws-min-depth", type=int)
    parser.add_argument("--max-unmodelled", type=float)
    parser.add_argument("--k-values")
    parser.add_argument("--coverage", type=int)
    parser.add_argument("--niter", type=int)
    args = parser.parse_args()

    pipeline_root = Path(__file__).resolve().parents[1]
    config_path = (args.config or pipeline_root / "config" / "pipeline.env").resolve()
    if not config_path.exists():
        fail(f"configuration file not found: {config_path}", "edit config/pipeline.env or pass its path")
    settings: dict[str, str] = {}
    for raw_line in config_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        try:
            settings[key.strip()] = shlex.split(value.strip())[0]
        except (IndexError, ValueError):
            fail(f"invalid setting in {config_path}: {raw_line}", "use KEY=value with quoted paths when needed")

    def setting_path(key: str, default: Path) -> Path:
        value = settings.get(key)
        if not value:
            return default
        path = Path(value).expanduser()
        return path if path.is_absolute() else pipeline_root / path

    args.counts_dir = args.counts_dir or setting_path("OUTPUT_DIR", pipeline_root / "results") / "counts"
    args.targets = args.targets or setting_path("TARGETS", pipeline_root / "__missing_targets__")
    args.panel = args.panel or setting_path("POPULATION_PANEL", pipeline_root / "__missing_panel__")
    args.output_dir = args.output_dir or setting_path("OUTPUT_DIR", pipeline_root / "results") / "metrics"
    args.binommix = args.binommix or pipeline_root / "scripts" / "07_binommix.R"
    args.population = args.population or settings.get("POPULATION", "Global")
    args.fws_min_depth = args.fws_min_depth if args.fws_min_depth is not None else int(settings.get("FWS_MIN_DEPTH", "50"))
    args.max_unmodelled = args.max_unmodelled if args.max_unmodelled is not None else float(settings.get("MAX_UNMODELLED_FRACTION", "0.02"))
    args.k_values = args.k_values or settings.get("K_VALUES", "1,2,3,4,5")
    args.coverage = args.coverage if args.coverage is not None else int(settings.get("MOI_COVERAGE_THRESHOLD", "10"))
    args.niter = args.niter if args.niter is not None else int(settings.get("MOI_NITER", "1000"))
    if not args.targets.exists() or not args.panel.exists():
        fail("targets or population panel is missing", "check TARGETS and POPULATION_PANEL in config/pipeline.env")
    count_files = sorted(args.counts_dir.glob("*.tsv"))
    if not count_files:
        fail(f"no count tables found in {args.counts_dir}", "complete step 06 first")
    for index, counts_path in enumerate(count_files, 1):
        sample = counts_path.stem
        print(f"[7/9] [{index}/{len(count_files)}] {sample}: Fws + BinomMix MOI", flush=True)
        run_sample(args, sample, counts_path, args.output_dir / f"{sample}.moi_fws.tsv")
    print(f"[7/9] OK: metrics completed for {len(count_files)} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

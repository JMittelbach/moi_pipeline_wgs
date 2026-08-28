# Synthetic smoke-test data

This directory contains deterministic, entirely synthetic data for exercising
the WGS MOI/Fws pipeline. It is not a biological simulation and its MOI/Fws
values must not be scientifically interpreted.

The fixture has two paired-end samples and 24 fixed SNP targets:

- `synthetic_clonal`: every target is entirely REF or entirely ALT;
- `synthetic_mixed`: target allele fractions repeat 0.20, 0.35, 0.65, 0.80.

Every target is covered by 60 molecules with unique alignment endpoints. Reads
are 150 bp, fragments are 250 bp, and all bases have Phred quality 40. The
separate configuration at `config/pipeline.synthetic.env` keeps test outputs
under `test_runs/` and does not touch normal `processed/` or `results/` data.

The committed generator documents the fixture. To recreate the files, run it
and then BGZF-compress and index the target table:

```bash
./test_data/generate_synthetic_data.py
bgzip --force test_data/reference/synthetic.targets.tsv
tabix --force -s 1 -b 2 -e 2 test_data/reference/synthetic.targets.tsv.gz
```

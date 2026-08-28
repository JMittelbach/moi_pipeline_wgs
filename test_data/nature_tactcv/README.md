# Nature TACT-CV WGS smoke test

This directory contains one public paired-end WGS run from the study
“Population genomics and transcriptomics of *Plasmodium falciparum* in Cambodia
and Vietnam…” (Nature Communications 2024, DOI
[10.1038/s41467-024-54915-6](https://doi.org/10.1038/s41467-024-54915-6)). The
study deposits its raw genome reads under BioProject
[PRJNA1011501](https://www.ncbi.nlm.nih.gov/bioproject/1011501).

Selected run:

- Sample: `VN001-3-029`
- SRA run: `SRR25865294`
- Layout: Illumina paired-end, 797,386 read pairs (239,215,800 bases)
- FASTQ MD5: `63147928306cb02ffc7f58317e20c94b` (R1),
  `cd47edc3730ac632957df60734abe632` (R2)

The FASTQs are intentionally ignored by Git. They were downloaded from the
public ENA mirrors and can be recreated with the commands below:

```bash
mkdir -p test_data/nature_tactcv/fastq
curl -L --fail --output test_data/nature_tactcv/fastq/VN001-3-029_R1.fastq.gz \
  'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/094/SRR25865294/SRR25865294_1.fastq.gz'
curl -L --fail --output test_data/nature_tactcv/fastq/VN001-3-029_R2.fastq.gz \
  'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/094/SRR25865294/SRR25865294_2.fastq.gz'
```

`metadata.txt` and `samples.tsv` were generated from those files. To run the
full workflow, provide the Pf3D7 reference and the matching frozen Pf4 target
and population-panel resources, then run:

```bash
./scripts/00_run_pipeline.sh config/pipeline.nature_tactcv_test.env
```

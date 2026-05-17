# ddns-minion-vp1-pipeline

[![DOI](https://zenodo.org/badge/1226769149.svg)](https://zenodo.org/badge/latestdoi/1226769149)
[![Smoke test](https://github.com/adnanhaider81/ddns-minion-vp1-pipeline/actions/workflows/smoke-test.yml/badge.svg)](https://github.com/adnanhaider81/ddns-minion-vp1-pipeline/actions/workflows/smoke-test.yml)

Standalone ONT MinION VP1 analysis pipeline for DDNS stool-culture amplicon sequencing.

The pipeline benchmarks each barcode folder against a bundled enterovirus/poliovirus VP1 reference panel, generates shortlisted reference hits, attempts consensus calling for selected candidates, and writes tabular plus HTML reports for run review.

## Portfolio quick view

This repository is intended to show a complete DDNS review path: ONT barcode folders -> VP1 screening -> shortlisted reference evidence -> consensus attempt -> HTML and tabular reporting. It is useful for supervisors, collaborators, and hiring teams who want to see how sequencing evidence is converted into an auditable public-health report without exposing restricted run data.

```mermaid
flowchart LR
  A["Barcode FASTQ folders"] --> B["Read filtering"]
  B --> C["VP1 reference screening"]
  C --> D["Candidate reference ranking"]
  D --> E["Consensus attempt"]
  E --> F["QC tables"]
  F --> G["HTML report"]
```

## Public repository checklist

| Item | Status |
| --- | --- |
| README, license, citation metadata | Present |
| Reproducible environment | `environment.yml` and `requirements.txt` |
| Tests or smoke checks | `tests/validate_resources.py` plus script syntax checks |
| Example or synthetic data | Barcode-map template and bundled non-private references |
| Documentation | `docs/` plus this README |
| Output screenshot or report example | Planned after a safe synthetic run fixture is added |
| Container recipe | `Dockerfile` |
| Zenodo DOI | [10.5281/zenodo.20257879](https://doi.org/10.5281/zenodo.20257879) |

## Contents

- `bin/ont_amplicon_vp1_by_folder.sh` - user-facing wrapper.
- `bin/vp1_pipeline_internal.sh` - internal analysis engine.
- `bin/report_tables.py` and `bin/report_html.py` - report generation.
- `references/references.all.unique.fasta` - bundled VP1 reference panel.
- `resources/BarcodedPrimers.xlsx` - bundled barcode primer workbook.
- `resources/checksums.sha256` - SHA-256 checksums for bundled resources.
- `examples/barcodes.csv` - minimal barcode-map template.

No FASTQ data, run reports, sample sheets, or generated outputs are included.

## Install

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate ddns-minion-vp1-pipeline
```

Confirm the core tools are available:

```bash
minimap2 --version
samtools --version
cutadapt --version
filtlong --version
medaka --version
```

For lightweight validation without the full conda environment:

```bash
python3 -m pip install -r requirements.txt
make validate
```

## Input Layout

The pipeline expects a MinKNOW-style `fastq_pass` directory containing barcode folders:

```text
fastq_pass/
  barcode01/
    reads_1.fastq.gz
  barcode02/
    reads_1.fastq.gz
```

The barcode map is a CSV with two required columns:

```csv
barcode,sample
barcode01,stool_culture_001
barcode02,stool_culture_002
```

Use neutral sample IDs in public outputs. Do not use patient identifiers.

## Run

```bash
bash bin/ont_amplicon_vp1_by_folder.sh \
  --barcode-map examples/barcodes.csv \
  --fastq-pass /path/to/fastq_pass \
  --out /path/to/output_dir \
  --run-name stool_culture_vp1_run
```

Common options:

```bash
--threads 8
--min-len 1000
--max-len 1300
--min-cov 50
--screen-topN 20
--rank-map-topN 10
--final-topK 3
--run-report-html /path/to/minion_report.html
```

The wrapper fixes these resources internally:

- `--refs references/references.all.unique.fasta`
- `--primers-xlsx resources/BarcodedPrimers.xlsx`
- `--preset piranha-vp1`

## Outputs

The output directory contains:

- `summary.tsv`
- `screen_per_reference.tsv`
- `benchmark_per_reference.tsv`
- `overall_reference_leaderboard.tsv`
- `top_hits_report.tsv`
- `all_best_consensus.fasta`
- `Table_Sample summary information.csv`
- `Table_Composition of samples.csv`
- `report.html`
- `results_bundle.zip`

Large intermediates are kept inside the output directory and should not be committed.

## Smoke Tests

These checks do not require MinION FASTQ data:

```bash
python3 -m pip install -r requirements.txt
bash -n bin/ont_amplicon_vp1_by_folder.sh
bash -n bin/vp1_pipeline_internal.sh
python3 -m py_compile bin/report_tables.py bin/report_html.py
python3 bin/report_tables.py --help >/dev/null
python3 bin/report_html.py --help >/dev/null
python3 tests/validate_resources.py
```

## Data Policy

This repository is intended for code, fixed reference resources, and documentation only. Keep local run data outside git:

- raw FASTQ files
- barcode maps containing identifiable sample names
- MinKNOW run reports
- generated BAM/FASTQ/FASTA/report outputs

## Version

Current release: `v1.0.3`.

## Citation

Please cite the archived Zenodo release when using this workflow:

Haider, S. A. (2026). ddns-minion-vp1-pipeline (v1.0.3). Zenodo. https://doi.org/10.5281/zenodo.20257879

The all-version Zenodo concept DOI is https://doi.org/10.5281/zenodo.20257053.

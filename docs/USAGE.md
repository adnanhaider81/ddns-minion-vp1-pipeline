# Usage Notes

## Barcode Map

The barcode map must be a CSV file with `barcode` and `sample` columns. Barcode values can be written as `1`, `01`, `BC01`, or `barcode01`; the pipeline normalizes them to `barcode01`.

Use anonymized sample names in public or shared reports.

## Reference Screening

Each sample is screened once against the bundled reference panel. The highest-supported references are shortlisted and benchmarked independently. Consensus calling is attempted for references meeting the configured thresholds.

## Reporting

The final HTML report includes:

- run metadata when a MinKNOW report is provided
- sample summary table
- composition table
- top hits
- confirmed detections
- embedded FASTA links for recovered consensus sequences

## Operational Notes

This workflow is designed for stool-culture amplicon sequencing outputs from ONT MinION runs. Interpret low-read and high-N consensus calls carefully, especially when mapped reads are close to reporting thresholds.

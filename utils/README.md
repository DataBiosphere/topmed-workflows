# Content 

## Variant-caller VCF checker workflow: Python code

The variant caller checker workflow overlap between two VCFs processed by the variant calling workflow overlap of two VCFs is based on the CCDG pipeline standard and uses the GATK Genotype Concordance checker.

`run_concordance_docker.py` runs in `topmed-workflows/variant-caller/variant-caller-wdl-checker/topmed-variantcaller-checker.wdl`. The test directory contains several tests on the methods in that Python file. Be sure to read the note in the preamble of the testing class.

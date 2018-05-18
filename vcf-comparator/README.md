# vcf-comparator

### Purpose :
A CWL workflow to determine the variant overlapping between two VCFs processed by the workflows based on the [CCDG pipeline standards](https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md)

The workflow also includes a checker to validate the VCF comparison by using a specific sensitivity/precision threshold.

### Requirements/expectations 
- Two separated VCF generated for two different workflows (or same workflow in different formats).
- Reference genome must be Hg38 with ALT contigs
- Name of the output text file
- (ONLY FOR CHECKER WORKFLOW) threshold

### Outputs 
- A TXT (tab-separated) file including true-positive, false-positive, false-negative, sensitivity, precision metrics for matching SNP and INDELS in both VCFs.

### Software version requirements :
- GATK 4.beta.3 or later (Concordance tool) (see [gotc docker](https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/))

Runtime parameters are optimized for Broad's Google Cloud Platform implementation.

# vcf-comparator

### Purpose :
A CWL workflow to determine the variant overlapping between two VCFs processed by the workflows based on the [CCDG pipeline standards](https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md). One of the VCF is considered the truth in order to estimate the overlapping and returns true-positive, false-positive, false-negative, sensitivity, precision values for both SNP and INDELS.

The workflow also includes a checker to validate the VCF comparison by using a specific sensitivity/precision threshold. A numerical threshold with the expected percentage of overlapping variants between the two VCFs (e.g. 0.95) together with the two compared VCFs is provided. The checker returns success if all four metrics (SNP sensitivity, SNP precision, INDEL sensitivity and INDEL precision) are higher that the established threshold. Otherwise, the checker fails indicating the VCFs cannot be considered equivalent according to that threshold.

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

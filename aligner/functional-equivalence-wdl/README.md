# functional-equivalence-wdl

### Purpose :
A [WDL](https://github.com/openwdl/wdl) workflow based on the [CCDG pipeline standards](https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md) for processing high-throughput sequencing data.

### Requirements/expectations 
- Human whole-genome pair-end sequencing data in unmapped BAM (uBAM) format
- One or more read groups, one per uBAM file, all belonging to a single sample (SM)
- Input uBAM files must additionally comply with the following requirements:
	- filenames all have the same suffix (we use ".unmapped.bam")
	- files must pass validation by ValidateSamFile
	- reads are provided in query-sorted order
	- all reads must have an RG tag
- Reference genome must be Hg38 with ALT contigs

### Outputs 
- A CRAM file and its index.

### Software version requirements :
- GATK 4.beta.3 or later
- Picard 2.x
- BWA MEM
- Samtools (see [gotc docker](https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/))
- Python 2.7

Cromwell version support 
 - Successfully tested on v30.1
 - Tutorial to [run WDL workflows on Google Cloud](http://cromwell.readthedocs.io/en/develop/tutorials/PipelinesApi101/)
 - NOTE: Does not work on versions < v23 due to output syntax

Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
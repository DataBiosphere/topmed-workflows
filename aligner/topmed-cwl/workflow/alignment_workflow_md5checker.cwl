cwlVersion: v1.0

doc: |
    This workflow processes high-throughput sequencing data for downstream processing

    Requirements/expectations :
    - Human whole-genome pair-end sequencing data in unmapped BAM (uBAM) format
    - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
    - Input uBAM files must additionally comply with the following requirements:
      - filenames all have the same suffix (we use ".unmapped.bam")
      - files must pass validation by ValidateSamFile
      - reads are provided in query-sorted order
      - all reads must have an RG tag
    - Reference genome must be Hg38 with ALT contigs

class: Workflow
id: alignment_pipeline_md5checker
requirements:
  - class: ScatterFeatureRequirement
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  input_bams:
    type: File[]
  indexed_reference_fasta:
    type: File
    secondaryFiles: [.64.amb, .64.ann, .64.bwt, .64.pac, .64.sa, .64.alt,
    ^.dict, .fai]
  output_basename: string
  read_name_regex: string?
  script: File
  ref_dict: File
  known_indels_sites_VCFs:
    type:
      type: array
      items: File
    secondaryFiles:
      - ^.gz.tbi
  dbSNP_vcf:
    type: File
    secondaryFiles: [^.gz.tbi]
  compression_level: int
  refcram: File

outputs:
  duplicates_marked_bam:
    type: File
    outputSource: alignment_pipeline/duplicates_marked_bam
  duplicates_marked_metrics:
    type: File
    outputSource: alignment_pipeline/duplicates_marked_metrics
  sorted_bam:
    type: File
    outputSource: alignment_pipeline/sorted_bam
  bqsr_reports:
    type: File
    outputSource: alignment_pipeline/bqsr_reports
  bqsr_bam:
    type: File
    outputSource: alignment_pipeline/bqsr_bam
  cramfile:
    type: File
    outputSource: alignment_pipeline/cramfile
  check:
    type: int
    outputSource: checker_md5/check

steps:
  alignment_pipeline:
    run: ./alignment_workflow.cwl
    in:
      input_bams: input_bams
      indexed_reference_fasta: indexed_reference_fasta
      output_basename: output_basename
      read_name_regex: read_name_regex
      script: script
      ref_dict: ref_dict
      known_indels_sites_VCFs: known_indels_sites_VCFs
      dbSNP_vcf: dbSNP_vcf
      compression_level: compression_level
    out: [duplicates_marked_bam, duplicates_marked_metrics, sorted_bam, bqsr_reports, bqsr_bam, cramfile]

  checker_md5:
      run: ./checker_md5.yaml.cwl
      in:
        refcram: refcram
        reference: indexed_reference_fasta
        cram: alignment_pipeline/cramfile
      out: [check]

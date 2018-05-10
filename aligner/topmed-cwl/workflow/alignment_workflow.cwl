cwlVersion: v1.0
class: Workflow
id: alignment_pipeline
requirements:
  - class: ScatterFeatureRequirement
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
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

outputs:
  duplicates_marked_bam:
    type: File
    outputSource: MarkDuplicates/output_markduplicates_bam
  duplicates_marked_metrics:
    type: File
    outputSource: MarkDuplicates/output_markduplicates_metrics
  sorted_bam:
    type: File
    outputSource: SortSam/output_sorted_bam
  bqsr_reports:
    type: File
    outputSource: GatherBqsrReports/output_report
  bqsr_bam:
    type: File
    outputSource: SortBam/output_sorted_bam
  cramfile:
    type: File
    outputSource: ConvertToCram/output

steps:
  GetBwaVersion:
    run: ../tasks/GetBwaVersion.yaml
    in: []
    out: [version]

  SamToFastqAndBwaMemAndMba:
    run: ../tasks/SamToFastqAndBwaMemAndMba.yaml
    scatter: [input_bam]
    in:
      input_bam: input_bams
      indexed_reference_fasta: indexed_reference_fasta
      bwa_version: GetBwaVersion/version
    out: 
      [output]

  MarkDuplicates:
    run: ../tasks/MarkDuplicates.yaml
    in:
      input_bams: SamToFastqAndBwaMemAndMba/output
      output_bam:
        source: output_basename
        valueFrom: $(self + "_markdup.bam")
      metrics_filename:
        source: output_basename
        valueFrom: $(self + "_markdup.txt")
      read_name_regex: read_name_regex
    out: [output_markduplicates_bam, output_markduplicates_metrics]

  SortSam:
    run: ../tasks/SortSam.yaml
    in:
      input_bam: MarkDuplicates/output_markduplicates_bam
      output_bam:
        source: output_basename
        valueFrom: $(self + "_markdup_sort.bam")
    out: [output_sorted_bam]

  CreateSequenceGroupingTSV:
    run: ../tasks/CreateSequenceGroupingTSV.yaml
    in:
      script: script
      ref_dict: ref_dict
    out: [groups]

  Expression_createsequencegrouping:
    run: ../tasks/Expression_createsequencegrouping.yaml
    in:
      sequence_grouping_with_unmapped_tsv: CreateSequenceGroupingTSV/groups
    out: [sequence_grouping_array]

  BaseRecalibrator:
    run: ../tasks/BaseRecalibrator.yaml
    in:
      input_bam: SortSam/output_sorted_bam
      known_indels_sites_VCFs: known_indels_sites_VCFs
      dbSNP_vcf: dbSNP_vcf
      ref_fasta: indexed_reference_fasta
      sequence_group_interval: Expression_createsequencegrouping/sequence_grouping_array
    scatter: [sequence_group_interval]
    out: [recalibration_report]

  GatherBqsrReports:
    run: ../tasks/GatherBqsrReports.yaml
    in:
      input_bqsr_reports: BaseRecalibrator/recalibration_report
      output_report_filename: 
        source: output_basename
        valueFrom: $(self + "BqsrReports.csv")
    out: [output_report]

  ApplyBQSR:
    run: ../tasks/ApplyBQSR.yaml
    in:
      compression_level: compression_level
      ref_fasta: indexed_reference_fasta
      input_bam: SortSam/output_sorted_bam
      output_bam:
        source: output_basename
        valueFrom: $(self + "_Bqsr.seg.bam")
      recalibration_report: GatherBqsrReports/output_report
      sequence_group_interval: Expression_createsequencegrouping/sequence_grouping_array
    scatter: [sequence_group_interval]
    out: [recalibrated_bam]

  GatherBamFiles:
    run: ../tasks/GatherBamFiles.yaml
    in:
      input_bam: ApplyBQSR/recalibrated_bam
      output_bam_name:
        source: output_basename
        valueFrom: $(self + "_Bqsr.bam")
      compression_level: compression_level
    out: [output_bam]

  SortBam:
    run: ../tasks/SortSam.yaml
    in:
      input_bam: GatherBamFiles/output_bam
      output_bam:
        source: output_basename
        valueFrom: $(self + "_BQSR.sort.bam")
    out: [output_sorted_bam]

  ConvertToCram:
    run: ../tasks/ConvertToCram.yaml
    in:
      reference: indexed_reference_fasta
      input_bam: SortBam/output_sorted_bam
    out: [output]

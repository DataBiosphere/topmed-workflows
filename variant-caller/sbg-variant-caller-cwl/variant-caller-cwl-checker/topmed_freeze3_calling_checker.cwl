class: Workflow
cwlVersion: v1.0
id: vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/checker-workflow/12
label: Checker Workflow
inputs:
  - id: reference_file
    type: File
  - id: reference
    type: File
  - id: pedigree_file
    type: File
  - id: num_of_jobs
    type: int?
  - id: genotype_unit
    type: int
  - id: discover_unit
    type: int
  - id: chromosomes
    type: 'string[]'
  - id: inputTruthVCFFile
    type: File
  - id: reference_genome
    type:
      type: enum
      symbols:
        - hg38
        - GRCh37
      name: reference_genome
  - id: bam_cram_file
    type: File
outputs:
  - id: TestVCFFile
    outputSource:
      - topmed_variant_calling_pipeline_cwl1/genotypes
    type: File
steps:
  - id: topmed_variant_calling_pipeline_cwl1
    in:
      - id: reference
        source:
          - reference
      - id: reference_file
        source:
          - reference_file
      - id: pedigree_file
        source:
          - pedigree_file
      - id: num_of_jobs
        source:
          - num_of_jobs
      - id: genotype_unit
        source:
          - genotype_unit
      - id: discover_unit
        source:
          - discover_unit
      - id: chromosomes
        source:
          - chromosomes
      - id: reference_genome
        source:
          - reference_genome
      - id: bam_cram_file
        source:
          - bam_cram_file
    out:
      - id: called_variant_sites
      - id: genotypes
      - id: makefile_log
      - id: vcf_output
      - id: vcf_index_output
    run: ../topmed_variant_calling_pipeline.cwl
    label: TOPMed Variant Calling Pipeline CWL1
  - id: topmed_variantcaller_checker
    in:
      - id: inputTruthVCFFile
        source:
          - inputTruthVCFFile
      - id: inputTestVCFFile
        source:
          - topmed_variant_calling_pipeline_cwl1/genotypes
    out: []
    run: topmed-variantcaller-checker.cwl
    label: topmed-variantcaller-checker
requirements:
  - class: SubworkflowFeatureRequirement
$namespaces:
  sbg: 'https://sevenbridges.com/'

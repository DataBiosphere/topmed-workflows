class: Workflow
cwlVersion: v1.0
id: vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/checker-workflow/3
label: Checker Workflow
inputs:
  - id: Expected_result
    'sbg:fileTypes': TAR.GZ
    type: File
    'sbg:x': -6
    'sbg:y': 102
  - id: reference_file
    'sbg:fileTypes': TGZ
    type: File
    'sbg:x': -465.4285583496094
    'sbg:y': -429.0484924316406
  - id: reference
    'sbg:fileTypes': FA
    type: File
    'sbg:x': -467
    'sbg:y': -315.5714416503906
  - id: pedigree_file
    'sbg:fileTypes': PED
    type: File
    'sbg:x': -466.5714416503906
    'sbg:y': -199.85714721679688
  - id: bam_cram_file
    'sbg:fileTypes': 'BAM, CRAM'
    type: 'File[]'
    'sbg:x': -476.4285583496094
    'sbg:y': 428.8571472167969
  - id: bai_crai_file
    'sbg:fileTypes': 'BAI, CRAI'
    type: 'File[]'
    'sbg:x': -475.4285888671875
    'sbg:y': 545.4285888671875
  - id: num_of_jobs
    type: int?
    'sbg:x': -469
    'sbg:y': -81
  - id: genotype_unit
    type: int
    'sbg:x': -473.28570556640625
    'sbg:y': 60.85714340209961
  - id: discover_unit
    type: int
    'sbg:x': -474.1428527832031
    'sbg:y': 314.71429443359375
  - id: reference_genome_1
    type:
      type: enum
      symbols:
        - hg38
        - GRCh37
      name: reference_genome
    'sbg:x': -462.28570556640625
    'sbg:y': -543.4285888671875
  - id: chromosomes
    type: 'string[]'
    'sbg:x': -473.2857360839844
    'sbg:y': 199.85714721679688
outputs:
  - id: Match
    outputSource:
      - validation_tool/Match
    type: string?
    'sbg:x': 381
    'sbg:y': -2
steps:
  - id: validation_tool
    in:
      - id: Expected_result
        source:
          - Expected_result
      - id: Workflow_result
        source:
          - topmed_variant_calling_pipeline_cwl1/genotypes
    out:
      - id: Match
    run: topmed-variantcaller-validation-tool.cwl
    label: Validation_tool
    'sbg:x': 220
    'sbg:y': -3
  - id: topmed_variant_calling_pipeline_cwl1
    in:
      - id: reference
        source:
          - reference
      - id: bam_cram_file
        source:
          - bam_cram_file
      - id: bai_crai_file
        source:
          - bai_crai_file
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
      - id: reference_genome_1
        source:
          - reference_genome_1
    out:
      - id: called_variant_sites
      - id: genotypes
      - id: makefile_log
    run: ../topmed_variant_calling_pipeline.cwl
    label: TOPMed Variant Calling Pipeline CWL1
    'sbg:x': -52
    'sbg:y': -79
requirements:
  - class: SubworkflowFeatureRequirement
$namespaces:
  sbg: 'https://sevenbridges.com'
'sbg:appVersion':
  - v1.0
'sbg:contributors':
  - vladimir_obucina
'sbg:createdBy': vladimir_obucina
'sbg:createdOn': 1527589051
'sbg:id': vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/checker-workflow/3
'sbg:image_url': >-
  https://igor.sbgenomics.com/ns/brood/images/vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/checker-workflow/3.png
'sbg:latestRevision': 3
'sbg:modifiedBy': vladimir_obucina
'sbg:modifiedOn': 1527589527
'sbg:project': vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline
'sbg:projectName': TOPMed Freeze 3a Variant Calling Pipeline
'sbg:publisher': sbg
'sbg:revision': 3
'sbg:revisionNotes': ''
'sbg:revisionsInfo':
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1527589051
    'sbg:revision': 0
    'sbg:revisionNotes': null
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1527589334
    'sbg:revision': 1
    'sbg:revisionNotes': 'UPDATE: First Version'
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1527589457
    'sbg:revision': 2
    'sbg:revisionNotes': ''
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1527589527
    'sbg:revision': 3
    'sbg:revisionNotes': ''
'sbg:sbgMaintained': false
'sbg:validationErrors': []

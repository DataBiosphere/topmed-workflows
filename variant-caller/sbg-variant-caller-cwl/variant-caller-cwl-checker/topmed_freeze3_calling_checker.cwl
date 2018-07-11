class: Workflow
cwlVersion: v1.0
id: vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/checker-workflow/12
label: Checker Workflow
inputs:
  - id: reference_file
    'sbg:fileTypes': TGZ
    type: File
    'sbg:x': -484.54437255859375
    'sbg:y': -419.2278137207031
  - id: reference
    'sbg:fileTypes': FA
    type: File
    'sbg:x': -480.9602966308594
    'sbg:y': -298.2476806640625
  - id: pedigree_file
    'sbg:fileTypes': PED
    type: File
    'sbg:x': -481.3563232421875
    'sbg:y': -179.4556121826172
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
    'sbg:x': -471.79205322265625
    'sbg:y': -61.45560836791992
  - id: genotype_unit
    type: int
    'sbg:x': -473.28570556640625
    'sbg:y': 60.85714340209961
  - id: discover_unit
    type: int
    'sbg:x': -474
    'sbg:y': 186.1717987060547
  - id: reference_genome_1
    type:
      type: enum
      symbols:
        - hg38
        - GRCh37
      name: reference_genome
    'sbg:x': -471.3960266113281
    'sbg:y': -547.1881103515625
  - id: chromosomes
    type: 'string[]'
    'sbg:x': -474.9367980957031
    'sbg:y': 303.5867004394531
  - id: inputTruthVCFFile
    'sbg:fileTypes': TAR.GZ
    type: File
    'sbg:x': 86.9299087524414
    'sbg:y': -197.2978973388672
outputs:
  - id: genotypes
    outputSource:
      - topmed_variant_calling_pipeline_cwl1/genotypes
    type: File
    'sbg:x': 364
    'sbg:y': 101.92414093017578
steps:
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
      - id: vcf_output
      - id: vcf_index_output
    run: ../topmed_variant_calling_pipeline.cwl
    label: TOPMed Variant Calling Pipeline CWL1
    'sbg:x': -53.936790466308594
    'sbg:y': 66.25931549072266
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
    'sbg:x': 361.55609130859375
    'sbg:y': -103.00116729736328
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
'sbg:id': vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/checker-workflow/12
'sbg:image_url': >-
  https://igor.sbgenomics.com/ns/brood/images/vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/checker-workflow/12.png
'sbg:latestRevision': 12
'sbg:modifiedBy': vladimir_obucina
'sbg:modifiedOn': 1530282177
'sbg:project': vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline
'sbg:projectName': TOPMed Freeze 3a Variant Calling Pipeline
'sbg:publisher': sbg
'sbg:revision': 12
'sbg:revisionNotes': 'UPDATE: Checker back to revision 13'
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
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1527672587
    'sbg:revision': 4
    'sbg:revisionNotes': Changed validation tool to use python script for checkinf file size
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1527690481
    'sbg:revision': 5
    'sbg:revisionNotes': >-
      UPDATE: Making log.stderr and log.stdout files, with exit message in case
      of failure.
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1527691926
    'sbg:revision': 6
    'sbg:revisionNotes': ''
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1529919626
    'sbg:revision': 7
    'sbg:revisionNotes': >-
      UPDATE: Changed Checker tool, and update CWL1 which now has different
      output format
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1529919769
    'sbg:revision': 8
    'sbg:revisionNotes': ''
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1530277807
    'sbg:revision': 9
    'sbg:revisionNotes': 'UPDATE!!!: changed checker tool'
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1530279332
    'sbg:revision': 10
    'sbg:revisionNotes': >-
      UPDATE: added output to check if everything works well. Last revision (13)
      was working!!!
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1530279360
    'sbg:revision': 11
    'sbg:revisionNotes': ''
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:modifiedOn': 1530282177
    'sbg:revision': 12
    'sbg:revisionNotes': 'UPDATE: Checker back to revision 13'
'sbg:sbgMaintained': false
'sbg:validationErrors': []

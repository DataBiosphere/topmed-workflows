class: Workflow
cwlVersion: v1.0
id: >-
  vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/topmed-variant-calling-pipeline-cwl1/11
label: TOPMed Variant Calling Pipeline CWL1
inputs:
  - id: reference
    'sbg:fileTypes': FA
    type: File
    'sbg:x': -460
    'sbg:y': -266.8571472167969
  - id: bam_cram_file
    'sbg:fileTypes': 'BAM, CRAM'
    type: 'File[]'
    'sbg:x': -457.8571472167969
    'sbg:y': -143.42857360839844
  - id: bai_crai_file
    'sbg:fileTypes': 'BAI, CRAI'
    type: 'File[]'
    'sbg:x': -456.4285888671875
    'sbg:y': -17.285715103149414
  - id: reference_file
    'sbg:fileTypes': TGZ
    type: File
    'sbg:x': -119
    'sbg:y': -629.1428833007812
  - id: pedigree_file
    'sbg:fileTypes': PED
    type: File
    'sbg:x': -119.71428680419922
    'sbg:y': -515.8571166992188
  - id: num_of_jobs
    type: int?
    'sbg:x': -121.57142639160156
    'sbg:y': -401.71429443359375
  - id: genotype_unit
    type: int
    'sbg:x': -127
    'sbg:y': -50.85714340209961
  - id: discover_unit
    type: int
    'sbg:x': -133
    'sbg:y': 73
  - id: chromosomes
    type: 'string[]'
    'sbg:x': -128.14285278320312
    'sbg:y': 198.2857208251953
  - id: reference_genome_1
    type:
      type: enum
      symbols:
        - hg38
        - GRCh37
      name: reference_genome
    'sbg:x': -459.2857360839844
    'sbg:y': -385.71429443359375
outputs:
  - id: called_variant_sites
    outputSource:
      - topmed_freeze3_calling/called_variant_sites
    type: File
    'sbg:x': 418.68548583984375
    'sbg:y': -43.3526725769043
  - id: genotypes
    outputSource:
      - topmed_freeze3_calling/genotypes
    type: File
    'sbg:x': 418.1114501953125
    'sbg:y': -197.72366333007812
  - id: makefile_log
    outputSource:
      - topmed_freeze3_calling/makefile_log
    type: File?
    'sbg:x': 423.53741455078125
    'sbg:y': -332.6839599609375
steps:
  - id: verifybamid_cwl1
    in:
      - id: bai_crai_file
        source:
          - bai_crai_file
      - id: bam_cram_file
        source:
          - bam_cram_file
      - id: reference
        source:
          - reference
      - id: reference_genome
        source:
          - reference_genome_1
    out:
      - id: output_index_file
    run: steps/verifybamid/verifybamid.cwl
    label: VerifyBamID_CWL1
    scatter:
      - bam_cram_file
    'sbg:x': -233.57144165039062
    'sbg:y': -197.14285278320312
  - id: topmed_freeze3_calling
    in:
      - id: bai_crai_files
        source:
          - bai_crai_file
      - id: bam_cram_files
        source:
          - bam_cram_file
      - id: chromosomes
        default: []
        source:
          - chromosomes
      - id: discover_unit
        source:
          - discover_unit
      - id: genotype_unit
        source:
          - genotype_unit
      - id: index_files
        source:
          - verifybamid_cwl1/output_index_file
      - id: num_of_jobs
        source:
          - num_of_jobs
      - id: pedigree_file
        source:
          - pedigree_file
      - id: reference_file
        source:
          - reference_file
      - id: reference_genome
        source:
          - reference_genome_1
    out:
      - id: called_variant_sites
      - id: genotypes
      - id: makefile_log
    run: steps/topmed_freeze3_calling/topmed_freeze3_calling.cwl
    label: Topmed_freeze3_CWL1
    'sbg:x': 157.14285278320312
    'sbg:y': -198
requirements:
  - class: ScatterFeatureRequirement
'sbg:modifiedOn': 1527500740
'sbg:latestRevision': 11
'sbg:sbgMaintained': false
'sbg:validationErrors': []
'sbg:id': >-
  vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/topmed-variant-calling-pipeline-cwl1/11
'sbg:createdOn': 1526996458
'sbg:contributors':
  - vladimir_obucina
'sbg:revisionsInfo':
  - 'sbg:revision': 0
    'sbg:modifiedOn': 1526996458
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': null
  - 'sbg:revision': 1
    'sbg:modifiedOn': 1526996882
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': 'Firste Version with CWL1 tools, scatter on VerifyBAMId is of type none'
  - 'sbg:revision': 2
    'sbg:modifiedOn': 1526997137
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': Exposed ports
  - 'sbg:revision': 3
    'sbg:modifiedOn': 1526997206
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': ''
  - 'sbg:revision': 4
    'sbg:modifiedOn': 1526997265
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': ''
  - 'sbg:revision': 5
    'sbg:modifiedOn': 1527001305
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': Separated bam and bai inputs for VerifyBAMId and Topmed_freeze3
  - 'sbg:revision': 6
    'sbg:modifiedOn': 1527007399
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': Added output to VerifyBamID
  - 'sbg:revision': 7
    'sbg:modifiedOn': 1527060490
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': 'UPDATE: Removed symlinks'
  - 'sbg:revision': 8
    'sbg:modifiedOn': 1527060551
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': ''
  - 'sbg:revision': 9
    'sbg:modifiedOn': 1527071783
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': >-
      UPDATE: changed VerifyBamID, error was lack od \n sign at the end of each
      output index file.
  - 'sbg:revision': 10
    'sbg:modifiedOn': 1527085565
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': Removed index file from output
  - 'sbg:revision': 11
    'sbg:modifiedOn': 1527500740
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': 'UPDATE: GRCh37 insted of hg19'
$namespaces:
  sbg: 'https://sevenbridges.com'
'sbg:modifiedBy': vladimir_obucina
'sbg:image_url': >-
  https://igor.sbgenomics.com/ns/brood/images/vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/topmed-variant-calling-pipeline-cwl1/11.png
'sbg:publisher': sbg
'sbg:projectName': TOPMed Freeze 3a Variant Calling Pipeline
'sbg:revisionNotes': 'UPDATE: GRCh37 insted of hg19'
'sbg:appVersion':
  - v1.0
'sbg:createdBy': vladimir_obucina
'sbg:project': vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline
'sbg:revision': 11

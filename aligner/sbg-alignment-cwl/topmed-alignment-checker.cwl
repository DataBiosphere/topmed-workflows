class: Workflow
cwlVersion: v1.0
id: topmed_alignment_checker
label: TOPMed Alignment - checker
$namespaces:
  sbg: 'https://sevenbridges.com'
inputs:
  - id: reference_genome
    'sbg:fileTypes': 'FA, FASTA'
    type: File
    label: Reference for output CRAM compressing
    'sbg:x': -400
    'sbg:y': -215
  - id: input_file
    'sbg:fileTypes': CRAM
    type: File
    label: Input CRAM file
    'sbg:x': -541
    'sbg:y': -111
  - id: bwa_index
    'sbg:fileTypes': TAR
    type: File
    label: BWA index
    'sbg:x': -366
    'sbg:y': 164.0259246826172
  - id: dbsnp
    'sbg:fileTypes': 'VCF, VCF.GZ'
    type: File?
    label: dbSNP VCF file
    'sbg:x': -486
    'sbg:y': 112
  - id: expected_md5
    type: string
    'sbg:exposed': true
  - id: decomp_ref
    'sbg:fileTypes': 'FASTA, FA'
    type: File?
    label: Reference for input CRAM decompressing
    'sbg:x': -582
    'sbg:y': 6
outputs:
  - id: aligned_out
    outputSource:
      - topmed_alignment/output
    'sbg:fileTypes': CRAM
    type: File?
    label: CRAM output
    'sbg:x': 29
    'sbg:y': 125
  - id: stdout
    outputSource:
      - alignment_validation/stdout
    type: File?
    'sbg:x': 189
    'sbg:y': -201
  - id: stderr
    outputSource:
      - alignment_validation/stderr
    type: File?
    'sbg:x': 179
    'sbg:y': 48
steps:
  - id: topmed_alignment
    in:
      - id: input_file
        source: input_file
      - id: bwa_index
        source: bwa_index
      - id: reference_genome
        source: reference_genome
      - id: dbsnp
        source: dbsnp
      - id: decomp_ref
        source: decomp_ref
    out:
      - id: output
    run: topmed-alignment.cwl
    label: TOPMed Alignment
    'sbg:x': -175
    'sbg:y': -16
  - id: alignment_validation
    in:
      - id: cram
        source: topmed_alignment/output
      - id: expected_md5
        source: expected_md5
      - id: reference
        source: reference_genome
    out:
      - id: stdout
      - id: stderr
    run: steps/alignment-validation.cwl
    label: Validation
    'sbg:x': 30
    'sbg:y': -107
hints:
  - class: 'sbg:AWSInstanceType'
    value: c4.4xlarge;ebs-gp2;64
requirements:
  - class: SubworkflowFeatureRequirement
'dct:creator':
  'foaf:mbox': 'mailto:support@sbgenomics.com'
  'foaf:name': Seven Bridges

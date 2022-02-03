class: Workflow
cwlVersion: v1.0
id: topmed_alignment
doc: >-
  A CWL wrapper of the TopMed alignment workflow described here:
  https://github.com/statgen/docker-alignment Tool Author: Hyun Min Kang
  (hmkang@umich.edu) and Adrian Tan (atks@umich.edu) Wrapper Author: Marko
  Zecevic (marko.zecevic@sbgenomics.com)
label: TOPMed Alignment
$namespaces:
  sbg: 'https://sevenbridges.com'
inputs:
  - id: input_file
    'sbg:fileTypes': CRAM
    type: File
    label: Input CRAM file
    'sbg:x': -253.4855499267578
    'sbg:y': 25.186986923217773
  - id: bwa_index
    'sbg:fileTypes': TAR
    type: File
    label: BWA Index
    'sbg:x': 114.04255676269531
    'sbg:y': 179.9574432373047
  - id: reference_genome
    'sbg:fileTypes': 'FA, FASTA'
    type: File
    label: Reference for output CRAM compressing
    'sbg:x': -106.77165222167969
    'sbg:y': -187.8012237548828
  - id: dbsnp
    'sbg:fileTypes': 'VCF, VCF.GZ'
    type: File?
    label: dbSNP VCF file
    'sbg:x': 400
    'sbg:y': 222
  - id: decomp_ref
    'sbg:fileTypes': 'FASTA, FA'
    type: File?
    label: Reference for input CRAM decompressing
    'sbg:x': -152.35092163085938
    'sbg:y': 139.45030212402344
  - id: threads
    type: int?
    label: Number of threads
  - id: ram_min
    type: int?
    label: Minimum amount of RAM - pre-align/sort/post-align
    'sbg:x': -90.64398956298828
    'sbg:y': 270.19915771484375
  - id: cores_min
    type: int?
    label: Minimum number of cores - pre-align/sort/post-align
    'sbg:x': -47.17603302001953
    'sbg:y': -313.5997314453125
  - id: cores_min_1
    type: int?
    label: Minimum number of cores - alignment
    'sbg:x': 10.978857040405273
    'sbg:y': -417.4698791503906
  - id: ram_min_1
    type: int?
    label: Minimum amount of RAM - alignment
    'sbg:x': -40.91050720214844
    'sbg:y': 391.3500061035156
outputs:
  - id: output
    outputSource:
      - topmed_post_align/output
    'sbg:fileTypes': CRAM
    type: File?
    label: Output CRAM file
    'sbg:x': 861
    'sbg:y': -72
steps:
  - id: topmed_pre_align
    in:
      - id: input_file
        source: input_file
      - id: decomp_ref
        source: decomp_ref
      - id: comp_ref
        source: reference_genome
      - id: threads
        default: 1
        source: threads
      - id: ram_min
        default: 7500
        source: ram_min
      - id: cores_min
        default: 2
        source: cores_min
    out:
      - id: fastq
      - id: list
    run: steps/pre-align.cwl
    label: Pre-align 1.0
    'sbg:x': 150.79986572265625
    'sbg:y': -54.310768127441406
  - id: topmed_align
    in:
      - id: reference
        source: bwa_index
      - id: fastq
        source: topmed_pre_align/fastq
      - id: list
        source: topmed_pre_align/list
      - id: ram_min
        default: 14000
        source: ram_min_1
      - id: cores_min
        default: 8
        source: cores_min_1
    out:
      - id: cram
    run: steps/align.cwl
    label: Align 1.0
    scatter:
      - fastq
    scatterMethod: dotproduct
    'sbg:x': 307
    'sbg:y': 63.021278381347656
  - id: samtools_sort
    in:
      - id: reference
        source: reference_genome
      - id: input_file
        source: topmed_align/cram
      - id: threads
        default: 1
        source: threads
      - id: ram_min
        default: 7500
        source: ram_min
      - id: cores_min
        default: 2
        source: cores_min
    out:
      - id: output
    run: steps/samtools-sort.cwl
    label: SAMtools Sort
    scatter:
      - input_file
    scatterMethod: dotproduct
    'sbg:x': 482.89361572265625
    'sbg:y': -1.872340440750122
  - id: topmed_post_align
    in:
      - id: reference
        source: reference_genome
      - id: dbsnp
        source: dbsnp
      - id: alignment_files
        source:
          - samtools_sort/output
      - id: threads
        default: 1
        source: threads
      - id: input_cram
        source: input_file
      - id: ram_min
        default: 7500
        source: ram_min
      - id: cores_min
        default: 2
        source: cores_min
    out:
      - id: output
    run: steps/post-align.cwl
    label: Post-align
    'sbg:x': 668.5106811523438
    'sbg:y': 8.829793930053711
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: '8'
requirements:
  - class: ScatterFeatureRequirement
'dct:creator':
  'foaf:mbox': 'mailto:support@sbgenomics.com'
  'foaf:name': Seven Bridges
'sbg:categories':
  - Alignment
'sbg:links':
  - id: 'https://github.com/statgen/docker-alignment'
    label: github
'sbg:toolAuthor': Hyun Min Kang (hmkang@umich.edu) and Adrian Tan (atks@umich.edu)
'sbg:wrapperAuthor': Marko Zecevic (marko.zecevic@sbgenomics.com)

class: Workflow
cwlVersion: v1.0
id: topmed_alignment
doc: >-
  A CWL wrapper of the TopMed alignment workflow described here:
  https://github.com/statgen/docker-alignment
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
    out:
      - id: fastq
      - id: list
    run: steps/pre-align.cwl
    label: Pre-align 1.0
    'sbg:x': 130.828125
    'sbg:y': 0
  - id: topmed_align
    in:
      - id: reference
        source: bwa_index
      - id: fastq
        source: topmed_pre_align/fastq
      - id: list
        source: topmed_pre_align/list
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
      - id: input_cram
        source: input_file
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

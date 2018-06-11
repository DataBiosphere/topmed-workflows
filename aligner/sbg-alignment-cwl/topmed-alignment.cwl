class: Workflow
cwlVersion: v1.0
id: topmed_alignment
label: TOPMed Alignment
$namespaces:
  sbg: 'https://sevenbridges.com'
inputs:
  - id: input_file
    'sbg:fileTypes': BAM
    type: File
    label: Input BAM file
    'sbg:x': 0
    'sbg:y': 60.5
  - id: bwa_index
    'sbg:fileTypes': TAR
    type: File
    label: BWA Index
    'sbg:x': 114.04255676269531
    'sbg:y': 179.9574432373047
  - id: reference_genome
    'sbg:fileTypes': 'FA, FASTA'
    type: File
    label: Reference genome sequence
    'sbg:x': 239.34042358398438
    'sbg:y': -142.0425567626953
  - id: dbsnp
    'sbg:fileTypes': 'VCF, VCF.GZ'
    type: File?
    label: dbSNP VCF file
    'sbg:x': 400
    'sbg:y': 222
outputs:
  - id: output
    outputSource:
      - topmed_post_align/output
    'sbg:fileTypes': CRAM
    type: File?
    label: Output CRAM file
    'sbg:x': 800.9786987304688
    'sbg:y': 10.361701965332031
steps:
  - id: topmed_pre_align
    in:
      - id: input_file
        source: input_file
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
    out:
      - id: output
    run: steps/post-align.cwl
    label: Post-align
    'sbg:x': 668.5106811523438
    'sbg:y': 8.829793930053711
requirements:
  - class: ScatterFeatureRequirement
'dct:creator':
  'foaf:mbox': 'mailto:support@sbgenomics.com'
  'foaf:name': Seven Bridges

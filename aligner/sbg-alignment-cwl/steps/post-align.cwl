class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://sevenbridges.com'
id: marko_zecevic/topmed-alignment/topmed-post-align/3
baseCommand:
  - samtools
  - merge
inputs:
  - id: reference
    type: File
  - id: dbsnp
    type: File
    inputBinding:
      position: 3
      prefix: '--dbsnp'
      shellQuote: false
  - id: alignment_files
    type: 'File[]'
    'sbg:fileTypes': BAM
  - id: threads
    type: int?
    label: Number of threads
outputs:
  - id: output
    type: File?
    outputBinding:
      glob: '*cram'
label: Post-align
arguments:
  - position: 1
    prefix: ''
    shellQuote: false
    valueFrom: >-
      -c merged.bam *.sorted.bam && bam-non-primary-dedup dedup_LowMem
      --allReadNames --binCustom --binQualS 0:2,3:3,4:4,5:5,6:6,7:10,13:20,23:30
      --log dedup_lowmem.metrics --recab --in merged.bam --out -.ubam
  - position: 2
    prefix: '--refFile'
    shellQuote: false
    valueFrom: |-
      ${
          return inputs.reference.path
      }
  - position: 4
    prefix: ''
    shellQuote: false
    valueFrom: '| samtools view -h -C'
  - position: 5
    prefix: '-T'
    shellQuote: false
    valueFrom: |-
      ${
          return inputs.reference.path
      }
  - position: 6
    prefix: ''
    shellQuote: false
    valueFrom: '-o output.cram'
  - position: 0
    prefix: '--threads'
    shellQuote: false
    valueFrom: |-
      ${
       if (inputs.threads) {
           return inputs.threads
       } else {
           return 1
       }
      }
  - position: 7
    prefix: '--threads'
    shellQuote: false
    valueFrom: |-
      ${
       if (inputs.threads) {
           return inputs.threads
       } else {
           return 1
       }
      }
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 7500
    coresMin: 8
  - class: DockerRequirement
    dockerPull: 'statgen/alignment:1.0.0'
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.alignment_files)
  - class: InlineJavascriptRequirement
hints:
  - class: 'sbg:AWSInstanceType'
    value: c4.4xlarge;ebs-gp2;512
'sbg:appVersion':
  - v1.0
'sbg:contributors':
  - marko_zecevic
'sbg:createdBy': marko_zecevic
'sbg:createdOn': 1525523301
'sbg:id': marko_zecevic/topmed-alignment/topmed-post-align/3
'sbg:image_url': >-
  https://igor.sbgenomics.com/ns/brood/images/marko_zecevic/topmed-alignment/topmed-post-align/3.png
'sbg:latestRevision': 3
'sbg:modifiedBy': marko_zecevic
'sbg:modifiedOn': 1526052274
'sbg:project': marko_zecevic/topmed-alignment
'sbg:projectName': TOPMed alignment
'sbg:publisher': sbg
'sbg:revision': 3
'sbg:revisionNotes': catch cram on output
'sbg:revisionsInfo':
  - 'sbg:modifiedBy': marko_zecevic
    'sbg:modifiedOn': 1525523301
    'sbg:revision': 0
    'sbg:revisionNotes': Copy of marko_zecevic/topmed-align/post-align/2
  - 'sbg:modifiedBy': marko_zecevic
    'sbg:modifiedOn': 1525968584
    'sbg:revision': 1
    'sbg:revisionNotes': single output
  - 'sbg:modifiedBy': marko_zecevic
    'sbg:modifiedOn': 1526043700
    'sbg:revision': 2
    'sbg:revisionNotes': shellquote on arg3 off
  - 'sbg:modifiedBy': marko_zecevic
    'sbg:modifiedOn': 1526052274
    'sbg:revision': 3
    'sbg:revisionNotes': catch cram on output
'sbg:sbgMaintained': false
'sbg:validationErrors': []

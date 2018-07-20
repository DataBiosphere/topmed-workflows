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
    listing: |-
        ${ 
            var out = []
            out.push(inputs.reference)
            for (var i = 0; i < inputs.alignment_files.length; i++) { 
                out.push(inputs.alignment_files[i]);
            }
            return out
            
        }
  - class: InlineJavascriptRequirement
hints:
  - class: 'sbg:AWSInstanceType'
    value: c4.4xlarge;ebs-gp2;512
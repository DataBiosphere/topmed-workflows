class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://sevenbridges.com'
id: marko_zecevic_topmed_alignment_topmed_post_align_3
baseCommand:
  - samtools
  - merge
inputs:
  - id: reference
    type: File
  - id: dbsnp
    type: File?
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
  - id: input_cram
    type: File
    label: Input CRAM
    doc: Input CRAM is passed to Post-align for output naming.
    'sbg:fileTypes': CRAM
outputs:
  - id: output
    type: File?
    outputBinding:
      glob: |-
        ${
          var outname = inputs.input_cram.path.split('/').slice(-1)[0]
          outname = outname.split('.').slice(0, -1).join('.')
          outname = outname.concat(".recab.cram")
          return outname
        }
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
    prefix: '-o'
    shellQuote: false
    valueFrom: |-
      ${
        var outname = inputs.input_cram.path.split('/').slice(-1)[0]
        outname = outname.split('.').slice(0, -1).join('.')
        outname = outname.concat(".recab.cram")
        return outname
      }
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
  - position: 8
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
    coresMin: 2
  - class: DockerRequirement
    dockerPull: 'statgen/alignment:1.0.0'
  - class: InitialWorkDirRequirement
    listing: |-
        ${ 
            var out = []
            out.push(inputs.reference)
            out.push(inputs.input_cram)
            for (var i = 0; i < inputs.alignment_files.length; i++) { 
                out.push(inputs.alignment_files[i]);
            }
            return out
            
        }
  - class: InlineJavascriptRequirement
hints:
  - class: 'sbg:AWSInstanceType'
    value: m5.large;ebs-gp2;700

class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://sevenbridges.com'
id: marko_zecevic_validation_app_topmed_alignment_topmed_pre_align_7
baseCommand:
  - samtools
  - view
inputs:
  - id: input_file
    type: File
    inputBinding:
      position: 2
      shellQuote: false
    label: Input CRAM file
    'sbg:fileTypes': CRAM
  - 'sbg:category': Input file
    id: decomp_ref
    type: File?
    label: Reference for input CRAM decompressing
    'sbg:fileTypes': 'FASTA, FA'
  - 'sbg:category': Input file
    id: comp_ref
    type: File
    label: Reference for output CRAM compressing
    'sbg:fileTypes': 'FASTA, FA'
  - 'sbg:toolDefaultValue': '1'
    id: threads
    type: int?
    label: Number of threads
outputs:
  - id: fastq
    type: 'File[]?'
    outputBinding:
      glob: '*.gz'
      outputEval: '$(inheritMetadata(self, inputs.input_file))'
  - id: list
    type: File?
    outputBinding:
      glob: '*.list'
      outputEval: '$(inheritMetadata(self, inputs.input_file))'
label: Pre-align 1.0
arguments:
  - position: 3
    prefix: >-
      | bam-ext-mem-sort-manager squeeze --in -.ubam --keepDups --rmTags
      AS:i,BD:Z,BI:Z,XS:i,MC:Z,MD:Z,NM:i,MQ:i --out -.ubam | samtools sort -l 1
      -@
    shellQuote: false
    valueFrom: |-
      ${
       if (inputs.threads) {
           return inputs.threads
       } else {
           return 1
       }
      }
  - position: 4
    prefix: '-n -T'
    shellQuote: false
    valueFrom: |-
      ${
          return inputs.input_file.nameroot + '.samtools_sort_tmp'
      }
  - position: 5
    prefix: ''
    shellQuote: false
    valueFrom: '- | samtools fixmate - - | bam-ext-mem-sort-manager bam2fastq --in -.bam'
  - position: 6
    prefix: '--outBase'
    shellQuote: false
    valueFrom: |-
      ${
          return inputs.input_file.nameroot
      }
  - position: 7
    prefix: ''
    shellQuote: false
    valueFrom: '--maxRecordLimitPerFq 20000000 --sortByReadNameOnTheFly --readname --gzip'
  - position: 1
    prefix: '-uh -F 0x900 --threads'
    shellQuote: false
    valueFrom: |-
      ${
       if (inputs.threads) {
           return inputs.threads
       } else {
           return 1
       }
      }
  - position: 2
    prefix: '-T'
    valueFrom: |-
      ${
       if (inputs.decomp_ref) {
           return inputs.decomp_ref.path
       } else {
           return inputs.comp_ref.path
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
      - $(inputs.comp_ref)
  - class: InlineJavascriptRequirement
    expressionLib:
      - |-

        var setMetadata = function(file, metadata) {
            if (!('metadata' in file))
                file['metadata'] = metadata;
            else {
                for (var key in metadata) {
                    file['metadata'][key] = metadata[key];
                }
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };
hints:
  - class: 'sbg:AWSInstanceType'
    value: c4.4xlarge;ebs-gp2;512

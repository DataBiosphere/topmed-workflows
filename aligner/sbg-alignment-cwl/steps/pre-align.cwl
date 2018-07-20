class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://sevenbridges.com'
id: marko_zecevic/validation-app-topmed-alignment/topmed-pre-align/7
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
  - id: output_name
    type: string?
    label: Output name
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
       if (inputs.out_name) {
           return inputs.out_name + '.samtools_sort_tmp'  
       } else {
           return inputs.input_file.nameroot + '.samtools_sort_tmp'
       }
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
       if (inputs.out_name) {
           return inputs.out_name
       } else {
           return inputs.input_file.nameroot
       }
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
'sbg:appVersion':
  - v1.0
'sbg:contributors':
  - marko_zecevic
'sbg:createdBy': marko_zecevic
'sbg:createdOn': 1527680399
'sbg:id': marko_zecevic/validation-app-topmed-alignment/topmed-pre-align/7
'sbg:image_url': >-
  https://igor.sbgenomics.com/ns/brood/images/marko_zecevic/validation-app-topmed-alignment/topmed-pre-align/7.png
'sbg:latestRevision': 7
'sbg:modifiedBy': marko_zecevic
'sbg:modifiedOn': 1530655047
'sbg:project': marko_zecevic/validation-app-topmed-alignment
'sbg:projectName': Validation app - TOPMed alignment
'sbg:publisher': sbg
'sbg:revision': 7
'sbg:revisionNotes': ''
'sbg:revisionsInfo':
  - 'sbg:modifiedBy': marko_zecevic
    'sbg:modifiedOn': 1527680399
    'sbg:revision': 0
    'sbg:revisionNotes': Copy of marko_zecevic/topmed-alignment/topmed-pre-align/0
  - 'sbg:modifiedBy': marko_zecevic
    'sbg:modifiedOn': 1527680399
    'sbg:revision': 1
    'sbg:revisionNotes': Copy of marko_zecevic/topmed-alignment/topmed-pre-align/1
  - 'sbg:modifiedBy': marko_zecevic
    'sbg:modifiedOn': 1530267373
    'sbg:revision': 2
    'sbg:revisionNotes': threads
  - 'sbg:modifiedBy': marko_zecevic
    'sbg:modifiedOn': 1530277627
    'sbg:revision': 3
    'sbg:revisionNotes': back to default instance
  - 'sbg:modifiedBy': marko_zecevic
    'sbg:modifiedOn': 1530653426
    'sbg:revision': 4
    'sbg:revisionNotes': back to default instance
  - 'sbg:modifiedBy': marko_zecevic
    'sbg:modifiedOn': 1530654447
    'sbg:revision': 5
    'sbg:revisionNotes': larger instance
  - 'sbg:modifiedBy': marko_zecevic
    'sbg:modifiedOn': 1530654713
    'sbg:revision': 6
    'sbg:revisionNotes': first is singlethreaded
  - 'sbg:modifiedBy': marko_zecevic
    'sbg:modifiedOn': 1530655047
    'sbg:revision': 7
    'sbg:revisionNotes': ''
'sbg:sbgMaintained': false
'sbg:validationErrors': []

class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://sevenbridges.com'
id: marko_zecevic/topmed-alignment/samtools-sort/0
baseCommand:
  - samtools
  - sort
inputs:
  - id: reference
    type: File
    inputBinding:
      position: 0
      prefix: '--reference'
      shellQuote: false
  - id: input_file
    type: File
    inputBinding:
      position: 4
      shellQuote: false
outputs:
  - id: output
    type: File?
    outputBinding:
      glob: '*.sorted.bam'
      outputEval: '$(inheritMetadata(self, inputs.input_file))'
label: SAMtools Sort
arguments:
  - position: 1
    prefix: '--threads'
    shellQuote: false
    valueFrom: '1'
  - position: 3
    prefix: '-o'
    shellQuote: false
    valueFrom: |-
      ${
          input_filename = inputs.input_file.path.split('/').pop()
          return input_filename.slice(0,input_filename.lastIndexOf('.')) + '.sorted.bam'
      }
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 3000
    coresMin: 2
  - class: DockerRequirement
    dockerPull: images.sbgenomics.com/marko_zecevic/topmed_alignment
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
'sbg:appVersion':
  - v1.0
'sbg:contributors':
  - marko_zecevic
'sbg:copyOf': marko_zecevic/topmed-align/samtools-sort/1
'sbg:createdBy': marko_zecevic
'sbg:createdOn': 1525523240
'sbg:id': marko_zecevic/topmed-alignment/samtools-sort/0
'sbg:image_url': >-
  https://igor.sbgenomics.com/ns/brood/images/marko_zecevic/topmed-alignment/samtools-sort/0.png
'sbg:latestRevision': 0
'sbg:modifiedBy': marko_zecevic
'sbg:modifiedOn': 1525523240
'sbg:project': marko_zecevic/topmed-alignment
'sbg:projectName': TOPMed alignment
'sbg:publisher': sbg
'sbg:revision': 0
'sbg:revisionNotes': Copy of marko_zecevic/topmed-align/samtools-sort/1
'sbg:revisionsInfo':
  - 'sbg:modifiedBy': marko_zecevic
    'sbg:modifiedOn': 1525523240
    'sbg:revision': 0
    'sbg:revisionNotes': Copy of marko_zecevic/topmed-align/samtools-sort/1
'sbg:sbgMaintained': false
'sbg:validationErrors': []

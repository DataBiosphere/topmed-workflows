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
          var input_filename = inputs.input_file.path.split('/').pop()
          return input_filename.slice(0,input_filename.lastIndexOf('.')) + '.sorted.bam'
      }
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 7000
    coresMin: 2
  - class: DockerRequirement
    dockerPull: 'statgen/alignment:1.0.0'
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.reference)
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

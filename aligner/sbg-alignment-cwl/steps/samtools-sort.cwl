class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://sevenbridges.com'
id: marko_zecevic_topmed_alignment_samtools_sort_0
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
  - default: 1
    id: threads
    type: int?
    label: Number of threads
  - default: 7000
    'sbg:toolDefaultValue': '7500'
    id: ram_min
    type: int?
    label: Minimum amount of RAM
  - default: 2
    'sbg:toolDefaultValue': '2'
    id: cores_min
    type: int?
    label: Minimum number of cores
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
    valueFrom: $(inputs.threads)
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
    ramMin: $(inputs.ram_min)
    coresMin: $(inputs.cores_min)
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

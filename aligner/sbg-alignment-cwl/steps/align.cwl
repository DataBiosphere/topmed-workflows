class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://sevenbridges.com'
id: marko_zecevic/topmed-alignment/topmed-align/0
baseCommand:
  - chmod
inputs:
  - 'sbg:category': Input files
    format: 'FASTA, FA'
    id: reference
    type: File
  - format: 'FASTQ, FQ, FASTQ.GZ, FQ.GZ'
    id: fastq
    type: File
    inputBinding:
      position: 8
      shellQuote: false
  - 'sbg:category': Input files
    format: LIST
    id: list
    type: File
    inputBinding:
      position: 9
      shellQuote: false
outputs:
  - id: cram
    type: File?
    outputBinding:
      glob: '*.cram'
      outputEval: |-
        ${
            return inheritMetadata(self, inputs.fastq)

        }
    format: CRAM
label: Align 1.0
arguments:
  - position: 1
    separate: false
    shellQuote: false
    valueFrom: +x
  - position: 2
    separate: false
    shellQuote: false
    valueFrom: align.sh
  - position: 3
    separate: false
    shellQuote: false
    valueFrom: '&&'
  - position: 4
    separate: false
    shellQuote: false
    valueFrom: tar
  - position: 5
    separate: false
    shellQuote: false
    valueFrom: '-xf'
  - position: 6
    shellQuote: false
    valueFrom: '&& ./align.sh'
  - position: 5
    shellQuote: false
    valueFrom: |-
      ${
          return inputs.reference.path
      }
  - position: 7
    shellQuote: false
    valueFrom: |-
      ${
          reference_file = inputs.reference.path.split('/')[inputs.reference.path.split('/').length - 1]
          name = reference_file.slice(0, -4) // cut .tar extension     
          return name
      }
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 12000
    coresMin: 2
  - class: DockerRequirement
    dockerPull: images.sbgenomics.com/marko_zecevic/topmed_alignment
  - class: InitialWorkDirRequirement
    listing:
      - entryname: align.sh
        entry: |-
          ${
           s = '#!/bin/bash\n'

           s += 'ref_path=$1\n'
           s += 'input_path=$2\n'
           s += 'list=$3\n'

           s += 'line=$(grep $(basename $input_path) < $list)\n'

           s += 'line_rg=$(echo $line | cut -d \' \' -f 4- | sed -e \"s\/ \/\\\t\/g\")\n'
           s += 'input_filename=$(basename $input_path)\n'
           s += 'output_filename=$(basename $input_filename \".fastq.gz\").cram\n'

           s += 'paired_flag=\"\"\n'
           s += 'if [[ $input_file_name =~ interleaved\.fastq\.gz$ ]]\n'
           s += 'then\n'
           s += '\tpaired_flag=\"-p\"\n'
           s += 'fi\n'

           s += 'bwa mem -t 32 -K 100000000 -Y ${paired_flag} -R \"$line_rg\" $ref_path $input_path | samblaster -a --addMateTags | samtools view -@ 32 -T $ref_path -C -o $output_filename -'
          return s
          }
  - class: InlineJavascriptRequirement
    expressionLib:
      - |-
        var updateMetadata = function(file, key, value) {
            file['metadata'][key] = value;
            return file;
        };


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

        var toArray = function(file) {
            return [].concat(file);
        };

        var groupBy = function(files, key) {
            var groupedFiles = [];
            var tempDict = {};
            for (var i = 0; i < files.length; i++) {
                var value = files[i]['metadata'][key];
                if (value in tempDict)
                    tempDict[value].push(files[i]);
                else tempDict[value] = [files[i]];
            }
            for (var key in tempDict) {
                groupedFiles.push(tempDict[key]);
            }
            return groupedFiles;
        };

        var orderBy = function(files, key, order) {
            var compareFunction = function(a, b) {
                if (a['metadata'][key].constructor === Number) {
                    return a['metadata'][key] - b['metadata'][key];
                } else {
                    var nameA = a['metadata'][key].toUpperCase();
                    var nameB = b['metadata'][key].toUpperCase();
                    if (nameA < nameB) {
                        return -1;
                    }
                    if (nameA > nameB) {
                        return 1;
                    }
                    return 0;
                }
            };

            files = files.sort(compareFunction);
            if (order == undefined || order == "asc")
                return files;
            else
                return files.reverse();
        };
hints:
  - class: 'sbg:AWSInstanceType'
    value: c4.2xlarge;ebs-gp2;64
'sbg:latestRevision': 0
'sbg:revisionsInfo':
  - 'sbg:revisionNotes': Copy of marko_zecevic/topmed-align/aligner-converted-to-cwl1/2
    'sbg:modifiedBy': marko_zecevic
    'sbg:revision': 0
    'sbg:modifiedOn': 1525523285
'sbg:publisher': sbg
'sbg:modifiedOn': 1525523285
'sbg:id': marko_zecevic/topmed-alignment/topmed-align/0
'sbg:validationErrors': []
'sbg:createdBy': marko_zecevic
'sbg:image_url': >-
  https://igor.sbgenomics.com/ns/brood/images/marko_zecevic/topmed-alignment/topmed-align/0.png
'sbg:copyOf': marko_zecevic/topmed-align/aligner-converted-to-cwl1/2
'sbg:revision': 0
'sbg:createdOn': 1525523285
'sbg:modifiedBy': marko_zecevic
'sbg:projectName': TOPMed alignment
'sbg:project': marko_zecevic/topmed-alignment
'sbg:appVersion':
  - v1.0
'sbg:revisionNotes': Copy of marko_zecevic/topmed-align/aligner-converted-to-cwl1/2
'sbg:cmdPreview': >-
  chmod +x align.sh && tar -xf  /path/to/reference.fasta.tar  &&
  ./align.sh  reference.fasta  /path/to/fastq.ext  /path/to/list.ext
'sbg:contributors':
  - marko_zecevic
'sbg:sbgMaintained': false
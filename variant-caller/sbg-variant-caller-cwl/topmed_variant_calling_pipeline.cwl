class: Workflow
cwlVersion: v1.0
id: >-
  vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/topmed-variant-calling-pipeline-cwl1/17
label: TOPMed Variant Calling Pipeline CWL1
inputs:
  - id: reference
    type: File
    secondaryFiles:
      - .fai
  - id: reference_file
    type: File
  - id: pedigree_file
    type: File?
  - id: num_of_jobs
    type: int?
  - id: genotype_unit
    type: int
  - id: discover_unit
    type: int
  - id: chromosomes
    type: 'string[]'
  - id: reference_genome
    type:
      type: enum
      symbols:
        - hg38
        - GRCh37
      name: reference_genome
  - id: bam_cram_file
    type: File
    secondaryFiles:
      - |-
        ${
            return (self.basename + self.nameext.replace('m','i'))
        }
outputs:
  - id: called_variant_sites
    outputSource:
      - topmed_freeze3_calling/called_variant_sites
    type: File
  - id: genotypes
    outputSource:
      - topmed_freeze3_calling/genotypes
    type: File
  - id: makefile_log
    outputSource:
      - topmed_freeze3_calling/makefile_log
    type: File?
  - id: vcf_output
    outputSource:
      - topmed_freeze3_calling/vcf_output
    type: 'File[]?'
  - id: vcf_index_output
    outputSource:
      - topmed_freeze3_calling/vcf_index_output
    type: 'File[]?'
steps:
  - id: verifybamid_cwl1
    in:
      - id: bam_cram_file
        source:
          - bam_cram_file
      - id: reference
        source:
          - reference
      - id: reference_genome
        source:
          - reference_genome
    out:
      - id: output_index_file
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: >-
        vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/verifybamid_cwl1/10
      baseCommand: []
      inputs:
        - format: 'BAM,CRAM'
          id: bam_cram_file
          type: File
          inputBinding:
            position: 1
            valueFrom: |-
              ${
                  return ''
              }
          label: BAM/CRAM Files
          doc: Bam or Cram file for the sample
          secondaryFiles:
            - |-
              ${
                  return (self.basename + self.nameext.replace('m','i'))
              }
        - format: FA
          id: reference
          type: File
          inputBinding:
            position: 1
            shellQuote: false
            valueFrom: |-
              ${

                  return ""

              }
          label: Reference
          doc: Reference file
          secondaryFiles:
            - .fai
        - id: reference_genome
          type:
            type: enum
            symbols:
              - hg38
              - GRCh37
            name: reference_genome
          label: Reference genome
      outputs:
        - id: output_index_file
          doc: Output for topmed freeze3 pipeline
          label: Output Index File
          type: File?
          outputBinding:
            glob: |-
              ${
                  return inputs.bam_cram_file.path.split("/").pop().split(".").shift() + ".index"

              }
          format: INDEX
      label: VerifyBamID_CWL1
      arguments:
        - position: 0
          prefix: ''
          separate: false
          shellQuote: false
          valueFrom: |-
            ${
                if (inputs.reference_genome == 'GRCh37') {
                    var UDPath = "/VerifyBamID/resource/1000g.phase3.100k.b37.vcf.gz.dat.UD"
                    var BedPath = "/VerifyBamID/resource/1000g.phase3.100k.b37.vcf.gz.dat.bed"
                    var MeanPath = "/VerifyBamID/resource/1000g.phase3.100k.b37.vcf.gz.dat.mu"
                } else if (inputs.reference_genome == 'hg38') {
                    var UDPath = "/VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat.UD"
                    var BedPath = "/VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat.bed"
                    var MeanPath = "/VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat.mu"
                }

                var comm = "export PATH=$PATH:/VerifyBamID/bin/ && VerifyBamID \
             --UDPath " + UDPath + " \
             --BedPath " + BedPath + " \
             --MeanPath " + MeanPath + " \
             --Reference " + inputs.reference.path + " --BamFile " + inputs.bam_cram_file.basename
                comm += " && python make_index.py --file "
                comm += inputs.bam_cram_file.path.split("/").pop().split(".").shift()
                comm += " --path " + inputs.bam_cram_file.path + " --result result.out"

                return comm
            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 1000
          coresMin: 1
        - class: DockerRequirement
          dockerPull: 'images.sbgenomics.com/vladimir_obucina/topmed:VerifyBamID'
        - class: InitialWorkDirRequirement
          listing:
            - entryname: make_index.py
              entry: >-
                import sys

                import glob, os

                import argparse


                parser = argparse.ArgumentParser()



                parser.add_argument("--file", help="File name.", type=str)

                parser.add_argument("--path", help="Path to file.", type=str)

                parser.add_argument("--result", help="Result from ", type=str)


                args = parser.parse_args()


                with open(args.result, 'rt') as in_file:
                    contents = in_file.read()
                    broj = contents.split("Alpha:",4)
                    contamination = broj[1].rstrip()



                output_file = args.file + ".index"

                f = open (output_file, "w+")


                f.write(args.file + '\t' + args.path + '\t' + contamination +
                '\n')


                f.close()
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
      $namespaces:
        sbg: 'https://sevenbridges.com'
      
cwlVersion: v1.0

class: CommandLineTool

id: checkcram_md5

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135
  - class: InlineJavascriptRequirement

baseCommand: []
arguments:
  - position: 0
    shellQuote: false
    valueFrom: |
      samtools view $(inputs.refcram.path) -T $(inputs.reference.path)| md5sum > sum_file.txt &&
      samtools view $(inputs.cram.path) -T $(inputs.reference.path)| md5sum --check sum_file.txt

inputs:
  cram:
    type: File
  refcram:
    type: File
  reference:
    type: File
    secondaryFiles: [.fai] 

stdout: check

outputs:
  check:
    type: int
    outputBinding:
      glob: check
      loadContents: true
      outputEval: |
        ${
          var message = "The output cram file is not identical"
          var flag = 1
          if (self[0].contents.trim() == "-: OK"){
            message = "The output cram file is identical" 
            flag = 0 
            }
          return flag
         }


class: CommandLineTool
cwlVersion: v1.0
id: vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/validation_tool/15
baseCommand: []
inputs:
  - 'sbg:category': Input Files
    id: Expected_result
    type: File
    label: Expected result
    'sbg:fileTypes': TAR.GZ
  - 'sbg:category': Input Files
    id: Workflow_result
    type: File
    'sbg:fileTypes': TAR.GZ
outputs:
  - id: Match
    type: string
    outputBinding:
      loadContents: true
      glob: match.vcf
      outputEval: |-
        ${
            if (self[0].size) {
                return "Failure"
            } else {
                return "Success"
            }
        }
  - id: stderr
    type: File
    outputBinding:
      glob: log.stderr
  - id: stdout
    type: File
    outputBinding:
      glob: log.stdout
label: Validation_tool
arguments:
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          var comm ="{ "
          comm += "mkdir workflow_result expected_result"
          comm += " && tar -xzvf " + inputs.Workflow_result.path + " --strip-components 1 -C workflow_result/"
          comm += " && tar -xzvf " + inputs.Expected_result.path + " --strip-components 1 -C expected_result/"
          comm += " && mkdir match"
          comm += " && bcftools isec -p match expected_result/*.vcf.gz workflow_result/*.vcf.gz "
          comm += " && vcf-concat match/0000.vcf match/0001.vcf >> match2.vcf"
          comm += " && sed '/^#/d' match2.vcf >> match.vcf"
          comm += " && python validation.py; } >> log.stdout 2>> log.stderr"
          
          return comm
          
      }
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 100
    coresMin: 1
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/vladimir_obucina/topmed:validation_tool'
  - class: InitialWorkDirRequirement
    listing:
      - entryname: validation.py
        entry: |-
          import os
          import sys

          if os.stat("match.vcf").st_size == 0:
              sys.exit(0)
          else:
              sys.exit("Checking failed, output file is not as expected!")
  - class: InlineJavascriptRequirement
'sbg:createdBy': vladimir_obucina
'sbg:modifiedBy': vladimir_obucina
'sbg:revision': 15
'sbg:project': vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline
'sbg:image_url': >-
  https://igor.sbgenomics.com/ns/brood/images/vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/validation_tool/15.png
'sbg:modifiedOn': 1527691907
'sbg:validationErrors': []
'sbg:contributors':
  - vladimir_obucina
'sbg:revisionNotes': ''
'sbg:publisher': sbg
'sbg:sbgMaintained': false
'sbg:appVersion':
  - v1.0
'sbg:createdOn': 1527513179
'sbg:latestRevision': 15
'sbg:id': vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/validation_tool/15
'sbg:revisionsInfo':
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 0
    'sbg:revisionNotes': null
    'sbg:modifiedOn': 1527513179
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 1
    'sbg:revisionNotes': First version with bcftools and vcftools
    'sbg:modifiedOn': 1527522399
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 2
    'sbg:revisionNotes': ''
    'sbg:modifiedOn': 1527522451
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 3
    'sbg:revisionNotes': ''
    'sbg:modifiedOn': 1527522628
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 4
    'sbg:revisionNotes': ''
    'sbg:modifiedOn': 1527522900
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 5
    'sbg:revisionNotes': ''
    'sbg:modifiedOn': 1527523085
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 6
    'sbg:revisionNotes': ''
    'sbg:modifiedOn': 1527585970
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 7
    'sbg:revisionNotes': 'UPDATE: working locally'
    'sbg:modifiedOn': 1527587057
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 8
    'sbg:revisionNotes': ''
    'sbg:modifiedOn': 1527587374
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 9
    'sbg:revisionNotes': ''
    'sbg:modifiedOn': 1527588145
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 10
    'sbg:revisionNotes': ''
    'sbg:modifiedOn': 1527588853
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 11
    'sbg:revisionNotes': ''
    'sbg:modifiedOn': 1527595087
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 12
    'sbg:revisionNotes': Added Python script to check file size and provide exit status
    'sbg:modifiedOn': 1527672240
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 13
    'sbg:revisionNotes': >-
      UPDATE: Making log.stderr and log.stdout files, with exit message in case
      of failure.
    'sbg:modifiedOn': 1527690442
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 14
    'sbg:revisionNotes': 'UPDATE: Added stderr and stdout to output ports.'
    'sbg:modifiedOn': 1527690923
  - 'sbg:modifiedBy': vladimir_obucina
    'sbg:revision': 15
    'sbg:revisionNotes': ''
    'sbg:modifiedOn': 1527691907
$namespaces:
  sbg: 'https://sevenbridges.com'
'sbg:projectName': TOPMed Freeze 3a Variant Calling Pipeline

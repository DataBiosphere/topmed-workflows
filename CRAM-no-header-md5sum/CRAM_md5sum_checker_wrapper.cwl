cwlVersion: v1.0
class: Workflow

dct:creator:
  '@id':  https://orcid.org/0000-0001-5173-4627
  foaf:name: Walter Shands
  foaf:mbox: jshands@ucsc.edu

requirements:
  - class: SubworkflowFeatureRequirement

inputs:
  reference_file: File
  reference_index_file: File
  input_file: File
  input_index_file: File?
  truth_md5sum_file: File

outputs:
  workflow_output_file:
    type: File
    outputSource: checker/results_file

steps:
  md5sum:
    run: md5sum/CRAM_md5sum.cwl
    in:
      reference_file: reference_file
      reference_index_file: reference_index_file
      input_file: input_file
      input_index_file: input_index_file
    out: [output_file]


  checker:
    run: checker/CRAM_md5sum_checker.cwl
    in:
      input_file_1: md5sum/output_file
      input_file_2: truth_md5sum_file
    out: [results_file]

doc: |
  This wraps the md5sum tool with a checker workflow that runs both the tool and a tool that performs verification of results

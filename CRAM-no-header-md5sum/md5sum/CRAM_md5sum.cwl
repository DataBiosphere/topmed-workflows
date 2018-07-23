cwlVersion: v1.0
class: Workflow

dct:creator:
  '@id':  https://orcid.org/0000-0001-5173-4627
  foaf:name: Walter Shands
  foaf:mbox: jshands@ucsc.edu

inputs:
  reference_file: File
  reference_index_file: File
  input_file: File
  input_index_file: File?

outputs:
  output_file:
    type: File
    outputSource: CRAM_md5sum/output_file

steps:
  CRAM_samtools:
    run: CRAM_tool_samtools.cwl
    in:
      reference_file: reference_file
      reference_index_file: reference_index_file
      input_file: input_file
      input_index_file: input_index_file

    out:
      [output_file]

  CRAM_md5sum:
    run: CRAM_tool_md5sum.cwl
    in:
      input_file: CRAM_samtools/output_file

    out: [output_file]

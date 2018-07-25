#!/usr/bin/env cwl-runner

class: CommandLineTool
id: CRAM_md5sum
label: md5sum tool
cwlVersion: v1.0

$namespaces:
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/

dct:creator:
  '@id':  https://orcid.org/0000-0001-5173-4627
  foaf:name: Walter Shands
  foaf:mbox: jshands@ucsc.edu

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 512

requirements:
  DockerRequirement:
    dockerPull: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135"

inputs:
  input_file:
    type: stdin
    doc: The SAM file that will have its md5sum calculated.

outputs:
  output_file:
    type: stdout
    format: http://edamontology.org/data_3671
    doc: A text file containing the md5sum of the input SAM file.

stdout: $(inputs.input_file["basename"])_md5sum.txt
baseCommand: [md5sum]

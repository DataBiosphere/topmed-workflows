#!/usr/bin/env cwl-runner

class: CommandLineTool
id: Md5sumWorkflowChecker
label: A tool that checks the md5sum workflow
cwlVersion: v1.0

hints:
  ResourceRequirement:
    # The command really requires very little resources.
    coresMin: 1
    ramMin: 1024
    outdirMin: 512

requirements:
  DockerRequirement:
    dockerPull: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135"

inputs:
  input_file_1:
    type: File
    inputBinding:
      position: 1
    doc: The first file that contains an md5sum.

  input_file_2:
    type: File
    inputBinding:
      position: 2
    doc: The second file that contains an md5sum.

outputs:
  results_file:
    type: stdout
    doc: A file that contains the result of the diff.

stdout: md5sum_diff.txt
baseCommand: [diff]

#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
doc: |
    GATK Concordance tool to compare overlapping variants in two VCFs.
id: GATKConcordance

requirements:
  - class: ShellCommandRequirement
  - class: ScatterFeatureRequirement
  - class: DockerRequirement
    dockerPull: us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135
  - class: ResourceRequirement
    ramMin: 6000

baseCommand: [/usr/gitc/gatk4/gatk-launch, Concordance]

arguments:
  - valueFrom: "-XX:GCHeapFreeLimit=10
      -XX:+PrintFlagsFinal
      -XX:+PrintGCTimeStamps
      -XX:+PrintGCDateStamps
      -XX:+PrintGCDetails
      -Xloggc:gc_log.log
      -Xms4000m"
    position: 0
    prefix: "--javaOptions"
    separate: true

inputs:
 
  - id: reference
    type: File
    secondaryFiles: [^.dict, .fai]
    inputBinding:
      prefix: -R

  - id: eval
    type: File
    secondaryFiles: [.tbi]
    inputBinding:
      prefix: -eval

  - id: truth
    type: File
    secondaryFiles: [.tbi]
    inputBinding:
      prefix: --truth      

  - id: summary
    type: string
    inputBinding:
      prefix: --summary

outputs:
  - id: concordance_summary
    type: File
    outputBinding:
      glob: $(inputs.summary)


#!/usr/bin/env cwl-runner

class: Workflow

doc: |
    Workflow including checker to compare VCFs using GATK Concordance tool.

cwlVersion: v1.0

id: concordance-test-workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: ShellCommandRequirement

inputs:
 
  - id: reference
    type: File
    secondaryFiles: [^.dict, .fai]

  - id: eval
    type: File
    secondaryFiles: [.tbi]

  - id: truth
    type: File
    secondaryFiles: [.tbi]

  - id: summary
    type: string

  - id: threshold
    type: float

outputs:

  - id: concordance_summary
    type: File
    outputSource: GATKConcordance/concordance_summary

steps:

  - id: GATKConcordance
    run: ./GATKConcordance.yaml.cwl
    in:
       - id: reference
         source: reference
       - id: eval
         source: eval
       - id: truth
         source: truth
       - id: summary
         source: summary

    out:
       - id: concordance_summary

  - id: ConcordanceChecker
    run: ./GATKConcordanceChecker.yaml.cwl
    in:
       - id: summary
         source: GATKConcordance/concordance_summary
       - id: threshold
         source: threshold

    out: []

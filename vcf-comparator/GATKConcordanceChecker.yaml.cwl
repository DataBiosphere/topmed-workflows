#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: GATKConcordanceChecker
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement

inputs:
 
  - id: summary
    type: File
    inputBinding:
      loadContents: true

  - id: threshold
    type: float

outputs: []

baseCommand: []

arguments:
  - position: 0
    shellQuote: false
    valueFrom: |    
      ${
          var lines = inputs.summary.contents.split("\n");
          var flag = 0;
          var message = "The VCFs can be considered identical.";
          for (var i=1; i < lines.length; i++) {
              var sensitivity = lines[i].split("\t")[4];
              if (sensitivity < inputs.threshold) {
                  message = "The VCFs do not have enough overlapping.";
                  flag = 1;
              }
              var precision = parseFloat(lines[i].split("\t")[5]);
              if (precision < inputs.threshold) {
                  message = "The VCFs do not have enough overlapping.";
                  flag = 1;
              }
          }
          var output = "printf \"" + message + "\\n\"; " + "exit " + flag + ";";
          return output
      }

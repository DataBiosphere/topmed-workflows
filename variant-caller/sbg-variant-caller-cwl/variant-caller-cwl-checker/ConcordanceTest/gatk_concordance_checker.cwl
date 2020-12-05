class: CommandLineTool
cwlVersion: v1.0
id: >-
  vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/gatkconcordancechecker/2
baseCommand: []
inputs:
  - id: summary
    type: File
    inputBinding:
      position: 0
      loadContents: true
  - id: threshold
    type: float
outputs:
  - id: output
    type: string?
    outputBinding:
      glob: Summary
label: gatkconcordancechecker
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
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
'sbg:publisher': sbg
'sbg:modifiedOn': 1529067768
'sbg:sbgMaintained': false
'sbg:contributors':
  - vladimir_obucina
'sbg:appVersion':
  - v1.0
$namespaces:
  sbg: 'https://sevenbridges.com/'
'sbg:revision': 2
'sbg:revisionsInfo':
  - 'sbg:revision': 0
    'sbg:modifiedOn': 1529059731
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': null
  - 'sbg:revision': 1
    'sbg:modifiedOn': 1529059781
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': First Version
  - 'sbg:revision': 2
    'sbg:modifiedOn': 1529067768
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': Added output
'sbg:latestRevision': 2
'sbg:createdOn': 1529059731
'sbg:projectName': TOPMed Freeze 3a Variant Calling Pipeline
'sbg:validationErrors': []
'sbg:createdBy': vladimir_obucina
'sbg:image_url': >-
  https://igor.sbgenomics.com/ns/brood/images/vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/gatkconcordancechecker/2.png
'sbg:modifiedBy': vladimir_obucina
'sbg:revisionNotes': Added output
'sbg:id': >-
  vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/gatkconcordancechecker/2
'sbg:project': vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline

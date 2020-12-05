class: Workflow
cwlVersion: v1.0
id: >-
  vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/concordancecheckerworkflow/4
label: ConcordanceCheckerWorkflow
inputs:
  - id: threshold
    type: float
    'sbg:x': -105
    'sbg:y': -229
  - id: truth
    'sbg:fileTypes': 'VCF, VCF.GZ'
    type: File
    'sbg:x': -458
    'sbg:y': -238
  - id: summary
    type: string
    'sbg:x': -457
    'sbg:y': -113
  - id: reference
    type: File
    'sbg:x': -456
    'sbg:y': 8
  - id: eval
    'sbg:fileTypes': 'VCF, VCF.GZ'
    type: File
    'sbg:x': -455
    'sbg:y': 130
outputs:
  - id: concordance_summary
    outputSource:
      - gatkconcordance/concordance_summary
    type: File
    'sbg:x': 45
    'sbg:y': 42
  - id: output
    outputSource:
      - gatkconcordancechecker/output
    type: string?
    'sbg:x': 276.625
    'sbg:y': -132.5
steps:
  - id: gatkconcordance
    in:
      - id: reference
        source:
          - reference
      - id: eval
        source:
          - eval
      - id: truth
        source:
          - truth
      - id: summary
        source:
          - summary
    out:
      - id: concordance_summary
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: >-
        vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/gatkconcordance/1
      baseCommand:
        - /usr/gitc/gatk4/gatk-launch
        - Concordance
      inputs:
        - id: reference
          type: File
          inputBinding:
            position: 0
            prefix: '-R'
          secondaryFiles:
            - ^.dict
            - .fai
        - id: eval
          type: File
          inputBinding:
            position: 0
            prefix: '-eval'
          secondaryFiles:
            - .tbi
        - id: truth
          type: File
          inputBinding:
            position: 0
            prefix: '--truth'
          secondaryFiles:
            - .tbi
        - id: summary
          type: string
          inputBinding:
            position: 0
            prefix: '--summary'
      outputs:
        - id: concordance_summary
          type: File
          outputBinding:
            glob: $(inputs.summary)
      doc: |
        GATK Concordance tool to compare overlapping variants in two VCFs.
      label: GATKConcordance
      arguments:
        - position: 0
          prefix: '--javaOptions'
          valueFrom: >-
            -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps
            -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:gc_log.log
            -Xms4000m
      requirements:
        - class: ScatterFeatureRequirement
        - class: ResourceRequirement
          ramMin: 6000
        - class: DockerRequirement
          dockerPull: 'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135'
        - class: InlineJavascriptRequirement
      'sbg:publisher': sbg
      'sbg:revisionNotes': First Version
      'sbg:image_url': >-
        https://igor.sbgenomics.com/ns/brood/images/vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/gatkconcordance/1.png
      'sbg:sbgMaintained': false
      $namespaces:
        sbg: 'https://sevenbridges.com/'
      'sbg:modifiedBy': vladimir_obucina
      'sbg:appVersion':
        - v1.0
      'sbg:contributors':
        - vladimir_obucina
      'sbg:revision': 1
      'sbg:revisionsInfo':
        - 'sbg:revision': 0
          'sbg:modifiedBy': vladimir_obucina
          'sbg:modifiedOn': 1529059805
          'sbg:revisionNotes': null
        - 'sbg:revision': 1
          'sbg:modifiedBy': vladimir_obucina
          'sbg:modifiedOn': 1529059837
          'sbg:revisionNotes': First Version
      'sbg:latestRevision': 1
      'sbg:createdOn': 1529059805
      'sbg:projectName': TOPMed Freeze 3a Variant Calling Pipeline
      'sbg:validationErrors': []
      'sbg:createdBy': vladimir_obucina
      'sbg:modifiedOn': 1529059837
      'sbg:project': vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline
      'sbg:id': >-
        vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/gatkconcordance/1
    label: GATKConcordance
    'sbg:y': -97
    'sbg:x': -253
  - id: gatkconcordancechecker
    in:
      - id: summary
        source:
          - gatkconcordance/concordance_summary
      - id: threshold
        source:
          - threshold
    out:
      - id: output
    run:
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
      'sbg:image_url': >-
        https://igor.sbgenomics.com/ns/brood/images/vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/gatkconcordancechecker/2.png
      'sbg:sbgMaintained': false
      $namespaces:
        sbg: 'https://sevenbridges.com/'
      'sbg:modifiedBy': vladimir_obucina
      'sbg:appVersion':
        - v1.0
      'sbg:contributors':
        - vladimir_obucina
      'sbg:revision': 2
      'sbg:revisionsInfo':
        - 'sbg:revision': 0
          'sbg:modifiedBy': vladimir_obucina
          'sbg:modifiedOn': 1529059731
          'sbg:revisionNotes': null
        - 'sbg:revision': 1
          'sbg:modifiedBy': vladimir_obucina
          'sbg:modifiedOn': 1529059781
          'sbg:revisionNotes': First Version
        - 'sbg:revision': 2
          'sbg:modifiedBy': vladimir_obucina
          'sbg:modifiedOn': 1529067768
          'sbg:revisionNotes': Added output
      'sbg:latestRevision': 2
      'sbg:createdOn': 1529059731
      'sbg:projectName': TOPMed Freeze 3a Variant Calling Pipeline
      'sbg:validationErrors': []
      'sbg:createdBy': vladimir_obucina
      'sbg:revisionNotes': Added output
      'sbg:project': vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline
      'sbg:id': >-
        vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/gatkconcordancechecker/2
    label: GATKConcordanceChecker
    'sbg:y': -135
    'sbg:x': 114
requirements: []
'sbg:publisher': sbg
'sbg:modifiedOn': 1529919260
'sbg:image_url': >-
  https://igor.sbgenomics.com/ns/brood/images/vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/concordancecheckerworkflow/4.png
'sbg:sbgMaintained': false
'sbg:contributors':
  - vladimir_obucina
'sbg:appVersion':
  - v1.0
$namespaces:
  sbg: 'https://sevenbridges.com/'
'sbg:revision': 4
'sbg:revisionsInfo':
  - 'sbg:revision': 0
    'sbg:modifiedOn': 1529060084
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': null
  - 'sbg:revision': 1
    'sbg:modifiedOn': 1529063811
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': First Version
  - 'sbg:revision': 2
    'sbg:modifiedOn': 1529067109
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': Changed order in which tools are run
  - 'sbg:revision': 3
    'sbg:modifiedOn': 1529067834
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': Added output to checker
  - 'sbg:revision': 4
    'sbg:modifiedOn': 1529919260
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': ''
'sbg:latestRevision': 4
'sbg:createdOn': 1529060084
'sbg:projectName': TOPMed Freeze 3a Variant Calling Pipeline
'sbg:createdBy': vladimir_obucina
'sbg:modifiedBy': vladimir_obucina
'sbg:validationErrors': []
'sbg:revisionNotes': ''
'sbg:id': >-
  vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/concordancecheckerworkflow/4
'sbg:project': vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline

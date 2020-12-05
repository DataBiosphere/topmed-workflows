class: CommandLineTool
cwlVersion: v1.0
id: vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/gatkconcordance/1
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
label: gatkconcordance
arguments:
  - position: 0
    prefix: '--javaOptions'
    valueFrom: >-
      -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps
      -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:gc_log.log -Xms4000m
requirements:
  - class: ScatterFeatureRequirement
  - class: ResourceRequirement
    ramMin: 6000
  - class: DockerRequirement
    dockerPull: 'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135'
  - class: InlineJavascriptRequirement
'sbg:publisher': sbg
'sbg:modifiedOn': 1529059837
'sbg:contributors':
  - vladimir_obucina
'sbg:sbgMaintained': false
'sbg:appVersion':
  - v1.0
$namespaces:
  sbg: 'https://sevenbridges.com/'
'sbg:revision': 1
'sbg:revisionsInfo':
  - 'sbg:revision': 0
    'sbg:modifiedOn': 1529059805
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': null
  - 'sbg:revision': 1
    'sbg:modifiedOn': 1529059837
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': First Version
'sbg:latestRevision': 1
'sbg:createdOn': 1529059805
'sbg:projectName': TOPMed Freeze 3a Variant Calling Pipeline
'sbg:validationErrors': []
'sbg:createdBy': vladimir_obucina
'sbg:image_url': >-
  https://igor.sbgenomics.com/ns/brood/images/vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/gatkconcordance/1.png
'sbg:modifiedBy': vladimir_obucina
'sbg:revisionNotes': First Version
'sbg:id': vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/gatkconcordance/1
'sbg:project': vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline

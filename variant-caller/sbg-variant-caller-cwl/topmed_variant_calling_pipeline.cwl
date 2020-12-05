class: Workflow
cwlVersion: v1.0
id: >-
  vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/topmed-variant-calling-pipeline-cwl1/17
doc: >-
    This is the CWL wrapper for U of Michigan's [TOPMed Freeze 3a Variant Calling Pipeline](https://github.com/statgen/topmed_freeze3_calling)
label: TOPMed Variant Calling Pipeline CWL1
inputs:
  - id: reference
    'sbg:fileTypes': FA
    type: File
    'sbg:x': -460
    'sbg:y': -266.8571472167969
    secondaryFiles:
      - .fai
  - id: reference_file
    'sbg:fileTypes': TGZ
    type: File
    'sbg:x': -119
    'sbg:y': -629.1428833007812
  - id: pedigree_file
    'sbg:fileTypes': PED
    type: File?
    'sbg:x': -119.71428680419922
    'sbg:y': -515.8571166992188
  - id: num_of_jobs
    type: int?
    'sbg:x': -121.57142639160156
    'sbg:y': -401.71429443359375
  - id: genotype_unit
    type: int
    'sbg:x': -127
    'sbg:y': -50.85714340209961
  - id: discover_unit
    type: int
    'sbg:x': -133
    'sbg:y': 73
  - id: chromosomes
    type: 'string[]'
    'sbg:x': -128.14285278320312
    'sbg:y': 198.2857208251953
  - id: reference_genome
    type:
      type: enum
      symbols:
        - hg38
        - GRCh37
      name: reference_genome
    'sbg:x': -459.2857360839844
    'sbg:y': -385.71429443359375
  - id: bam_cram_file
    type: File
    'sbg:x': -548.7529296875
    'sbg:y': -25.529401779174805
    secondaryFiles:
      - |-
        ${
            return (self.basename + self.nameext.replace('m','i'))
        }
outputs:
  - id: called_variant_sites
    outputSource:
      - topmed_freeze3_calling/called_variant_sites
    type: File
    'sbg:x': 418.68548583984375
    'sbg:y': -43.3526725769043
  - id: genotypes
    outputSource:
      - topmed_freeze3_calling/genotypes
    type: File
    'sbg:x': 418.1114501953125
    'sbg:y': -197.72366333007812
  - id: makefile_log
    outputSource:
      - topmed_freeze3_calling/makefile_log
    type: File?
    'sbg:x': 423.53741455078125
    'sbg:y': -332.6839599609375
  - id: vcf_output
    outputSource:
      - topmed_freeze3_calling/vcf_output
    'sbg:fileTypes': GZ
    type: File?
    'sbg:x': 421.19287109375
    'sbg:y': -622.8525390625
  - id: vcf_index_output
    outputSource:
      - topmed_freeze3_calling/vcf_index_output
    'sbg:fileTypes': TBI
    type: 'File[]?'
    'sbg:x': 424.314697265625
    'sbg:y': -474.9278869628906
steps:
  - id: verifybamid_cwl1
    in:
      - id: bam_cram_file
        source:
          - bam_cram_file
      - id: reference
        source:
          - reference
      - id: reference_genome
        source:
          - reference_genome
    out:
      - id: output_index_file
    run: steps/verifybamid/verifybamid.cwl
    label: VerifyBamID_CWL1
    'sbg:x': -233.57144165039062
    'sbg:y': -197.14285278320312
  - id: topmed_freeze3_calling
    in:
      - id: bam_cram_file
        source:
          - bam_cram_file
      - id: chromosomes
        default: []
        source:
          - chromosomes
      - id: discover_unit
        source:
          - discover_unit
      - id: genotype_unit
        source:
          - genotype_unit
      - id: index_files
        source:
          - verifybamid_cwl1/output_index_file
      - id: num_of_jobs
        source:
          - num_of_jobs
      - id: pedigree_file
        source:
          - pedigree_file
      - id: reference_file
        source:
          - reference_file
      - id: reference_genome
        source:
          - reference_genome
    out:
      - id: called_variant_sites
      - id: genotypes
      - id: makefile_log
      - id: vcf_output
      - id: vcf_index_output
    run: steps/topmed_freeze3_calling/topmed_freeze3_calling.cwl
    label: Topmed_freeze3_CWL1
    'sbg:x': 157.14285278320312
    'sbg:y': -198
requirements: 
    - class: InlineJavascriptRequirement
$namespaces:
  sbg: 'https://sevenbridges.com/'


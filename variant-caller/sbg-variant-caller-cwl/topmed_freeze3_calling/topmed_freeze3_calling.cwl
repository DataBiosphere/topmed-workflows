class: CommandLineTool
cwlVersion: v1.0
id: >-
  vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/topmed_freeze3_calling/15
baseCommand: []
inputs:
  - format: 'BAI,CRAI'
    id: bai_crai_files
    type: 'File[]'
    label: BAI/CRAI files
    'sbg:fileTypes': 'BAI, CRAI'
  - format: 'BAM, CRAM'
    id: bam_cram_files
    type: 'File[]'
    inputBinding:
      position: 1
      shellQuote: false
      valueFrom: |-
        ${
            return ''
        }
    label: BAM/CRAM Files
    'sbg:fileTypes': 'BAM, CRAM'
  - id: chromosomes
    type: 'string[]'
    label: Chromosome
  - 'sbg:category': Input parameter
    id: discover_unit
    type: int
    label: DiscoverUnit
  - 'sbg:category': Input parameter
    id: genotype_unit
    type: int
    label: GenotypeUnit
  - format: INDEX
    'sbg:category': Input Files
    id: index_files
    type: 'File[]?'
    label: Index files
    doc: Index files with contamination values
    'sbg:fileTypes': INDEX
  - 'sbg:toolDefaultValue': '4'
    id: num_of_jobs
    type: int?
    label: Number of jobs
  - format: PED
    'sbg:category': Input Files
    id: pedigree_file
    type: File
    label: Pedigree File
    'sbg:fileTypes': PED
  - format: TGZ
    'sbg:category': Input Files
    id: reference_file
    type: File
    label: Reference File
    doc: Reference genome
    'sbg:fileTypes': TGZ
  - 'sbg:category': Input file
    id: reference_genome
    type:
      type: enum
      symbols:
        - hg38
        - GRCh37
      name: reference_genome
    label: Reference genome
outputs:
  - id: called_variant_sites
    label: Called variant sites
    type: File
    outputBinding:
      glob: called_variant_sites.tar.gz
  - id: genotypes
    label: Genotypes
    type: File
    outputBinding:
      glob: genotypes.tar.gz
  - id: makefile_log
    doc: Log when running make files
    type: File?
    outputBinding:
      glob: makefile.log
    format: LOG
label: Topmed_freeze3_CWL1
arguments:
  - position: 0
    prefix: ''
    separate: false
    shellQuote: false
    valueFrom: |-
      ${

          if (inputs.reference_genome == 'GRCh37') {
              var chr = inputs.chromosomes.join(' ')
          } else if (inputs.reference_genome == 'hg38'){
              for (var i = 0; i < inputs.chromosomes.length; i++) {
                  inputs.chromosomes[i] = "chr" + inputs.chromosomes[i];
              }
              var chr = inputs.chromosomes.join(' ')
          }
          var comm = "cat "

          for (var i = 0; i < inputs.index_files.length; i++) {
              comm += inputs.index_files[i].path + " ";
          }
          comm += "> trio_data.index && "
          /*
          for (var i = 0; i < inputs.bam_cram_files.length; i++) {
              comm += " ln -sf " + inputs.bam_cram_files[i].path + " /root/topmed_freeze3_calling/" + inputs.bam_cram_files[i].path.split("/").pop() + " && ";
          }
          for (var i = 0; i < inputs.bai_crai_files.length; i++) {
              comm += " ln -sf " + inputs.bai_crai_files[i].path + " /root/topmed_freeze3_calling/" + inputs.bai_crai_files[i].path.split("/").pop() + " && ";
          }
          */

          comm += "python /root/topmed_freeze3_calling/scripts/append_gcconfig.py -du " + inputs.discover_unit + " -gu " + inputs.genotype_unit + " -gen " + inputs.reference_genome
          comm += " && cp trio_data.index " + inputs.pedigree_file.path + " /root/topmed_freeze3_calling/data/ &&  cwd=\"$PWD\"  && cd /root/topmed_freeze3_calling/"
          comm += " && tar -xzvf " + inputs.reference_file.path + " -C /root/topmed_freeze3_calling/"

          // Step1
          comm += " && perl scripts/step1-detect-and-merge-variants.pl " + chr
          comm += " && sed -i '1s/^/SHELL:=\\/bin\\/bash\\n/' out/aux/Makefile"
          comm += " && make -f out/aux/Makefile"
          if (inputs.num_of_jobs) {
              comm += " -j " + inputs.num_of_jobs
          }
          comm += " >> makefile.log"

          // Step2
          comm += " && perl scripts/step2-joint-genotyping.pl " + chr
          comm += " && make -f out/paste/*.Makefile"
          if (inputs.num_of_jobs) {
              comm += " -j " + inputs.num_of_jobs
          }
          comm += " >> makefile.log"

          // Step3a
          comm += " && perl scripts/step3a-compute-milk-score.pl " + chr
          comm += " && make -f out/aux/milk/*.Makefile"
          if (inputs.num_of_jobs) {
              comm += " -j " + inputs.num_of_jobs
          }
          comm += " >> makefile.log"
          // Step3b
          comm += " && perl scripts/step3b-run-svm-milk-filter.pl " + chr
          comm += " >> makefile.log"
          comm += " && cd out"


          comm += " && tar -zcvf genotypes.tar.gz paste"
          comm += " && tar -zcvf called_variant_sites.tar.gz svm"
          comm += " && mv ../makefile.log genotypes.tar.gz called_variant_sites.tar.gz \"$cwd\""
          return comm
      }
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 8
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/vladimir_obucina/topmed:freeze3a'
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.pedigree_file)
      - entryname: gcconfig.pm
        entry: >+
          package gcconfig;


          use base qw/Exporter/;

          use Cwd qw(realpath);

          use File::Basename qw(dirname);

          use POSIX qw(pow sqrt);

          use FindBin;


          ## Variables and methods shared across the package

          our @EXPORT = qw($ref $md5 $bgzip $samtools $bcftools $bamUtil $tabix
          $index $pedf $out $vt $discoverUnit $genotypeUnit $omnivcf $hapmapvcf
          $dbsnp $invNorm $svmlearn $svmclassify $vcfsummary $vcfsummary2);


          ############################################################

          ### MODIFY THESE VARIABLES TO YOUR COMPUTING ENVIRONMENT

          our $index = "data/trio_data.index";

          our $pedf = "data/trio_data.ped";

          our $out = "out";

          ############################################################

          ### MODIFY THESE VARIABLES TO IF REFERENCE IS LOCATED ELSEWHERE

          our $refDir = "$FindBin::Bin/../gotcloud.ref";

          ############################################################

          ### MODIFY THESE VARIABLES TO IF EXTERNAL BINARIES ARE USED

          our $bgzip = "$FindBin::Bin/../htslib/bgzip";

          our $tabix = "$FindBin::Bin/../htslib/tabix";

          our $vt = "$FindBin::Bin/../vt/vt";

          our $samtools = "$FindBin::Bin/../samtools/samtools";

          our $bcftools = "$FindBin::Bin/../bcftools/bcftools";

          our $bamUtil = "$FindBin::Bin/../gotcloud/src/bin/bamUtil";

          our $svmlearn = "$FindBin::Bin/../gotcloud/src/bin/svm-train";

          our $svmclassify = "$FindBin::Bin/../gotcloud/src/bin/svm-predict";

          our $invNorm = "$FindBin::Bin/../gotcloud/src/bin/invNorm";

          our $vcfsummary = "$FindBin::Bin/vcf-summary";

          our $vcfsummary2 = "$FindBin::Bin/vcf-summary-v2";

      - entryname: append_gcconfig.py
        entry: >-
          import argparse



          parser = argparse.ArgumentParser()

          parser.add_argument("--discover_unit", "-du", help="discover unit
          size", type=str)

          parser.add_argument("--genotype_unit", "-gu", help="genotype unit
          size", type=str)

          parser.add_argument("--genome", "-gen", help="genome used", type=str)

          args = parser.parse_args()


          discover_unit = args.discover_unit

          genotype_unit = args.genotype_unit

          genome = args.genome


          f = open("/root/topmed_freeze3_calling/scripts/gcconfig.pm", "a")

          f.write("\n")

          if genome == "GRCh37":
              f.write("our $md5 = \"$refDir/md5/%2s/%2s/%s\";\n")
              f.write("our $ref = \"$refDir/hs37d5.fa\";\n")
              f.write("our $dbsnp = \"$refDir/dbsnp_142.b37.vcf.gz\";\n")
              f.write("our $hapmapvcf = \"$refDir/hapmap_3.3.b37.sites.vcf.gz\";\n")
              f.write("our $omnivcf = \"$refDir/1000G_omni2.5.b37.sites.PASS.vcf.gz\";\n")
          elif genome == "hg38":
              f.write("our $md5 = \"$refDir/md5/%2s/%2s/%s\";\n")
              f.write("our $ref = \"$refDir/hs38DH.fa\";\n")
              f.write("our $dbsnp = \"$refDir/dbsnp_142.b38.vcf.gz\";\n")
              f.write("our $hapmapvcf = \"$refDir/hapmap_3.3.b38.sites.vcf.gz\";\n")
              f.write("our $omnivcf = \"$refDir/1000G_omni2.5.b38.sites.PASS.vcf.gz\";\n")

          f.write("our $discoverUnit = " + discover_unit + ";\n")

          f.write("our $genotypeUnit = " + genotype_unit + ";\n1;")


          f.close()
  - class: InlineJavascriptRequirement
    expressionLib:
      - |-
        var updateMetadata = function(file, key, value) {
            file['metadata'][key] = value;
            return file;
        };


        var setMetadata = function(file, metadata) {
            if (!('metadata' in file))
                file['metadata'] = metadata;
            else {
                for (var key in metadata) {
                    file['metadata'][key] = metadata[key];
                }
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };

        var toArray = function(file) {
            return [].concat(file);
        };

        var groupBy = function(files, key) {
            var groupedFiles = [];
            var tempDict = {};
            for (var i = 0; i < files.length; i++) {
                var value = files[i]['metadata'][key];
                if (value in tempDict)
                    tempDict[value].push(files[i]);
                else tempDict[value] = [files[i]];
            }
            for (var key in tempDict) {
                groupedFiles.push(tempDict[key]);
            }
            return groupedFiles;
        };

        var orderBy = function(files, key, order) {
            var compareFunction = function(a, b) {
                if (a['metadata'][key].constructor === Number) {
                    return a['metadata'][key] - b['metadata'][key];
                } else {
                    var nameA = a['metadata'][key].toUpperCase();
                    var nameB = b['metadata'][key].toUpperCase();
                    if (nameA < nameB) {
                        return -1;
                    }
                    if (nameA > nameB) {
                        return 1;
                    }
                    return 0;
                }
            };

            files = files.sort(compareFunction);
            if (order == undefined || order == "asc")
                return files;
            else
                return files.reverse();
        };
hints:
  - class: 'sbg:useSbgFS'
    value: 'true'
'sbg:modifiedOn': 1527500690
'sbg:revision': 15
'sbg:image_url': >-
  https://igor.sbgenomics.com/ns/brood/images/vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/topmed_freeze3_calling/15.png
'sbg:revisionsInfo':
  - 'sbg:revision': 0
    'sbg:modifiedOn': 1525956103
    'sbg:modifiedBy': mikojicic
    'sbg:revisionNotes': null
  - 'sbg:revision': 1
    'sbg:modifiedOn': 1525956286
    'sbg:modifiedBy': mikojicic
    'sbg:revisionNotes': init
  - 'sbg:revision': 2
    'sbg:modifiedOn': 1525956409
    'sbg:modifiedBy': mikojicic
    'sbg:revisionNotes': name changed
  - 'sbg:revision': 3
    'sbg:modifiedOn': 1525957728
    'sbg:modifiedBy': mikojicic
    'sbg:revisionNotes': rev 71
  - 'sbg:revision': 4
    'sbg:modifiedOn': 1525957765
    'sbg:modifiedBy': mikojicic
    'sbg:revisionNotes': delete input
  - 'sbg:revision': 5
    'sbg:modifiedOn': 1525959177
    'sbg:modifiedBy': mikojicic
    'sbg:revisionNotes': javascript save
  - 'sbg:revision': 6
    'sbg:modifiedOn': 1525962377
    'sbg:modifiedBy': mikojicic
    'sbg:revisionNotes': use sbgfs
  - 'sbg:revision': 7
    'sbg:modifiedOn': 1525963786
    'sbg:modifiedBy': mikojicic
    'sbg:revisionNotes': default number of jobs fixed
  - 'sbg:revision': 8
    'sbg:modifiedOn': 1526053632
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': Removed file requirements
  - 'sbg:revision': 9
    'sbg:modifiedOn': 1526995956
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': Draft2 revision 98
  - 'sbg:revision': 10
    'sbg:modifiedOn': 1526996080
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': ''
  - 'sbg:revision': 11
    'sbg:modifiedOn': 1527013723
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': Removed making symlinks
  - 'sbg:revision': 12
    'sbg:modifiedOn': 1527068264
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': Added append_gcconfig.py script to create files
  - 'sbg:revision': 13
    'sbg:modifiedOn': 1527069675
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': ''
  - 'sbg:revision': 14
    'sbg:modifiedOn': 1527085513
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': Removed Index file from outputs.
  - 'sbg:revision': 15
    'sbg:modifiedOn': 1527500690
    'sbg:modifiedBy': vladimir_obucina
    'sbg:revisionNotes': 'UPDATE: GRCh37 insted of hg19'
$namespaces:
  sbg: 'https://sevenbridges.com'
'sbg:modifiedBy': vladimir_obucina
'sbg:publisher': sbg
'sbg:projectName': TOPMed Freeze 3a Variant Calling Pipeline
'sbg:revisionNotes': 'UPDATE: GRCh37 insted of hg19'
'sbg:appVersion':
  - v1.0
'sbg:latestRevision': 15
'sbg:sbgMaintained': false
'sbg:id': >-
  vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/topmed_freeze3_calling/15
'sbg:createdOn': 1525956103
'sbg:contributors':
  - mikojicic
  - vladimir_obucina
'sbg:validationErrors': []
'sbg:cmdPreview': >-
  cat /path/to/index_files-1.ext /path/to/index_files-2.ext > trio_data.index &&
  ln -sf /path/to/bam_cram_files-1.ext
  /root/topmed_freeze3_calling/bam_cram_files-1.ext &&  ln -sf
  /path/to/bam_cram_files-2.ext
  /root/topmed_freeze3_calling/bam_cram_files-2.ext &&  ln -sf
  /path/to/bai_crai_files-1.ext
  /root/topmed_freeze3_calling/bai_crai_files-1.ext &&  ln -sf
  /path/to/bai_crai_files-1.ext
  /root/topmed_freeze3_calling/bai_crai_files-1.ext && python
  /root/topmed_freeze3_calling/scripts/append_gcconfig.py -du 4 -gu 7 -gen hg38
  && cp trio_data.index /path/to/PedigreeFile.ext
  /root/topmed_freeze3_calling/data/ &&  cwd="$PWD"  && cd
  /root/topmed_freeze3_calling/ && tar -xzvf /path/to/ReferenceFile.ext -C
  /root/topmed_freeze3_calling/ && perl
  scripts/step1-detect-and-merge-variants.pl chrchromosome-string-value-1
  chrchromosome-string-value-2 && sed -i '1s/^/SHELL:=\/bin\/bash\n/'
  out/aux/Makefile && make -f out/aux/Makefile >> makefile.log && perl
  scripts/step2-joint-genotyping.pl chrchromosome-string-value-1
  chrchromosome-string-value-2 && make -f out/paste/*.Makefile >> makefile.log
  && perl scripts/step3a-compute-milk-score.pl chrchromosome-string-value-1
  chrchromosome-string-value-2 && make -f out/aux/milk/*.Makefile >>
  makefile.log && perl scripts/step3b-run-svm-milk-filter.pl
  chrchromosome-string-value-1 chrchromosome-string-value-2 >> makefile.log &&
  cd out && tar -zcvf genotypes.tar.gz paste && tar -zcvf
  called_variant_sites.tar.gz svm && mv ../makefile.log genotypes.tar.gz
  called_variant_sites.tar.gz "$cwd"
'sbg:createdBy': mikojicic
'sbg:project': vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline

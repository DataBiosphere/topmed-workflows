class: CommandLineTool
cwlVersion: v1.0
id: >-
  vladimir_obucina_topmed_freeze_3a_variant_calling_pipeline_topmed_freeze3_calling_29
baseCommand: []
inputs:
  - id: bam_cram_file
    type: File
    label: BAM/CRAM Files
    secondaryFiles:
      - |-
        ${
            return (self.basename + self.nameext.replace('m','i'))
        }
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
  - 'sbg:category': Input Files
    id: index_files
    type: File?
    label: Index files
    doc: Index files with contamination values
    'sbg:fileTypes': INDEX
  - 'sbg:toolDefaultValue': '4'
    id: num_of_jobs
    type: int?
    label: Number of jobs
  - 'sbg:category': Input Files
    id: pedigree_file
    type: File?
    label: Pedigree File
    'sbg:fileTypes': PED
  - 'sbg:category': Input Files
    id: reference_file
    type: File
    label: Reference File
    doc: Reference genome
    'sbg:fileTypes': TGZ
  - id: reference_genome
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
  - id: vcf_output
    doc: Output VCF file
    label: Output VCF file
    type: 'File[]?'
    outputBinding:
      glob: '*vcf.gz'
    'sbg:fileTypes': GZ
  - id: vcf_index_output
    doc: VCF index output file
    label: Output Index
    type: 'File[]?'
    outputBinding:
      glob: '*.tbi'
    'sbg:fileTypes': TBI
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

          comm += inputs.index_files.path + " ";
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
          if (inputs.pedigree_file != null) {
            comm += " && cp " + inputs.pedigree_file.path + " /root/topmed_freeze3_calling/data/"
          }
            
          comm += " && cp trio_data.index /root/topmed_freeze3_calling/data/"
          comm += " &&  cwd=\"$PWD\"  && cd /root/topmed_freeze3_calling/"
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
          if (inputs.pedigree_file) {
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
          }
          
          comm += " && cd out"
          comm += " && tar -zcvf genotypes.tar.gz paste"
          if (inputs.pedigree_file) {
              comm += " && tar -zcvf called_variant_sites.tar.gz svm"
              comm += " && mv called_variant_sites.tar.gz \"$cwd\""
          }
          comm += " && mv ../makefile.log genotypes.tar.gz paste/*vcf.gz paste/*.tbi \"$cwd\""
          //comm += " && mv ../makefile.log paste/*vcf.gz paste/*.tbi \"$cwd\""
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
$namespaces:
  sbg: 'https://sevenbridges.com/'

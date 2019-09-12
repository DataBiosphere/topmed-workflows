#! /bin/bash

# Set the exit code of a pipeline to that of the rightmost command
# to exit with a non-zero status, or zero if all commands of the pipeline exit
#set -o pipefail
# cause a bash script to exit immediately when a command fails
set ! -f
# cause the bash shell to treat unset variables as an error and exit immediately
set -u
# echo each line of the script to stdout so we can see what is happening
set -o xtrace
#to turn off echo do 'set +o xtrace'

   docker version
   cwltool --version
   java -version
   #java -jar cromwell-34.jar --version
   gsutil --version

   INPUT_DIR=inputs

   mkdir -p "$INPUT_DIR"


   [[ ! -f "$INPUT_DIR/truth_topmed_variant_caller_1_14_NWD176325_output.tar.gz" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller_checker/truth_topmed_variant_caller_1_14_NWD176325_output.tar.gz   ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/NWD176325.recab.cram" ]] && gsutil cp gs://topmed_workflow_testing/topmed_aligner/input_files/NWD176325.recab.cram  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/NWD176325.recab.cram.crai" ]] && gsutil cp gs://topmed_workflow_testing/topmed_aligner/input_files/NWD176325.recab.cram.crai  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1000G_omni2.5.b38.sites.PASS.vcf.gz" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1000G_omni2.5.b38.sites.PASS.vcf.gz  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1000G_omni2.5.b38.sites.PASS.vcf.gz.tbi" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1000G_omni2.5.b38.sites.PASS.vcf.gz.tbi  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr10.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr10.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr11_KI270927v1_alt.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr11_KI270927v1_alt.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr11.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr11.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr12.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr12.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr13.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr13.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr14_GL000009v2_random.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr14_GL000009v2_random.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr14_KI270846v1_alt.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr14_KI270846v1_alt.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr14.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr14.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr15.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr15.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr16.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr16.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr17_KI270857v1_alt.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr17_KI270857v1_alt.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr17_KI270862v1_alt.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr17_KI270862v1_alt.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr17_KI270909v1_alt.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr17_KI270909v1_alt.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr17.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr17.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr18.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr18.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr19_KI270938v1_alt.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr19_KI270938v1_alt.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr19.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr19.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr1_KI270706v1_random.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr1_KI270706v1_random.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr1_KI270766v1_alt.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr1_KI270766v1_alt.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr1.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr1.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr20.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr20.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr21.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr21.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr22_KI270879v1_alt.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr22_KI270879v1_alt.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr22_KI270928v1_alt.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr22_KI270928v1_alt.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr22.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr22.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr2_KI270773v1_alt.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr2_KI270773v1_alt.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr2_KI270894v1_alt.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr2_KI270894v1_alt.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr2.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr2.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr3.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr3.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr4_GL000008v2_random.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr4_GL000008v2_random.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr4.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr4.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr5.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr5.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr6.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr6.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr7_KI270803v1_alt.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr7_KI270803v1_alt.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr7.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr7.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr8_KI270821v1_alt.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr8_KI270821v1_alt.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr8.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr8.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chr9.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chr9.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chrUn_KI270742v1.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chrUn_KI270742v1.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/1kg.pilot_release.merged.indels.sites.hg38.chrX.vcf" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/1kg.pilot_release.merged.indels.sites.hg38.chrX.vcf  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/dbsnp_142.b38.vcf.gz" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/dbsnp_142.b38.vcf.gz  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/dbsnp_142.b38.vcf.gz.tbi" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/dbsnp_142.b38.vcf.gz.tbi  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/dbsnp.All.vcf.gz" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/dbsnp.All.vcf.gz  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/dbsnp.All.vcf.gz.tbi" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/dbsnp.All.vcf.gz.tbi  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hapmap_3.3.b38.sites.vcf.gz" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hapmap_3.3.b38.sites.vcf.gz  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hapmap_3.3.b38.sites.vcf.gz.tbi" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hapmap_3.3.b38.sites.vcf.gz.tbi  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH-bs.umfa" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH-bs.umfa  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.dict" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.dict  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa.alt" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.alt  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa.amb" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.amb  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa.ann" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.ann  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa.bwt" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.bwt  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa.fai" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.fai  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa.pac" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.pac  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa.sa" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.sa  ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.winsize100.gc" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.winsize100.gc  ./"$INPUT_DIR"

   pwd
   ls -al
   df -h
   #sudo du -hsx ./* | sort -n | head -100

   java -XX:+PrintCommandLineFlags
   java -Xms3G -Xmx7G -jar cromwell-34.jar run ../variant-caller/variant-caller-wdl-checker/topmed_freeze3_calling_checker.wdl -i ../variant-caller/variant-caller-wdl-checker/topmed_freeze3_calling_checker.wdl.local.json
   sudo rm -rf cromwell-executions

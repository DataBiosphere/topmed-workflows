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

   # if the reference genome file is not present then download it
   [[ ! -f "$INPUT_DIR/hs38DH.fa" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa.fai" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.fai ./"$INPUT_DIR"

   [[ ! -f "cromwell-34.jar" ]] && wget https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar

   [[ ! -f "$INPUT_DIR/hs38DH.fa.alt" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.alt ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa.bwt" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.bwt ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa.pac" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.pac ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa.ann" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.ann ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa.amb" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.amb ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa.sa" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.sa ./"$INPUT_DIR"

   [[ ! -f "$INPUT_DIR/NWD176325.0005.recab.cram" ]] && gsutil cp gs://topmed_workflow_testing/topmed_aligner/input_files/NWD176325.0005.recab.cram ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/NWD176325.0005.recab.cram.crai" ]] && gsutil cp gs://topmed_workflow_testing/topmed_aligner/input_files/NWD176325.0005.recab.cram.crai ./"$INPUT_DIR"

   [[ ! -f "$INPUT_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.gz" ]] && gsutil cp gs://topmed_workflow_testing/topmed_aligner/reference_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi" ]] && gsutil cp gs://topmed_workflow_testing/topmed_aligner/reference_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi ./"$INPUT_DIR"

   [[ ! -f "$INPUT_DIR/truth_NWD176325.0005.recab.cram" ]] && gsutil cp gs://topmed_workflow_testing/topmed_aligner_checker/truth_NWD176325.0005.recab.cram ./"$INPUT_DIR"

   pwd
   ls -al
   df -h
   #sudo du -hsx ./* | sort -n | head -100

   java -XX:+PrintCommandLineFlags
   java -Xms3G -Xmx7G -jar cromwell-34.jar run ../aligner/u_of_michigan_aligner-checker/u_of_michigan_aligner_checker.wdl -i ../aligner/u_of_michigan_aligner-checker/u_of_michigan_aligner_checker.local.json
   rm -rf cromwell-executions

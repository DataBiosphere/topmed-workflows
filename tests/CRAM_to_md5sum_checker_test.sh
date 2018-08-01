#! /bin/bash

# Set the exit code of a pipeline to that of the rightmost command
# to exit with a non-zero status, or zero if all commands of the pipeline exit
set -o pipefail
# cause a bash script to exit immediately when a command fails
set -e
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

   [[ ! -f "$INPUT_DIR/NWD119836.0005.recab.cram" ]] && gsutil cp gs://topmed_workflow_testing/topmed_aligner/input_files/NWD119836.0005.recab.cram ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/NWD119836.0005.recab.cram" ]] && gsutil cp gs://topmed_workflow_testing/topmed_aligner/input_files/NWD119836.0005.recab.cram.crai ./"$INPUT_DIR"

   [[ ! -f "$INPUT_DIR/NWD119836.0005.recab.cram.md5sum.txt" ]] && gsutil cp gs://topmed_workflow_testing/CRAM_to_md5sum/NWD119836.0005.recab.cram.md5sum.txt ./"$INPUT_DIR"

   pwd
   ls -al
   df -h
   #sudo du -hsx ./* | sort -n | head -100
   #java -XX:+PrintCommandLineFlags

   cwltool ../CRAM-no-header-md5sum/CRAM_md5sum_checker_wrapper.cwl ../CRAM-no-header-md5sum/CRAM_md5sum_checker_wrapper.cwl.local.json
   java -jar cromwell-34.jar run ../CRAM-no-header-md5sum/CRAM_md5sum_checker_wrapper.wdl -i ../CRAM-no-header-md5sum/CRAM_md5sum_checker_wrapper.wdl.local.json



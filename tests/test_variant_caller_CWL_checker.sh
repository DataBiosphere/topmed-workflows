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
   gsutil --version

   INPUT_DIR=inputs

   mkdir -p "$INPUT_DIR"

   # if the reference genome file is not present then download it
   [[ ! -f "$INPUT_DIR//NWD176325_genotypes.tar.gz" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller_checker/NWD176325_genotypes.tar.gz ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/NWD176325.0005.recab.cram" ]] && gsutil cp gs://topmed_workflow_testing/topmed_aligner/input_files/NWD176325.0005.recab.cram ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/NWD176325.0005.recab.cram.crai" ]] && gsutil cp gs://topmed_workflow_testing/topmed_aligner/input_files/NWD176325.0005.recab.cram.crai ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/"hs38DH-db142-v1.tgz ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH-db142-v1.tgz ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa ./"$INPUT_DIR"
   [[ ! -f "$INPUT_DIR/hs38DH.fa.fai" ]] && gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.fai ./"$INPUT_DIR"

   pwd
   ls -al
   df -h
   #sudo du -hsx ./* | sort -n | head -100

   echo "Running CWL aligner checker"
   cwltool --no-match-user --non-strict ../variant-caller/sbg-variant-caller-cwl/variant-caller-cwl-checker/topmed_freeze3_calling_checker.cwl ../variant-caller/sbg-variant-caller-cwl/variant-caller-cwl-checker/topmed_freeze3_calling_checker.local.json

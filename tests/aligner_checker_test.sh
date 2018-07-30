#! /bin/bash

# Set the exit code of a pipeline to that of the rightmost command
# to exit with a non-zero status, or zero if all commands of the pipeline exit
#set -o pipefail
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
   java -jar cromwell-34.jar --version
   gsutil --version

   # if the reference genome file is not present then download it
   if [ ! -f hs38DH.fa ]; then
     gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa .
   fi

   if [ ! -f hs38DH.fa.fai ]; then
     gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.fai .
   fi

   gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.alt .
   gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.bwt .
   gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.pac .
   gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.ann .
   gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.amb .
   gsutil cp gs://topmed_workflow_testing/topmed_variant_caller/reference_files/hg38/hs38DH.fa.sa .

   gsutil cp gs://topmed_workflow_testing/topmed_aligner/input_files/NWD176325.0005.recab.cram .
   gsutil cp gs://topmed_workflow_testing/topmed_aligner/input_files/NWD176325.0005.recab.cram.crai .

   gsutil cp gs://topmed_workflow_testing/topmed_aligner/reference_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz .
   gsutil cp gs://topmed_workflow_testing/topmed_aligner/reference_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi .

   gsutil cp gs://topmed_workflow_testing/topmed_aligner_checker/truth_NWD176325.0005.recab.cram .

   pwd
   ls -al
   df -h
   sudo du -hsx ./* | sort -n | head -100
   java -XX:+PrintCommandLineFlags

   travis_wait 55 java -Xms3G -Xmx5G -jar cromwell-34.jar run aligner/u_of_michigan_aligner-checker/u_of_michigan_aligner_checker.wdl -i aligner/u_of_michigan_aligner-checker/u_of_michigan_aligner_checker.local.json


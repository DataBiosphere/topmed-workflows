## This is the TopMed alignment workflow WDL for the workflow code located here:
## https://github.com/statgen/docker-alignment
##
## NOTE:
## The reference genome files to use are located here:
## ftp://share.sph.umich.edu/gotcloud/ref/hs38DH-db142-v1.tgz
##
## You can get the dbSNP reference files that the post align task uses at:
## gs://topmed-irc-share/resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz
## gs://topmed-irc-share/resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi
##

workflow TopMedAligner {

  File input_crai_file
  File input_cram_file

  String docker_image

  File ref_alt
  File ref_bwt
  File ref_pac
  File ref_ann
  File ref_amb
  File ref_sa

  File ref_fasta
  File ref_fasta_index

  File dbSNP_vcf
  File dbSNP_vcf_index

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  Int additional_disk = select_first([increase_disk_size, 20])

  # Sometimes the output is larger than the input, or a task can spill to disk. In these cases we need to account for the
  # input (1) and the output (1.5) or the input(1), the output(1), and spillage (.5).
  Float bwa_disk_multiplier = 2.5

  # Converting CRAM to fastq.gz takes extra disk space to store the fastq.gz files
  Float CRAM_to_fastqgz_multiplier = 2.5

  # Creating CRAM files from fastq.gz files increases the disk space needed
  Float fastq_gz_to_CRAM_multiplier = 1.5

  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data
  # so it needs more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a
  # larger multiplier
  Float sort_sam_disk_multiplier = 3.25

  # Get the size of the standard reference files as well as the additional reference files needed for BWA
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Float ref_extra_size = size(ref_alt, "GB") + size(ref_bwt, "GB") + size(ref_pac, "GB") + size(ref_ann, "GB") + size(ref_amb, "GB") + size(ref_sa, "GB")
  Float dbsnp_size = size(dbSNP_vcf, "GB") + size(dbSNP_vcf_index, "GB")
  Float cram_size = size(input_cram_file, "GB") + size(input_crai_file, "GB")
  Float fastq_gz_files_size = CRAM_to_fastqgz_multiplier * cram_size

  call PreAlign {
     input:
      input_crai = input_crai_file,
      input_cram = input_cram_file,
      disk_size = ref_size + (bwa_disk_multiplier * cram_size) + (sort_sam_disk_multiplier * cram_size) + cram_size + additional_disk + fastq_gz_files_size,
      docker_image = docker_image,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }

  call Align {
     input:
      input_list_file = PreAlign.output_list_file,
      input_fastq_gz_files = PreAlign.output_fastq_gz_files,

      disk_size = ref_size + ref_extra_size + (bwa_disk_multiplier * fastq_gz_files_size) + additional_disk,
      docker_image = docker_image,

      ref_alt = ref_alt,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_ann = ref_ann,
      ref_amb = ref_amb,
      ref_sa = ref_sa,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }

  Float CRAMS_files_size = fastq_gz_to_CRAM_multiplier * cram_size

  call PostAlign {
     input:
      input_cram_files = Align.output_cram_files,

      # The merged cram can be bigger than the summed sizes of the individual aligned crams,
      # so account for the output size by multiplying the input size by bwa disk multiplier.
      disk_size = ref_size + dbsnp_size + CRAMS_files_size + (sort_sam_disk_multiplier * CRAMS_files_size) + (bwa_disk_multiplier * CRAMS_files_size) + additional_disk,
      docker_image = docker_image,

      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,

      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index
  }

  output {
      File aligner_output = PostAlign.output_cram_file
  }
}

  task PreAlign {
     File input_crai
     File input_cram

     Float disk_size
     String docker_image

     File ref_fasta
     File ref_fasta_index

     # Assign a basename to the intermediate files
     String pre_output_base = "pre_output_base"

     command {

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

      echo "Running pre-alignment"

      samtools view -T ${ref_fasta} -uh -F 0x900 ${input_cram} \
        | bam-ext-mem-sort-manager squeeze --in -.ubam --keepDups --rmTags AS:i,BD:Z,BI:Z,XS:i,MC:Z,MD:Z,NM:i,MQ:i --out -.ubam \
        | samtools sort -l 1 -@ 1 -n -T ${pre_output_base}.samtools_sort_tmp - \
        | samtools fixmate - - \
        | bam-ext-mem-sort-manager bam2fastq --in -.bam --outBase ${pre_output_base} --maxRecordLimitPerFq 20000000 --sortByReadNameOnTheFly --readname --gzip

    }
     output {
      File output_list_file = "${pre_output_base}.list"
      # Capture all the files mentioned in the pre_output_base.list file
      # So they will be present for the Align task
      Array[File] output_fastq_gz_files = glob("${pre_output_base}.*")
    }
   runtime {
      memory: "10 GB"
      cpu: "32"
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      zones: "us-central1-a us-central1-b us-east1-d us-central1-c us-central1-f us-east1-c"
      docker: docker_image
    }
  }


  task Align {
     File input_list_file
     Array[File] input_fastq_gz_files

     Float disk_size
     String docker_image

     File ref_alt
     File ref_bwt
     File ref_pac
     File ref_ann
     File ref_amb
     File ref_sa

     File ref_fasta
     File ref_fasta_index

     # We have to use a trick to make Cromwell
     # skip substitution when using the bash ${<variable} syntax
     # This is necessary to get the <var>=$(<command>) sub shell
     # syntax to work and assign the value to a variable when
     # running in Cromwell
     # See https://gatkforums.broadinstitute.org/wdl/discussion/comment/44570#Comment_44570
     String dollar = "$"
     command <<<

      # Set the exit code of a pipeline to that of the rightmost command
      # to exit with a non-zero status, or zero if all commands of the pipeline exit
      # NOTE: Setting this will cause the pipeline to fail on Mac OS and Travis CI
      #       in some cases. It is commented out mainly so Travis CI will work.
      #       The failure was in samblaster
      #set -o pipefail
      # cause a bash script to exit immediately when a command fails
      set -e
      # cause the bash shell to treat unset variables as an error and exit immediately
      set -u
      # echo each line of the script to stdout so we can see what is happening
      set -o xtrace
      #to turn off echo do 'set +o xtrace'

      echo "Running alignment"

      # Get the Cromwell directory that is the input file location
      input_file_location=$(dirname ${input_fastq_gz_files[0]})

      while read line
      do
        line_rg=$(echo ${dollar}{line} | cut -d ' ' -f 4- | sed -e "s/ /\\\t/g")
        input_path=$(echo ${dollar}{line} | cut -f 2 -d ' ')
        input_filename=$(basename ${dollar}{input_path})
        output_filename=$(basename ${dollar}{input_filename} ".fastq.gz").cram

        # Prepend the path to the input file with the Cromwell input directory
        input_path=${dollar}{input_file_location}"/"${dollar}{input_filename}

        paired_flag=""
        if [[ ${dollar}{input_filename} =~ interleaved\.fastq\.gz$ ]]
        then
          paired_flag="-p"
        fi

        bwa mem -t 32 -K 100000000 -Y ${dollar}{paired_flag} -R ${dollar}{line_rg} ${ref_fasta} ${dollar}{input_path} | samblaster -a --addMateTags | samtools view -@ 32 -T ${ref_fasta} -C -o ${dollar}{output_filename} -
      done <<< "$(tail -n +2 ${input_list_file})"

    >>>
     output {
      Array[File] output_cram_files = glob("*.cram")
    }
   runtime {
      memory: "10 GB"
      cpu: "32"
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      zones: "us-central1-a us-central1-b us-east1-d us-central1-c us-central1-f us-east1-c"
      docker: docker_image
    }
  }

task PostAlign {
     Float disk_size
     String docker_image

     File ref_fasta
     File ref_fasta_index

     File dbSNP_vcf
     File dbSNP_vcf_index

     Array[File] input_cram_files

     # We have to use a trick to make Cromwell
     # skip substitution when using the bash ${<variable} syntax
     # This is necessary to get the <var>=$(<command>) sub shell 
     # syntax to work and assign the value to a variable when 
     # running in Cromwell
     # See https://gatkforums.broadinstitute.org/wdl/discussion/comment/44570#Comment_44570 
     String dollar = "$"

     command <<<
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

      echo "Running post alignment"

      # Get the Cromwell directory that is the input file location
      input_file_location=$(dirname ${input_cram_files[0]})

      rc=0
      for input_file in ${dollar}{input_file_location}"/"*.cram 
      do 
        # Put the output file in the local Cromwell working dir
        input_base_file_name=$(basename ${dollar}{input_file} ".cram")
        tmp_prefix=${dollar}{input_base_file_name}.tmp
        samtools sort --reference ${ref_fasta} --threads 1 -T $tmp_prefix -o ${dollar}{input_base_file_name}.sorted.bam ${dollar}{input_file}

#        tmp_prefix=${dollar}{input_file%.cram}.tmp
#        samtools sort --reference ${ref_fasta} --threads 1 -T $tmp_prefix -o ${dollar}{input_file%.cram}.sorted.bam ${dollar}{input_file}

        rc=$?
        [[ $rc != 0 ]] && break
#        rm -f ${dollar}{input_file} ${dollar}{tmp_prefix}*
        # Remove the tmp file; no need to remove the input file from the previous task
        rm -f ${dollar}{tmp_prefix}*
      done

      if [[ $rc == 0 ]]
      then 
        samtools merge --threads 1 -c merged.bam *.sorted.bam \
          && rm ./*.sorted.bam \
          && bam-non-primary-dedup dedup_LowMem --allReadNames --binCustom --binQualS 0:2,3:3,4:4,5:5,6:6,7:10,13:20,23:30 --log dedup_lowmem.metrics --recab --in merged.bam --out -.ubam --refFile ${ref_fasta} --dbsnp ${dbSNP_vcf} \
          | samtools view -h -C -T ${ref_fasta} -o output_file.cram --threads 1
        rc=$?
      fi
    >>>
     output {
      File output_cram_file = "output_file.cram"
    }
   runtime {
      memory: "10 GB"
      cpu: "32"
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      zones: "us-central1-a us-central1-b us-east1-d us-central1-c us-central1-f us-east1-c"
      docker: docker_image
    }
  }


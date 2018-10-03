workflow discoverAndMergeVariants {

  File? input_crai_file
  File input_cram_file

  File ref_fasta
  File ref_fasta_index

  String docker_image

  Int? CPUs
  Int CPUs_default = select_first([CPUs, 1])

  Int? memory
  Int memory_default = select_first([memory, 7])

  Int? preemptible_tries
  Int preemptible_tries_default = select_first([preemptible_tries, 3])

  Int? max_retries
  Int max_retries_default = select_first([max_retries, 3])

  Boolean? dynamically_calculate_file_size
  Boolean dynamically_calculate_disk_requirement = select_first([dynamically_calculate_file_size, true])

  Float? CRAM_file_max_disk_size_override
  Float CRAM_file_max_disk_size_override_default = select_first([CRAM_file_max_disk_size_override, 200.0])

  Float? CRAI_file_max_disk_size_override
  Float CRAI_file_max_disk_size_override_default = select_first([CRAI_file_max_disk_size_override, 10.0])

  Float? ReferenceGenome_disk_size_override
  Float ReferenceGenome_disk_size_override_default = select_first([ReferenceGenome_disk_size_override, 4.0])

  Float? ReferenceGenome_index_disk_size_override
  Float ReferenceGenome_index_disk_size_override_default = select_first([ReferenceGenome_index_disk_size_override, 1.0])


  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  Int additional_disk = select_first([increase_disk_size, 20])

  Float reference_size = if(dynamically_calculate_disk_requirement) then size(ref_fasta, "GB") + size(ref_fasta_index, "GB") else ReferenceGenome_index_disk_size_override_default + ReferenceGenome_disk_size_override_default

  # Account for size of generated CRAM index file
  Float cram_size = if(dynamically_calculate_disk_requirement) then size(input_cram_file, "GB") else CRAM_file_max_disk_size_override_default
  Float crai_size = if(dynamically_calculate_disk_requirement) then (if (defined(input_crai_file)) then size(input_crai_file, "GB") else (cram_size * 0.00003)) else CRAI_file_max_disk_size_override_default


  Array[String] detect_and_merge_targets_list
  File detect_and_merge_Makefile
  File gcconfig_pm
  File config_pm
  File trio_data_index


  scatter(target in detect_and_merge_targets_list) {

      # Get the Cromwell basename  of the CRAM file
      # This should have been the ID that is put in
      # the trio_data.index file and is present in the
      # target name as follows:
      # 'out/aux/individual/<ID>/<chr#>_<range>.sites.bcf.OK'
      #String base_name_wo_extension = basename(input_cram_file, ".cram")
      String base_name_wo_extension = sub(basename(input_cram_file), "\\..*$", "")
      # If substiting what the target name should be (accepting any region)
      # with an empty string results
      # in an empty string then this is a target for this CRAM
      #Boolean cram_target_valid = sub(target, ".*\\/out\\/aux\\/individual\\/${base_name_wo_extension}\\/.*.sites.bcf.OK", "") == ""
      String search_string = ".*out\\/aux\\/individual\\/${base_name_wo_extension}/.*\\.sites\\.bcf\\.OK"
      Boolean cram_target_valid = sub(target, search_string, "") == ""

      if ( cram_target_valid ) {
          call runDiscoveryAndMergeTarget as scatter_runDiscoveryAndMergeTarget {
              input:
                  input_cram = input_cram_file,
                  input_crai = input_crai_file,
                  ref_fasta = ref_fasta,
                  ref_fasta_index = ref_fasta_index,

                  target = target,
                  sample_id = base_name_wo_extension,
                  trio_data_index = trio_data_index,
                  gcconfig_pm = gcconfig_pm,
                  config_pm = config_pm,
                  detect_and_merge_Makefile = detect_and_merge_Makefile,


                  disk_size = cram_size + crai_size + reference_size + additional_disk,
                  memory = memory_default,
                  CPUs = CPUs_default,
                  preemptible_tries = preemptible_tries_default,
                  max_retries = max_retries_default,
                  docker_image = docker_image
          }
      }

  }

  Array[Pair[String, File]?] discovery_ID_to_BCF_file = scatter_runDiscoveryAndMergeTarget.discovery_ID_to_BCF_file
  Array[Pair[String, File]?] discovery_ID_to_log_file = scatter_runDiscoveryAndMergeTarget.discovery_ID_to_log_file

  output {
      Array[Pair[String, File]?] discovery_ID_to_BCF_file_output = discovery_ID_to_BCF_file
      Array[Pair[String, File]?] discovery_ID_to_log_file_output = discovery_ID_to_log_file
  }
}

  task runDiscoveryAndMergeTarget {
     File input_cram
     File? input_crai

     File ref_fasta
     File ref_fasta_index

     File detect_and_merge_Makefile
     File gcconfig_pm
     File config_pm
     File trio_data_index
     String target
     String sample_id

     Float memory
     Float disk_size
     Int CPUs
     Int preemptible_tries
     Int max_retries
     String docker_image

     String CRAM_basename = basename(input_cram)
     String output_crai_file_name = "${CRAM_basename}.crai"

     String output_BCF_file_name = sub(target, "\\.sites\\.bcf\\.OK", "\\.sites\\.bcf")
     String output_log_file_name = sub(target, "\\.sites\\.bcf\\.OK", "\\.discover2\\.log")


     # We have to use a trick to make Cromwell
     # skip substitution when using the bash ${<variable} syntax
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

      echo "Running variant discovery ${target} for ${input_cram}"

      # If there is no CRAM index file generate it
      if [[ -z "${input_crai}" ]]
      then
          printf "Creating index for ${input_cram}"
          /root/topmed_freeze3_calling/samtools/samtools index ${input_cram}
      fi

      # Symlink the config files to where the variant caller scripts expect them to be
      ln -s ${gcconfig_pm} /root/topmed_freeze3_calling/scripts/gcconfig.pm
      ln -s ${config_pm} /root/topmed_freeze3_calling/scripts/config.pm
      ln -s ${detect_and_merge_Makefile} /root/topmed_freeze3_calling/out/aux/Makefile
      ln -s ${trio_data_index} /root/topmed_freeze3_calling/data/trio_data.index

      ln -s ${ref_fasta}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hs38DH.fa
      ln -s ${ref_fasta_index}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hs38DH.fa.fai

      # Run discovery and merge on the target
      make SHELL='/bin/bash' -f /root/topmed_freeze3_calling/out/aux/Makefile ${target}

      >>>
        output {
          Pair[String, File] discovery_ID_to_BCF_file = (sample_id, "${output_BCF_file_name}")
          Pair[String, File] discovery_ID_to_log_file = (sample_id, "${output_log_file_name}")
       }
      runtime {
         memory: sub(memory, "\\..*", "") + " GB"
         cpu: sub(CPUs, "\\..*", "")
         maxRetries: max_retries
         preemptible: preemptible_tries
         disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
         zones: "us-central1-a us-central1-b us-east1-d us-central1-c us-central1-f us-east1-c"
         docker: docker_image
       }
  }



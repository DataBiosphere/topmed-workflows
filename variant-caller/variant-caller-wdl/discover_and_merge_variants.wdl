workflow discoverAndMergeVariants {

  File? input_crai_file
  File? input_cram_file
  Array[Map[String, Array[File]]]? BCFFiles

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
  Float cram_size = if(dynamically_calculate_disk_requirement) then (if (defined(input_cram_file)) then size(input_cram_file, "GB") else CRAM_file_max_disk_size_override_default) else CRAM_file_max_disk_size_override_default
  Float crai_size = if(dynamically_calculate_disk_requirement) then (if (defined(input_crai_file)) then size(input_crai_file, "GB") else (cram_size * 0.00003)) else CRAI_file_max_disk_size_override_default


  Array[String] detect_and_merge_targets_list
  File detect_and_merge_Makefile
  File gcconfig_pm
  File config_pm
  File trio_data_index

#  scatter(target in detect_and_merge_targets_list) {
#
#      # Get the Cromwell basename  of the CRAM file
#      # This should have been the ID that is put in
#      # the trio_data.index file and is present in the
#      # target name as follows:
#      # 'out/aux/individual/<ID>/<chr#>_<range>.sites.bcf.OK'
#      #String base_name_wo_extension = basename(input_cram_file, ".cram")
#      String base_name_wo_extension = sub(basename(input_cram_file), "\\..*$", "")
#      # If substiting what the target name should be (accepting any region)
#      # with an empty string results
#      # in an empty string then this is a target for this CRAM
#      #Boolean cram_target_valid = sub(target, ".*\\/out\\/aux\\/individual\\/${base_name_wo_extension}\\/.*.sites.bcf.OK", "") == ""
#      String search_string = ".*out\\/aux\\/individual\\/${base_name_wo_extension}/.*\\.sites\\.bcf\\.OK"
#      Boolean cram_target_valid = sub(target, search_string, "") == ""
#
#      if ( cram_target_valid ) {
#          call runDiscoveryAndMergeTarget as scatter_runDiscoveryAndMergeTarget {
#              input:
#                  input_cram = input_cram_file,
#                  input_crai = input_crai_file,
#                  ref_fasta = ref_fasta,
#                  ref_fasta_index = ref_fasta_index,
#
#                  target = target,
#                  sample_id = base_name_wo_extension,
#                  trio_data_index = trio_data_index,
#                  gcconfig_pm = gcconfig_pm,
#                  config_pm = config_pm,
#                  detect_and_merge_Makefile = detect_and_merge_Makefile,
#
#
#                  disk_size = cram_size + crai_size + reference_size + additional_disk,
#                  memory = memory_default,
#                  CPUs = CPUs_default,
#                  preemptible_tries = preemptible_tries_default,
#                  max_retries = max_retries_default,
#                  docker_image = docker_image
#          }
#      }
#
#  }
#
#  Array[Pair[String, File]?] discovery_ID_to_BCF_file = scatter_runDiscoveryAndMergeTarget.discovery_ID_to_BCF_file
#  Array[Pair[String, File]?] discovery_ID_to_log_file = scatter_runDiscoveryAndMergeTarget.discovery_ID_to_log_file
#
#  output {
#      Array[Pair[String, File]?] discovery_ID_to_BCF_file_output = discovery_ID_to_BCF_file
#      Array[Pair[String, File]?] discovery_ID_to_log_file_output = discovery_ID_to_log_file
#  }


  # Get the Cromwell basename  of the CRAM file
  # This should have been the ID that is put in
  # the trio_data.index file and is present in the
  # target name as follows:
  # 'out/aux/individual/<ID>/<chr#>_<range>.sites.bcf.OK'
  String base_name_wo_extension = if (defined(input_cram_file)) then sub(basename(input_cram_file), "\\..*$", "") else ""

  call runDiscoverVariants {
          input:
              sampleBCFs = BCFFiles,
              input_cram = input_cram_file,
              input_crai = input_crai_file,
              ref_fasta = ref_fasta,
              ref_fasta_index = ref_fasta_index,

              all_sample_targets = detect_and_merge_targets_list,
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

  output {
      Map[String, Array[File]] discovery_ID_to_BCF_file_output = runDiscoverVariants.discovery_ID_to_BCF_files
      #Array[File] discovery_ID_to_BCF_file_output = runDiscoverVariants.discovery_ID_to_BCF_files

      #Pair[String, Array[File]] discovery_ID_to_BCF_file_output = runDiscoverVariants.discovery_ID_to_BCF_files
      #Pair[String, Array[File]] discovery_ID_to_log_file_output = runDiscoverVariants.discovery_ID_to_log_files
  }
}


  task runDiscoverVariants {
     File? input_cram
     File? input_crai
     Array[Map[String,Array[File]]]? sampleBCFs

     File ref_fasta
     File ref_fasta_index

     File detect_and_merge_Makefile
     File gcconfig_pm
     File config_pm
     File trio_data_index
     Array[String] all_sample_targets
     String sample_id

     Float memory
     Float disk_size
     Int CPUs
     Int preemptible_tries
     Int max_retries
     String docker_image

     #String CRAM_basename = basename(input_cram)
     #String output_crai_file_name = "${CRAM_basename}.crai"

     #String output_BCF_files = "out/aux/individual/${sample_id}/*.sites.bcf"
     #String output_log_files = "out/aux/individual/${sample_id}/*.sites.log"
     String output_BCF_files = "out/aux/individual/${sample_id}/*"
     #String output_log_files = "out/aux/individual/${sample_id}/*.sites.log"
     String output_merged_BCF_files = "out/aux/union/*"


     # We have to use a trick to make Cromwell
     # skip substitution when using the bash ${<variable} syntax
     # See https://gatkforums.broadinstitute.org/wdl/discussion/comment/44570#Comment_44570
     String dollar = "$"

     command <<<

      # Make the log directory so the variant caller can output the logs there
      mkdir -p out/log
      # Make the aux directorys so the variant caller can create the sample dirs there
      mkdir -p out/aux/individual
      mkdir -p out/aux/union
      mkdir -p out/aux/sites
      mkdir -p out/aux/evaluation
      mkdir -p out/paste



      python3.5 <<CODE

      import csv
      import os
      from shutil import copy
      import sys
      import errno
      import json
      import errno

      if len("${input_cram}") == 0:
          all_BCFS_json_string = ""
          all_BCFs_json = []
          with open("${write_json(sampleBCFs)}", 'r') as all_BCFs_json_file:
              all_BCFs_json = json.load(all_BCFs_json_file)
              print("BCF JSON is:{}".format(all_BCFs_json))
          #all_BCFs_json = json.loads(all_BCFs_json_string)
          for sample_BCF_output in all_BCFs_json:
              for sampleID, BCF_output_array in sample_BCF_output.items():
                  for BCF_file in BCF_output_array:
                      BCF_file_basename = os.path.basename(BCF_file)
                      symlink_path = "out/aux/individual/" + sampleID + "/" + BCF_file_basename

                      # Create the directory to hold the BCFs for the sample
                      # and don't throw an exception if it already exists
                      # https://stackoverflow.com/questions/16029871/how-to-run-os-mkdir-with-p-option-in-python
                      directory_name = os.path.dirname(symlink_path)
                      print("Creating directory {}".format(directory_name))
                      try:
                          os.makedirs(directory_name)
                      except OSError as exc:
                          if exc.errno == errno.EEXIST and os.path.isdir(directory_name):
                              pass

                      print("Creating symlink for {} as {}".format(BCF_file, symlink_path))
                      os.symlink(BCF_file, symlink_path)







      CODE



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

      printf "Running variant discovery"

      # If a CRAM file was input
      if [[ -n "${input_cram}" ]]; then

          printf "Creating symlink for ${input_cram}"
          cram_basename=$(basename "${input_cram}")
          # Symlink the CRAM file so that the script can find it
          # We assume the CRAM name was put into the trio_data.index file as
          # the basename of the input CRAM file in the setupConfigFiles task
          # The CRAI file will be in the current working dir whether
          # it was input or created in this task
          ln -s ${input_cram} ${dollar}{cram_basename}

          output_crai_file_name=${dollar}{cram_basename}.crai
          # If there is no CRAM index file generate it
          if [[ -z "${input_crai}" ]]
          then
              printf "Creating index for ${input_cram}"
              /root/topmed_freeze3_calling/samtools/samtools index ${input_cram}
              # Symlink the CRAI file from the current working directory to the
              # input dir where samtools creates it
              # based on the location of the input CRAM file
              # The trio_data.index file lists the CRAM in the current working
              # directory so samtools expects the CRAI to be in the CWD
              ln -s ${input_cram}.crai ${dollar}{output_crai_file_name}
          else
              ln -s ${input_crai} ${dollar}{output_crai_file_name}
          fi
      fi


      # Symlink the config files to where the variant caller scripts expect them to be
      # Use the -f switch to unlink the existing config files so we can create the link
      # to the config files with the variant caller updated paths
      ln -fs ${gcconfig_pm} /root/topmed_freeze3_calling/scripts/gcconfig.pm
      ln -fs ${config_pm} /root/topmed_freeze3_calling/scripts/config.pm
      ln -fs ${trio_data_index} /root/topmed_freeze3_calling/data/trio_data.index

      # Make sure the directory where the reference files are supposed to be
      # located exists in the container
      mkdir -p /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38

      ln -s ${ref_fasta}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hs38DH.fa
      ln -s ${ref_fasta_index}  /root/topmed_freeze3_calling/data/local.org/ref/gotcloud.ref/hg38/hs38DH.fa.fai

      # Make the log directory so the variant caller can output the logs there
      mkdir -p out/log
      # Make the aux directorys so the variant caller can create the sample dirs there
      mkdir -p out/aux/individual/${sample_id}
      mkdir -p out/aux/union
      mkdir -p out/aux/sites
      mkdir -p out/aux/evaluation
      mkdir -p out/paste

      printf "Search ${ sep=',' all_sample_targets } for matches to ${sample_id}"
      # Get the list of Makefile targets for the input CRAM file
      ALL_CRAM_MAKEFILE_TARGETS="${ sep=',' all_sample_targets }"
      #CRAM_MAKEFILE_TARGETS=$(grep -o "^.*\/out\/aux\/individual\/${sample_id}\/chr[X_0-9]*.sites.bcf.OK" ${dollar}{ALL_CRAM_MAKEFILE_TARGETS})


      # If there is no CRAM file then variant discovery should
      # have been done on all CRAMs and the next step is to 
      # merge the results. In this case we need to run the Makefile
      # targets that put the merged variants in the 'union' folder
      if [[ -z "${input_cram}" ]]
      then
          REG_EX="out\/aux\/union\/chr[X_0-9]*\.sites\.bcf(\.csi)?\.OK"
      else
          REG_EX=".*out\/aux\/individual\/"${sample_id}"\/chr[X_0-9]*\.sites\.bcf\.OK"
      fi

      #REG_EX=".*out\/aux\/individual\/"${sample_id}"\/chr[X_0-9]*.sites.bcf.OK"
      CRAM_MAKEFILE_TARGETS=$(echo ${dollar}{ALL_CRAM_MAKEFILE_TARGETS} | awk -v pat="$REG_EX" 'BEGIN{RS=","} {where = match($1, pat); if (where !=0 ) printf "%s ",$1 }')

      printf "Running discovery with targets: ${dollar}{CRAM_MAKEFILE_TARGETS}"
      # Run discovery and merge on the target
      make SHELL='/bin/bash' -f ${detect_and_merge_Makefile} ${dollar}{CRAM_MAKEFILE_TARGETS}

      >>>
        output {
          Map[String, Array[File]] discovery_ID_to_BCF_files = if (defined("${input_cram}")) then {sample_id : glob("${output_BCF_files}")} else {sample_id : glob("${output_merged_BCF_files}")}

          #Array[File] discovery_ID_to_BCF_files = if (defined("${input_cram}")) then glob("${output_BCF_files}") else glob("${output_merged_BCF_files}")

          #Pair[String, Array[File]] discovery_ID_to_BCF_files = (sample_id, glob("${output_BCF_files}"))
          #Pair[String, Array[File]] discovery_ID_to_log_files = (sample_id, glob("${output_log_files}"))
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



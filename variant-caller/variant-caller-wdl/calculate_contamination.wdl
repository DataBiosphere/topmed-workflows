## Calculates DNA contamination by using the Docker image built from the Dockerfile
# at https://github.com/DataBiosphere/topmed-workflows/blob/<version>/variant-caller/variant-caller-wdl/verifybamid/Dockerfile
## which is based on the VerifyBamID tool at https://github.com/griffan/VerifyBamID

workflow calulateDNAContamination {

  File? input_crai_file
  File input_cram_file

  File ref_fasta
  File ref_fasta_index

  Int? CalcContamination_CPUs
  Int CPUs_default = select_first([CalcContamination_CPUs, 1])

  Int? CalcContamination_memory
  Int memory_default = select_first([CalcContamination_memory, 7])

  Int? CalcContamination_preemptible_tries
  Int preemptible_tries_default = select_first([CalcContamination_preemptible_tries, 3])

  Int? CalcContamination_max_retries
  Int max_retries_default = select_first([CalcContamination_max_retries, 3])

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

  String? reference_genome_version
  String reference_genome = select_first([reference_genome_version, 'hg38'])

  String? docker_image = "quay.io/ucsc_cgl/verifybamid:1.30.0"

  call VerifyBamID {
     input:
      input_crai = input_crai_file,
      input_cram = input_cram_file,

      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,

      reference_genome = reference_genome,
      preemptible_tries = preemptible_tries_default,

      memory = memory_default,
      CPUs = CPUs_default,
      max_retries = max_retries_default,
      disk_size = cram_size + crai_size + reference_size + additional_disk,
      docker_image = docker_image
  }

  output {
      Array[String] calculate_DNA_contamination_output = VerifyBamID.DNA_contamination_output_files

  }
}

  task VerifyBamID {
     File? input_crai
     File input_cram

     File ref_fasta
     File ref_fasta_index

     String reference_genome

     Int preemptible_tries
     Int CPUs
     Int memory
     Float disk_size
     Int max_retries
     String docker_image

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

      echo "Running VerifyBamID"

      # If there is no CRAM index file generate it
      if [[ -z "${input_crai}" ]]
      then 
          printf "Creating index for ${input_cram}"
          samtools index ${input_cram}
      fi

      if [[ ${reference_genome} == 'hg37' ]]
      then
             printf "VerifyBamID: Using hg37 genome\n"
             UDPath="/opt/verifybamid/VerifyBamID/resource/1000g.phase3.100k.b37.vcf.gz.dat.UD"
             BedPath="/opt/verifybamid/VerifyBamID/resource/1000g.phase3.100k.b37.vcf.gz.dat.bed"
             MeanPath="/opt/verifybamid/VerifyBamID/resource/1000g.phase3.100k.b37.vcf.gz.dat.mu"
      elif [[ ${reference_genome} == 'hg38' ]]
      then
             printf "VerifyBamID: Using hg38 genome\n"
             UDPath="/opt/verifybamid/VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat.UD"
             BedPath="/opt/verifybamid/VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat.bed"
             MeanPath="/opt/verifybamid/VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat.mu"
      else
          printf "ERROR: Invalid reference genome version string: %s. It should be hg37 or hg38\n" ${reference_genome}
          exit 1
      fi

      export PATH=$PATH:/VerifyBamID/bin/ && VerifyBamID \
          --UDPath ${dollar}{UDPath} \
          --BedPath ${dollar}{BedPath} \
          --MeanPath ${dollar}{MeanPath} \
          --Reference ${ref_fasta} --BamFile ${input_cram}

      # Get the contamination value from the result file
      while read line; do
          if [[ ${dollar}{line} =~ Alpha: ]]
          then
              contamination_white_space=$(echo ${dollar}{line}| cut -d':' -f 2)
              # Remove leading and training whitespace from variable
              # https://stackoverflow.com/questions/369758/how-to-trim-whitespace-from-a-bash-variable
              contamination="$(echo -e "${dollar}{contamination_white_space}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
              printf "Contamination is ${dollar}{contamination}"
              echo "${dollar}{contamination}" > contamination.txt
              break
          fi
      done <result.out

    >>>
     output {
       Array[String] DNA_contamination_output_files = [input_cram, read_string("contamination.txt")]
    }
   runtime {
      preemptible: preemptible_tries
      maxRetries: max_retries
      memory: sub(memory, "\\..*", "") + " GB"
      cpu: sub(CPUs, "\\..*", "")
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      zones: "us-central1-a us-central1-b us-east1-d us-central1-c us-central1-f us-east1-c"
      docker: docker_image
    }
  }

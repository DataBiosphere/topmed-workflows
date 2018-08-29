## Calculates DNA contamination by using the Docker image built from the Dockerfile
## https://github.com/sbg/sbg_dockstore_tools/blob/master/topmed-workflows/variant-caller/verifybamid/Dockerfile
## which is based on the VerifyBamID tool at https://github.com/griffan/VerifyBamID

workflow calulateDNAContamination {

  File? input_crai_file
  File input_cram_file

  File ref_fasta
  File ref_fasta_index

  Int? CalcContamination_CPUs
  Int CalcContamination_CPUs_default = select_first([CalcContamination_CPUs, 1])

  Int? preemptible_tries
  Int preemptible_tries_default = select_first([preemptible_tries, 3])

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  Int additional_disk = select_first([increase_disk_size, 20])

  Float reference_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  # Account for size of generated CRAM index file
  Float cram_size = size(input_cram_file, "GB")
  Float crai_size = if (defined(input_crai_file)) then size(input_crai_file, "GB") else (cram_size * 0.00003)

  String? reference_genome_version
  String reference_genome = select_first([reference_genome_version, 'hg38'])

#  String? docker_image = "images.sbgenomics.com/vladimir_obucina/topmed:VerifyBamID"
  String? docker_image = "walts-verify-bamid:latest"

  call VerifyBamID {
     input:
      input_crai = input_crai_file,
      input_cram = input_cram_file,

      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,

      reference_genome = reference_genome,
      preemptible_tries = preemptible_tries_default,

      CalcContamination_mem = CalcContamination_mem_default,
      CalcContamination_CPUs = CalcContamination_CPUs_default,
      disk_size = cram_size + crai_size + reference_size + additional_disk,
      docker_image = docker_image      
  }

  output {
      Array[File] calculate_DNA_contamination_output = VerifyBamID.DNA_contamination_output_files

  }
}

  task VerifyBamID {
     File? input_crai
     File input_cram

     File ref_fasta
     File ref_fasta_index

     String reference_genome

     Int preemptible_tries
     Int CalcContamination_CPUs
     Int CalcContamination_mem
     Float disk_size
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

    >>>
     output {
       Array[File] DNA_contamination_output_files = [input_cram, "result.out"]
    }
   runtime {
      preemptible: preemptible_tries
      memory: sub(CalcContamination_mem, "\\..*", "") + " GB"
      cpu: sub(CalcContamination_CPUs, "\\..*", "")
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      zones: "us-central1-a us-central1-b us-east1-d us-central1-c us-central1-f us-east1-c"
      docker: docker_image
    }
  }

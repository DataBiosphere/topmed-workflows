## Calculates DNA contamination by using the Docker image built from the Dockerfile
## https://github.com/sbg/sbg_dockstore_tools/blob/master/topmed-workflows/variant-caller/verifybamid/Dockerfile
## which is based on the VerifyBamID tool at https://github.com/griffan/VerifyBamID

workflow calulateDNAContamination {

  #File input_crai_file
  File input_cram_file

  File ref_fasta
  File ref_fasta_index

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  Int additional_disk = select_first([increase_disk_size, 20])

  Float reference_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  #Float cram_size = size(input_cram_file, "GB") + size(input_crai_file, "GB")
  Float cram_size = size(input_cram_file, "GB")

  String? reference_genome = 'hg38'

  String? docker_image = "images.sbgenomics.com/vladimir_obucina/topmed:VerifyBamID"

  call VerifyBamID {
     input:
      #input_crai = input_crai_file,
      input_cram = input_cram_file,

      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,

      reference_genome = reference_genome,

      disk_size = cram_size + reference_size +  + additional_disk,
      docker_image = docker_image
 
      
  }

  task VerifyBamID {
     #File input_crai
     File input_cram

     File ref_fasta
     File ref_fasta_index

     String reference_genome

     Float disk_size
     String docker_image

     # We have to use a trick to make Cromwell
     # skip substitution when using the bash ${<variable} syntax
     # This is necessary to get the <var>=$(<command>) sub shell 
     # syntax to work and assign the value to a variable when 
     # running in Cromwell
     # See https://gatkforums.broadinstitute.org/wdl/discussion/comment/44570#Comment_44570 
     String dollar = "$"

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

      echo "Running VerifyBamID"

      if [ ${dollar}{reference_genome}" == 'hg37' ]; then 
             UDPath="/VerifyBamID/resource/1000g.phase3.100k.b37.vcf.gz.dat.UD"
             BedPath="/VerifyBamID/resource/1000g.phase3.100k.b37.vcf.gz.dat.bed"
             MeanPath="/VerifyBamID/resource/1000g.phase3.100k.b37.vcf.gz.dat.mu"
      elif [ ${dollar}{reference_genome} == 'hg38']; then
             UDPath="/VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat.UD"
             BedPath="/VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat.bed"
             MeanPath="/VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat.mu"
      else
          printf "ERROR: Invalid reference genome version string: %s. It should be hg37 or hg38\n" "$reference_genome"
          exit 1
      fi
       
      export PATH=$PATH:/VerifyBamID/bin/ && VerifyBamID \
          --UDPath ${dollar}{UDPath} \
          --BedPath ${dollar}{BedPath} \
          --MeanPath ${dollar}{MeanPath} \
          --Reference ${ref_fasta} --BamFile ${input_cram}

    }
     output {
      File DNA_contamination_output_file = "result.out"
    }
   runtime {
      memory: "10 GB"
      cpu: "32"
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      zones: "us-central1-a us-central1-b us-east1-d us-central1-c us-central1-f us-east1-c"
      docker: docker_image
    }
  }

  output {
      File calculate_DNA_contamination_output = VerifyBamID.DNA_contamination_output_file
  }
}

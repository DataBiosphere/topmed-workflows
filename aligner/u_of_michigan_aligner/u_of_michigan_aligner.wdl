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

  File ref_fasta
  File ref_fasta_index
  File dbSNP_vcf
  File dbSNP_vcf_index

  # Get the size of the standard reference files as well as the additional reference files needed for BWA
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Float dbsnp_size = size(dbSNP_vcf, "GB") + size(dbSNP_vcf_index, "GB")


  call PreAlign {
     input:
      input_crai = input_crai_file,
      input_cram = input_cram_file,
      disk_size = ref_size,
      docker_image = docker_image,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }

  call Align {
     input:
      input_list_file = PreAlign.output_list_file,
      disk_size = ref_size,
      docker_image = docker_image,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }

   call PostAlign {
     input:
      num_output_files = Align.num_output_files,
      disk_size = ref_size + dbsnp_size,
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

     String pre_output_base = "pre_output_base"

     command {

      set -o pipefail
      set -e

      #echo each line of the script to stdout so we can see what is happening
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

     Float disk_size
     String docker_image

     File ref_fasta
     File ref_fasta_index


     command {

      set -o pipefail
      set -e

      #echo each line of the script to stdout so we can see what is happening
      set -o xtrace
      #to turn off echo do 'set +o xtrace'

      echo "Running alignment"
      num_output_files=0

      while read line
      do
        line_rg=$(echo $line | cut -d ' ' -f 4- | sed -e "s/ /\\\t/g")
        input_path=$(echo $line | cut -f 2 -d ' ')
        input_filename=$(basename $input_path)
        output_filename=$(basename $input_filename ".fastq.gz").cram
      
        paired_flag=""
        if [[ $input_file_name =~ interleaved\.fastq\.gz$ ]]
        then
          paired_flag="-p"
        fi
      
        bwa mem -t 32 -K 100000000 -Y $paired_flag -R "$line_rg" ${ref_fasta} $input_path | samblaster -a --addMateTags | samtools view -@ 32 -T ${ref_fasta} -C -o $output_filename -

        num_output_files=$((num_output_files+1))
 
      done <<< "$(tail -n +2 ${input_list_file})"


    }
     output {
      Int num_output_files = "$num_output_files"
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

     Int num_output_files

     command {

      set -o pipefail
      set -e

      #echo each line of the script to stdout so we can see what is happening
      set -o xtrace
      #to turn off echo do 'set +o xtrace'

      echo "Running alignment"
    
      INPUT_DIR="."
      rc=0
      for input_file in $INPUT_DIR/*.cram 
      do 
        tmp_prefix=${input_file%.cram}.tmp
        samtools sort --reference ${ref_fasta} --threads 1 -T $tmp_prefix -o ${input_file%.cram}.sorted.bam $input_file
        rc=$?
        [[ $rc != 0 ]] && break
        rm -f $input_file ${tmp_prefix}*
      done
      
      if [[ $rc == 0 ]]
      then 
        samtools merge --threads 1 -c $INPUT_DIR/merged.bam $INPUT_DIR/*.sorted.bam \
          && rm $INPUT_DIR/*.sorted.bam \
          && bam-non-primary-dedup dedup_LowMem --allReadNames --binCustom --binQualS 0:2,3:3,4:4,5:5,6:6,7:10,13:20,23:30 --log $INPUT_DIR/dedup_lowmem.metrics --recab --in $INPUT_DIR/merged.bam --out -.ubam --refFile ${ref_fasta} --dbsnp ${dbNSP_vcf} \
          | samtools view -h -C -T ${ref_fasta} -o $INPUT_DIR/output.cram --threads 1
        rc=$?
      fi
    }
     output {
      File output_cram_file = "$INPUT_DIR/output.cram"
    }
   runtime {
      memory: "10 GB"
      cpu: "32"
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      zones: "us-central1-a us-central1-b us-east1-d us-central1-c us-central1-f us-east1-c"
      docker: docker_image
    }
  }



task checkerTask {
  File inputCRAMFile
  File referenceFasta
  Int expectedNumofReads
  File docker_image
  Int total_size

  command {
    printf "The CRAM file is ${inputCRAMFile}"
    samtools flagstat ${inputCRAMFile} > cram_flagstat.txt
    ##cram_reads=$(( $(grep -i in[[:space:]]*total cram_flagstat.txt | cut -d " " -f1) - $(grep secondary cram_flagstat.txt | cut -d " " -f1) - $(grep supplementary cram_flagstat.txt | cut -d " " -f1) ))
    total_reads=$(grep -i in[[:space:]]*total cram_flagstat.txt | cut -d " " -f1)
    secondary_reads=$(grep secondary cram_flagstat.txt | cut -d " " -f1)
    supplementary_reads=$(grep supplementary cram_flagstat.txt | cut -d " " -f1)
    cram_reads=$((total_reads - secondary_reads - supplementary_reads))
    printf "The number of reads in the CRAM file is $cram_reads\n"
    exepected_num_reads=${expectedNumofReads} 
    printf "The expected number of reads is $exepected_num_reads"
 
    if [ "$cram_reads" -eq "$exepected_num_reads" ]; then
      printf "There are an equal number of reads"
      exit 0;
    else
      # the test failed
      printf "The number of reads is different"
      exit 1;
    fi
  }

  runtime {
    docker: docker_image
    disks: "local-disk " + total_size + " HDD"
  }
}

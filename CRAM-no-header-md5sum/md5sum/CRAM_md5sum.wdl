task md5 {
  File inputCRAMFile
  File inputCRAMIndexFile
  String inputFileName
  File referenceFile
  File referenceIndexFile
 
  command {
    samtools view ${inputCRAMFile} -T ${referenceFile} -t ${referenceIndexFile} | md5sum > ${inputFileName}.md5sum.txt
  }
  output {
    File value = "${inputFileName}.md5sum.txt"
 }

 runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135"
    cpu: 1
    memory: "10 GB"
    disks: "local-disk 200 HDD"
  }
}

workflow CRAM_to_md5sum {
  File inputCRAMFile
  File inputCRAMIndexFile
  File referenceFile
  File referenceIndexFile

  String inputFileName = basename("${inputCRAMFile}")
  call md5 { input: inputCRAMFile=inputCRAMFile, inputCRAMIndexFile=inputCRAMIndexFile, 
                    inputFileName=inputFileName, referenceFile=referenceFile , 
                    referenceIndexFile=referenceIndexFile}
  output {
    File md5sum_value = md5.value
  }
}

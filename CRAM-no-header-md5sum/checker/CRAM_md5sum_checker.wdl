task checkmd5 {
  File truthMd5sumFile
  File inputMd5sumFile

  command {
    diff ${inputMd5sumFile} ${truthMd5sumFile}
  }
  output {
    File value = stdout()
 }

 runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135"
    cpu: 1
    memory: "10 GB"
    disks: "local-disk 200 HDD"
  }
}

workflow CRAM_to_md5sum_checker {
  File truthMd5sumFile
  File inputMd5sumFile

  call checkmd5 { input: 
                    inputMd5sumFile=inputMd5sumFile,
                    truthMd5sumFile=truthMd5sumFile}
  output {
    File result = checkmd5.value
  }
}

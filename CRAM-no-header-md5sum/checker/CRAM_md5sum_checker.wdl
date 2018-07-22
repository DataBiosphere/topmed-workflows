task checkmd5 {
  File truthMd5sumFile
  File inputMd5sumFile
  Float disk_size

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
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
}

workflow CRAM_to_md5sum_checker {
  File truthMd5sumFile
  File inputMd5sumFile

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  Int additional_disk = select_first([increase_disk_size, 20])

  call checkmd5 { input:
                    inputMd5sumFile=inputMd5sumFile,
                    truthMd5sumFile=truthMd5sumFile,
                    disk_size = additional_disk}
  output {
    File result = checkmd5.value
  }
}

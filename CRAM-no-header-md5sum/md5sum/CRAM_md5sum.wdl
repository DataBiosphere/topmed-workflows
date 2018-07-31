task md5 {
  File inputCRAMFile
  File? inputCRAMIndexFile
  String inputFileName
  File referenceFile
  File referenceIndexFile
  Float disk_size

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
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
}

workflow CRAM_to_md5sum {
  File inputCRAMFile
  File? inputCRAMIndexFile
  File referenceFile
  File referenceIndexFile

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  Int additional_disk = select_first([increase_disk_size, 20])

  # Get the size of the standard reference file
  Float ref_size = size(referenceFile, "GB") + size(referenceIndexFile, "GB")
  Float cram_size = size(inputCRAMFile, "GB") + size(inputCRAMIndexFile, "GB")

  String inputFileName = basename("${inputCRAMFile}")
  call md5 { input: inputCRAMFile = inputCRAMFile,
                    inputCRAMIndexFile = inputCRAMIndexFile,
                    inputFileName = inputFileName, referenceFile = referenceFile,
                    referenceIndexFile = referenceIndexFile,
                    disk_size = ref_size + cram_size + additional_disk}
  output {
    File md5sum_value = md5.value
  }
}

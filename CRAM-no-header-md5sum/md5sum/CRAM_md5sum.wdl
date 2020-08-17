version 1/0

task md5 {
 input {
  File inputCRAMFile
  File? inputCRAMIndexFile
  String inputFileName
  File referenceFile
  File referenceIndexFile
  Int disk_size
 }

  command {
    samtools view ${inputCRAMFile} -T ${referenceFile} -t ${referenceIndexFile} | md5sum > ${inputFileName}.md5sum.txt
  }
  output {
    File value = "${inputFileName}.md5sum.txt"
  }

 runtime {
    docker: "quay.io/nivasquez/final-dockstore-tool-md5sum:latest"
    cpu: 1
    memory: "10 GB"
    disks: "local-disk " + disk_size + " HDD"
    continueOnReturnCode: true
  }
}

workflow CRAM_to_md5sum {
 input {
  File inputCRAMFile
  File? inputCRAMIndexFile
  File referenceFile
  File referenceIndexFile
  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size
 }

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  Int additional_disk = select_first([increase_disk_size, 20])

  # Get the size of the standard reference file
  Int ref_size = ceil(size(referenceFile, "GB") + size(referenceIndexFile, "GB"))
  Int cram_size = ceil(size(inputCRAMFile, "GB") + size(inputCRAMIndexFile, "GB"))

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

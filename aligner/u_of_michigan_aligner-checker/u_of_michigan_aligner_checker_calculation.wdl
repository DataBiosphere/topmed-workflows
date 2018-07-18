task checkerTask {
  File inputCRAMFile
  File inputTruthCRAMFile
  File referenceFile

  String docker_image

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  Int additional_disk = select_first([increase_disk_size, 20])

  Float? est_size_inputCRAMFile
  Float size_inputCRAMFile = select_first([est_size_inputCRAMFile, 30.0])

  Float? est_size_inputTruthCRAMFile
  Float size_inputTruthCRAMFile = select_first([est_size_inputTruthCRAMFile, 1.0])

  Float? est_ref_size
  Float ref_size = select_first([est_ref_size, 4.0])

  Float disk_size = size_inputTruthCRAMFile + size_inputCRAMFile + ref_size + additional_disk
  #Float disk_size = size(inputTruthCRAMFile, "GB") + size(inputCRAMFile, "GB") + size(referenceFile, "GB") + additional_disk


  command {
     # The md5sums for the SAM files without headers created from the CRAM files should match
     # if the pipeline is working properly
     # I.e. strip the headers from the CRAM files and output SAM files and compare the md5sum
     # of the two SAM files
     samtools view ${inputTruthCRAMFile} -T ${referenceFile} | md5sum > sum_file.txt &&
     samtools view ${inputCRAMFile} -T ${referenceFile} | md5sum --check sum_file.txt

  }

  runtime {
    docker: docker_image
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
}

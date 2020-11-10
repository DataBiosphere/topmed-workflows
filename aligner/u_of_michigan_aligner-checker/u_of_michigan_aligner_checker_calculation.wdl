version 1.0

task checkerTask {
  input {
    File inputCRAMFile
    File inputTruthCRAMFile
    File referenceFile

    String docker_image

    # Optional input to increase all disk sizes in case of outlier sample
    # with strange size behavior
    Int? increase_disk_size

    # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
    # Cromwell error from asking for 0 disk when the input is less than 1GB
    Int additional_disk = select_first([increase_disk_size, 200])

    Float disk_size = additional_disk
   # The size function causes an error when a relative path is provided as input in the JSON
   # input file. Somehow Cromwell confuses where the file is for the size function in this case.
   #  Float disk_size = size(inputTruthCRAMFile, "GB") + size(inputCRAMFile, "GB") + size(referenceFile, "GB") + additional_disk
  }

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
    disks: "local-disk " + ceil(disk_size) + " HDD"
  }
}

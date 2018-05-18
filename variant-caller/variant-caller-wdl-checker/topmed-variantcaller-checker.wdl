task checkerTask {
  File inputCRAMFile
  Int expectedNumofReads
  String docker_image

  Float disk_size = size(inputCRAMFile, "GB")


  command <<<
    python <<CODE
    import sys
      sys.exit(0)
    CODE
  >>>

  runtime {
    docker: docker_image
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
}

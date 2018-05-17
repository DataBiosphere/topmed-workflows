task checkerTask {
  File inputCRAMFile
  Int expectedNumofReads
  String docker_image

  Float disk_size = size(inputCRAMFile, "GB")

  command <<<
    python <<CODE
    with open(filepath1, 'r') as default_file:
       good_data = []
       for line in default_file:
           line = line.strip()
           if not line.startswith('#'):
               good_data.append(line.split('\t'))

    with open(filepath2, 'r') as test_file:
       test_data = []
       for line in test_file:
           line = line.strip()
           if not line.startswith('#'):
               test_data.append(line.split('\t'))

    for i in range(len(test_data)):
       if test_data[i] != good_data[i]:
           for j in range(len(test_data[i])):
               # Only compare chromosome, position, ID, reference, and alts.
               # Quality score may vary (<1%) between systems because of
               # (assumed) rounding differences.  Same for the "info" sect.
               if j < 5:
                   assert test_data[i][j] == good_data[i][j], "File does not match: %r" % filepath1
    CODE
  >>>

  runtime {
    docker: docker_image
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
}

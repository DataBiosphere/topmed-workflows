task checkerTask {
  File inputTruthVCFFile
  File inputTestVCFFile
  String docker_image

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  Int additional_disk = select_first([increase_disk_size, 20])

  Float disk_size = size(inputTruthVCFFile, "GB") + size(inputTestVCFFile, "GB") + additional_disk


  command <<<
    python3 <<CODE
    from __future__ import print_function, division
    import sys, os, tarfile, gzip
    
    def read_and_compare_vcfs_from_tar_gz(tar_gz_test, tar_gz_truth):
        """
        Reads the VCF files from the tar gz file produced by the U of Michigan
        WDL variant caller and the truth targ gz file and compares each of them
        and returns 1 if any do not compare favorably and 0 if all of them do.
    
        """
    
        print("test file:{}    truth file:{}".format(tar_gz_test, tar_gz_truth))
        with tarfile.open(tar_gz_test, "r") as test_variant_caller_output, \
             tarfile.open(tar_gz_truth, "r") as truth_variant_caller_output:
    
             test_vcf_file_names = test_variant_caller_output.getnames()
             #print("vcf file names are:{}".format(test_vcf_file_names))

             # Check that the truth VCF tar file is not empty; if it is something is wrong
             truth_vcf_file_names = truth_variant_caller_output.getnames()
             if not truth_vcf_file_names or len(truth_vcf_file_names) == 0:
                 print("The truth tar gz file is empty", file=sys.stderr)
                 sys.exit(1)

 
             for truth_vcf_file_info in truth_variant_caller_output.getmembers():
                 #truth_vcf_file_name = os.path.basename(truth_vcf_file_info.name)
                 #print("Truth vcf file name is:{}".format(truth_vcf_file_name))           
    
                 if truth_vcf_file_info.isfile() and \
                    os.path.basename(truth_vcf_file_info.name).startswith("chr") and \
                    os.path.basename(truth_vcf_file_info.name).endswith("vcf.gz"):
    
                    print("Checking to see if truth vcf file {} is present in {}".format(truth_vcf_file_info.name, tar_gz_test))
                    # If a VCF file is missing in the test output then
                    # the VCFs are not the same and return error
                    if truth_vcf_file_info.name not in test_vcf_file_names:
                        print("VCF file {} is missing from variant caller output".format(truth_vcf_file_info.name), file=sys.stderr)
                        sys.exit(1)
   
                    # Get file like objects for the gzipped vcf files
                    test_vcf_file_info = test_variant_caller_output.getmember(truth_vcf_file_info.name) 
                    test_vcf_file = test_variant_caller_output.extractfile(test_vcf_file_info)
                    truth_vcf_file = truth_variant_caller_output.extractfile(truth_vcf_file_info)

                    vcfs_are_not_concordant = compare_vcf_files(truth_vcf_file, 
                           truth_vcf_file_info.name, test_vcf_file, test_vcf_file_info.name)
                    if vcfs_are_not_concordant:
                        print("VCF file {} is not concordant with the truth set VCF file".format(test_vcf_file.name), file=sys.stderr)
                        sys.exit(1)
        print("The VCF files are concordant")
        sys.exit(0)

    read_and_compare_vcfs_from_tar_gz("${inputTestVCFFile}", "${inputTruthVCFFile}")

    CODE
  >>>

  runtime {
    docker: docker_image
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
}

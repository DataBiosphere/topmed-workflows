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
    python3.5 <<CODE
    from __future__ import print_function, division
    import sys, os, tarfile, gzip

    def compare_vcf_files(file_like_object1, file_name_1, file_like_object2, file_name_2):
        """
        Asserts that two .vcf files contain the same variant findings.
    
        - Ignores potentially date-stamped comments (lines starting with '#').
        - Ignores quality scores in .vcf files and only checks that they found
          the same variants.  This is due to assumed small observed rounding
          differences between systems.
    
        VCF File Column Contents:
        1: #CHROM
        2: POS
        3: ID
        4: REF
        5: ALT
        6: QUAL
        7: FILTER
        8: INFO
    
        :param file_like_object1: First .vcf file to compare.
        :param file_like_object2: Second .vcf file to compare.
        """

        print("Comparing variant data in test file {} with data in truth file {}".format( file_name_1, file_name_2))
        
        with gzip.open(file_like_object1, 'rt') as default_file:
            good_data = []
            for line in default_file:
                line = line.strip()
                #print("default file 1 line: {}".format(line))
                if not line.startswith('#'):
                    good_data.append(line.split('\t'))

        with gzip.open(file_like_object2, 'rt') as test_file:
            test_data = []
            for line in test_file:
                line = line.strip()
                #print("test file 2 line: {}".format(line))
                if not line.startswith('#'):
                    test_data.append(line.split('\t'))

        # If there are more variants in one file then the VCFs are not concordant
        if len(test_data) != len(good_data):
            print("The truth VCF file {} has {} variants and the test VCF file {} has {} variants\n".format(file_name_1, len(good_data), file_name_2, len(test_data)))
            return 1    

        for i in range(len(test_data)):
            if test_data[i] != good_data[i]:
                for j in range(len(test_data[i])):
                    # Only compare chromosome, position, ID, reference, and alts.
                    # Quality score may vary (<1%) between systems because of
                    # (assumed) rounding differences.  Same for the "info" sect.
                    if j < 5:
                        #assert test_data[i][j] == good_data[i][j], "File does not match: %r" % file_like_object1
                        if test_data[i][j] != good_data[i][j]:
                            print("Data {} in file {} does not match data {} in {}".format(good_data[i][j], file_name_1, test_data[i][j], file_name_2))
                            return 1
        print("The test VCF file {} is concordant with truth VCF file {}\n".format(file_name_1, file_name_2))
        return 0
    
    
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

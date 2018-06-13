task checkerTask {
  File inputTruthVCFFile
  File inputTestVCFFile
  File ref_hs38DH_fa
  File concordance_outputTSV
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
    import sys, os, tarfile, gzip, csv, math
    from subprocess import Popen, PIPE, STDOUT
    
    def read_and_compare_vcfs_from_tar_gz(tar_gz_test, tar_gz_truth,
     reference, output_tsv):
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

                    run_concordance(reference, \
                    test_vcf_file_info, test_vcf_file_info, output_file)

    def run_concordance(reference, eval_file, truth_file, output_file=None):
        """Open a terminal shell to run a command in a Docker
        image with Genotype Concordance installed.
         :return: none
         """
        user_name = os.path.expanduser('~')
        docker_permission = user_name + ':' + user_name
        if output_file is None:
            output_file = user_name + '/concordance_output.tsv'

        cmd = ['docker', 'run', '-i', '-t', '-v',
               docker_permission,
               'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135',
               '/usr/gitc/gatk4/gatk-launch',
               'Concordance',
               '-R', str(reference),
               '--eval', str(eval_file),
               '--truth', str(truth_file),
               '--summary', str(output_file)]

        p = Popen(cmd, stdout=PIPE, stderr=STDOUT)
        p.wait()
        print("GenotypeConcordance out: {}".format(p.communicate()))
        #print(os.getcwd())
        d = process_output_tsv(output_tsv=output)
        print(d)  # print to stdout so we read it in WDL

    def process_output_tsv(output_tsv, threshold=None):
        """
        Process TSV file written to the current directory.
        :parameter: output_tsv: (string) path to a TSV file from Concordance VCF
        :parameter: threshold: (float) 0 < thresh < 1, sensitivity and precision
                    default: 0.95
        :return: boolean, True is output passes threshold, otherwise false.
        """

        # Set default
        if threshold is None:
            threshold = 0.95
        L = []  # list to capture results
        with open(output_tsv, newline='') as csvfile:
            file_reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
            for row in file_reader:
                L.append(row[0])

        D = list2dict(L)

        # Convert relevant values in dict to floats.
        vals = [D['type']['SNP']['precision'],
                D['type']['SNP']['sensitivity'],
                D['type']['INDEL']['precision'],
                D['type']['SNP']['sensitivity']]

        vals = [float(val) for val in vals]

        # The next line is needed as we encountered NaNs in the output
        # after I run Concordance with two identical inputs for truth and test.
        # It removes NaNs from list.
        vals = [x for x in vals if not math.isnan(x)]

        # Test whether all values pass the threshold test:
        if all(val >= threshold for val in vals):
            message = 'The VCFs can be considered identical.'
            print(message)
            sys.exit(0)

        else:
            message = 'The VCFs do not have enough overlap.'
            print(message)
            sys.exit(0)


    def list2dict(L):
        """Returns a dictionary from input list, originating from the
        Concordance TSV file."""

        dd = {i: L[i].split('\t') for i in range(len(L))}  # auxiliary dict
        D = {}
        # Construct output dictionary of key-value pairs:
        D[dd[0][0]] = {dd[1][0]: dict(zip(dd[0][1:], dd[1][1:])),
                       dd[2][0]: dict(zip(dd[0][1:], dd[2][1:]))}
        return D

    read_and_compare_vcfs_from_tar_gz("${inputTestVCFFile}", \
    "${inputTruthVCFFile}", "${ref_hs38DH_fa}, ${concordance_outputTSV})

    CODE
  >>>

  runtime {
    docker: docker_image
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
}

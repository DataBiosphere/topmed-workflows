task checkerTask {
  File inputTruthVCFFile
  File inputTestVCFFile

  File ref_1000G_omni2_5_b38_sites_PASS_vcf_gz
  File ref_1000G_omni2_5_b38_sites_PASS_vcf_gz_tbi
  File chr10_vcf
  File chr11_KI270927v1_alt_vcf
  File chr11_vcf
  File chr12_vcf
  File chr13_vcf
  File chr14_GL000009v2_random_vcf
  File chr14_KI270846v1_alt_vcf
  File chr14_vcf
  File chr15_vcf
  File chr16_vcf
  File chr17_KI270857v1_alt_vcf
  File chr17_KI270862v1_alt_vcf
  File chr17_KI270909v1_alt_vcf
  File chr17_vcf
  File chr18_vcf
  File chr19_KI270938v1_alt_vcf
  File chr19_vcf
  File chr1_KI270706v1_random_vcf
  File chr1_KI270766v1_alt_vcf
  File chr1_vcf
  File chr20_vcf
  File chr21_vcf
  File chr22_KI270879v1_alt_vcf
  File chr22_KI270928v1_alt_vcf
  File chr22_vcf
  File chr2_KI270773v1_alt_vcf
  File chr2_KI270894v1_alt_vcf
  File chr2_vcf
  File chr3_vcf
  File chr4_GL000008v2_random_vcf
  File chr4_vcf
  File chr5_vcf
  File chr6_vcf
  File chr7_KI270803v1_alt_vcf
  File chr7_vcf
  File chr8_KI270821v1_alt_vcf
  File chr8_vcf
  File chr9_vcf
  File chrUn_KI270742v1_vcf
  File chrX_vcf
  File ref_dbsnp_142_b38_vcf_gz
  File ref_dbsnp_142_b38_vcf_gz_tbi
  File ref_dbsnp_All_vcf_gz
  File ref_dbsnp_All_vcf_gz_tbi
  File ref_hapmap_3_3_b38_sites_vcf_gz
  File ref_hapmap_3_3_b38_sites_vcf_gz_tbi
  File ref_hs38DH_bs_umfa
  File ref_hs38DH_dict
  File ref_hs38DH_fa
  File ref_hs38DH_fa_alt
  File ref_hs38DH_fa_amb
  File ref_hs38DH_fa_ann
  File ref_hs38DH_fa_bwt
  File ref_hs38DH_fa_fai
  File ref_hs38DH_fa_pac
  File ref_hs38DH_fa_sa
  File ref_hs38DH_winsize100_gc

  String docker_concordance_image

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  Int additional_disk = select_first([increase_disk_size, 500])

  # The size function causes an error when a relative path is provided as input in the JSON
  # input file, which is done in the test input JSON 
  # topmed_freeze3_calling_checker.wdl.local.json, so bump up the additional_disk size.
  # Somehow Cromwell confuses where the file is for the size function in this case.
  #Float reference_size = (
  #size(ref_1000G_omni2_5_b38_sites_PASS_vcf_gz, "GB") +
  #size(ref_1000G_omni2_5_b38_sites_PASS_vcf_gz_tbi, "GB") +
  #size(chr10_vcf, "GB") +
  #size(chr11_KI270927v1_alt_vcf, "GB") +
  #size(chr11_vcf, "GB") +
  #size(chr12_vcf, "GB") +
  #size(chr13_vcf, "GB") +
  #size(chr14_GL000009v2_random_vcf, "GB") +
  #size(chr14_KI270846v1_alt_vcf, "GB") +
  #size(chr14_vcf, "GB") +
  #size(chr15_vcf, "GB") +
  #size(chr16_vcf, "GB") +
  #size(chr17_KI270857v1_alt_vcf, "GB") +
  #size(chr17_KI270862v1_alt_vcf, "GB") +
  #size(chr17_KI270909v1_alt_vcf, "GB") +
  #size(chr17_vcf, "GB") +
  #size(chr18_vcf, "GB") +
  #size(chr19_KI270938v1_alt_vcf, "GB") +
  #size(chr19_vcf, "GB") +
  #size(chr1_KI270706v1_random_vcf, "GB") +
  #size(chr1_KI270766v1_alt_vcf, "GB") +
  #size(chr1_vcf, "GB") +
  #size(chr20_vcf, "GB") +
  #size(chr21_vcf, "GB") +
  #size(chr22_KI270879v1_alt_vcf, "GB") +
  #size(chr22_KI270928v1_alt_vcf, "GB") +
  #size(chr22_vcf, "GB") +
  #size(chr2_KI270773v1_alt_vcf, "GB") +
  #size(chr2_KI270894v1_alt_vcf, "GB") +
  #size(chr2_vcf, "GB") +
  #size(chr3_vcf, "GB") +
  #size(chr4_GL000008v2_random_vcf, "GB") +
  #size(chr4_vcf, "GB") +
  #size(chr5_vcf, "GB") +
  #size(chr6_vcf, "GB") +
  #size(chr7_KI270803v1_alt_vcf, "GB") +
  #size(chr7_vcf, "GB") +
  #size(chr8_KI270821v1_alt_vcf, "GB") +
  #size(chr8_vcf, "GB") +
  #size(chr9_vcf, "GB") +
  #size(chrUn_KI270742v1_vcf, "GB") +
  #size(chrX_vcf, "GB") +
  #size(ref_dbsnp_142_b38_vcf_gz, "GB") +
  #size(ref_dbsnp_142_b38_vcf_gz_tbi, "GB") +
  #size(ref_dbsnp_All_vcf_gz, "GB") +
  #size(ref_dbsnp_All_vcf_gz_tbi, "GB") +
  #size(ref_hapmap_3_3_b38_sites_vcf_gz, "GB") +
  #size(ref_hapmap_3_3_b38_sites_vcf_gz_tbi, "GB") +
  #size(ref_hs38DH_bs_umfa, "GB") +
  #size(ref_hs38DH_dict, "GB") +
  #size(ref_hs38DH_fa, "GB") +
  #size(ref_hs38DH_fa_alt, "GB") +
  #size(ref_hs38DH_fa_amb, "GB") +
  #size(ref_hs38DH_fa_ann, "GB") +
  #size(ref_hs38DH_fa_bwt, "GB") +
  #size(ref_hs38DH_fa_fai, "GB") +
  #size(ref_hs38DH_fa_pac, "GB") +
  #size(ref_hs38DH_fa_sa, "GB") +
  #size(ref_hs38DH_winsize100_gc, "GB")
  #)

  #Float disk_size = reference_size + size(inputTruthVCFFile, "GB") + size(inputTestVCFFile, "GB") + additional_disk
  Float disk_size = additional_disk


  command <<<
    python3 <<CODE
    from __future__ import print_function, division
    import sys, os, tarfile, gzip, csv, math, shutil
    from subprocess import Popen, PIPE, STDOUT
    
    def read_and_compare_vcfs_from_tar_gz(tar_gz_truth, tar_gz_test, reference):
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
             test_vcf_file_basenames = [os.path.basename(file_name) for file_name in test_vcf_file_names]
             #print("vcf file basenames are:{}".format(test_vcf_file_basenames))

             # Check that the truth VCF tar file is not empty; if it is something is wrong
             truth_vcf_file_names = truth_variant_caller_output.getnames()
             if not truth_vcf_file_names or len(truth_vcf_file_names) == 0:
                 print("The truth tar gz file is empty", file=sys.stderr)
                 sys.exit(1)


             for truth_vcf_file_info in truth_variant_caller_output.getmembers():
                 truth_vcf_file_name = os.path.basename(truth_vcf_file_info.name)
                 #print("Truth vcf file name is:{}".format(truth_vcf_file_name))

                 if truth_vcf_file_info.isfile() and \
                    os.path.basename(truth_vcf_file_info.name).startswith("chr") and \
                    os.path.basename(truth_vcf_file_info.name).endswith("vcf.gz"):

                    #print("Checking to see if truth vcf file {} is present in {}".format(os.path.basename(truth_vcf_file_info.name), tar_gz_test))
                    # If a VCF file is missing in the test output then
                    # the VCFs are not the same and return error
                    if os.path.basename(truth_vcf_file_info.name) not in test_vcf_file_basenames:
                        print("VCF file {} is missing from variant caller output".format(os.path.basename(truth_vcf_file_info.name)), file=sys.stderr)
                        sys.exit(1)

                    # Get file like objects for the gzipped vcf files
                    for test_vcf_file_info in test_variant_caller_output.getmembers():
                       if os.path.basename(test_vcf_file_info.name) == os.path.basename(truth_vcf_file_info.name):
                           test_vcf_file = test_variant_caller_output.extractfile(test_vcf_file_info)
                           #print("Got test vcf file:{} with file name {}".format(test_vcf_file, test_vcf_file_info.name)) 
                           break;

                    truth_vcf_file = truth_variant_caller_output.extractfile(truth_vcf_file_info)

                    # The following code writes the file-like objects to the host disk.
                    # This is necessary for the Java executable to read the files.
                    fnames = ['truth.vcf', 'test.vcf']
                    file_like_objects = [truth_vcf_file, test_vcf_file]
                    cnt = 0
                    for fl_obj in file_like_objects:
                        with gzip.open(fl_obj, 'r') as f_in, open(fnames[cnt], 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                        cnt = cnt + 1

                    # Run the GATK VCF checker procedure with those input files.
                    run_concordance(reference, 'truth.vcf', 'test.vcf')

    def run_concordance(reference, truth_file, eval_file):
        """Open a terminal shell to run a command in a Docker
        image with Genotype Concordance installed.
         :return: none
        """

        # Create file to capture GATK Concordance on local host.
        output_file = 'concordance_outputTSV.tsv'

        cmd = ['/usr/gitc/gatk4/gatk-launch',
               'Concordance',
               '-R', str(reference),
               '--eval', eval_file,
               '--truth', truth_file,
               '--summary', str(output_file)]

        p = Popen(cmd, stdout=PIPE, stderr=STDOUT)

        # Show the output from inside the Docker on the host terminal.
        print("GenotypeConcordance out: {}".format(p.communicate()))

        d = process_output_tsv(output_tsv=output_file)
        print(d)  # print to stdout so we read it in WDL

        os.remove(outfile)

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

        try:
            with open(output_tsv, newline='') as csvfile:
                file_reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for row in file_reader:
                    L.append(row[0])
        except FileNotFoundError:
            print('no output TSV file found')

        D = list2dict(L)

        # Convert relevant values in dict to floats.
        vals = [D['type']['SNP']['precision'],
                D['type']['SNP']['sensitivity'],
                D['type']['INDEL']['precision'],
                D['type']['INDEL']['sensitivity']]

        if not vals:
            msg = 'GATK Concordance VCF checker output is empty - aborting.'
            print(msg)
            sys.exit(1)

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
            sys.exit(1)


    def list2dict(L):
        """Returns a dictionary from input list, originating from the
        Concordance TSV file."""

        dd = {i: L[i].split('\t') for i in range(len(L))}  # auxiliary dict
        D = {}
        # Construct output dictionary of key-value pairs:
        D[dd[0][0]] = {dd[1][0]: dict(zip(dd[0][1:], dd[1][1:])),
                       dd[2][0]: dict(zip(dd[0][1:], dd[2][1:]))}
        return D

    read_and_compare_vcfs_from_tar_gz("${inputTruthVCFFile}", \
    "${inputTestVCFFile}", "${ref_hs38DH_fa}")

    CODE
  >>>


  runtime {
    docker: docker_concordance_image
    cpu: "16"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
}

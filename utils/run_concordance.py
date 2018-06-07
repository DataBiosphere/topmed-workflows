#!/usr/bin/env python

import os
import tarfile
import sys
from subprocess import Popen, PIPE, STDOUT

def main(test_fn, truth_fn, reference, output):

    with tarfile.open(test_fn, 'r') as test_vcf, \
        tarfile.open(truth_fn, 'r') as truth_vcf:

        test_vcf_fnames = test_vcf.getnames()  # is a list
        #test_vcf_fnames = test_vcf.getmembers()  # is a list

        # Make sure truth VCF tar file is not empty; else something is wrong.
        truth_vcf_file_names = truth_vcf.getnames()
        if not truth_vcf_file_names or len(truth_vcf_file_names) == 0:
            print("The truth tar gz file is empty")
            #print("The truth tar gz file is empty", file=sys.stderr)
            sys.exit(1)

        for truth_vcf_info in truth_vcf.getmembers():
            #test_vcf_fname = os.path.basename(test_vcf_info.name)
            #print("Test VCF filename: {}".format(test_vcf_fname))

            if truth_vcf_info.isfile() and \
                    os.path.basename(
                        truth_vcf_info.name).startswith("chr") and \
                    os.path.basename(
                        truth_vcf_info.name).endswith("vcf.gz"):

                print("Checking to see if truth vcf file {} is present in {}".format(
                        truth_vcf_info.name, test_fn))
                # If a VCF file is missing in the test output then
                # the VCFs are not the same and return error
                if truth_vcf_info.name not in test_vcf_fnames:
                    # (checks whether string is in list)
                    print("VCF file {} is missing from test variant caller output".
                          format(truth_vcf_info.name))
                    #print("VCF file {} is missing from variant caller output".format(test_vcf_info.name), file=sys.stderr)
                    sys.exit(1)

                # Get file like objects for the gzipped vcf files
                test_vcf_file_info = test_vcf.getmember(truth_vcf_info.name)
                test_vcf_file = test_vcf.extractfile(test_vcf_file_info)

                #truth_vcf_file_info = truth_vcf.getmember(truth_vcf_info.name)
                truth_vcf_file = truth_vcf.extractfile(truth_vcf_info)
                print("Test VCF: {}".format(test_vcf_file))

                #p = Popen(['java', '-cp',
                p = Popen(['java', '-jar', '$HOME/bin/GenomeAnalysisTK.jar',
                           '-T', 'GenotypeConcordance',
                           '-R', str(reference),
                           '-eval', str(test_vcf_file),
                           '-comp', str(truth_vcf_file),
                           '-o', str(output)],
                          stdout=PIPE, stderr=STDOUT)



                print("p output arg from Java command: {}".format(p.communicate()))


test_fn = '/home/ubuntu/vcf_test/test.varcall.tar.gz'
truth_fn = '/home/ubuntu/vcf_test/truth.varcall.tar.gz'
reference = '/home/ubuntu/hg38/hs38DH.fa'
output = '/home/ubuntu/vcf_test/out.grp'

main(test_fn, truth_fn, reference, output)

if __name__=='__main__':
    main(test_fn, truth_fn, reference, output)

#!/usr/bin/env python

import os
import tarfile
import sys
from subprocess import Popen, PIPE, STDOUT

def run_concordance(test_fn, truth_fn, reference, output):

    with tarfile.open(test_fn, 'r') as test_vcf, \
        tarfile.open(truth_fn, 'r') as truth_vcf:

        test_vcf_fnames = test_vcf.getmembers()
        truth_vcf_fnames = truth_vcf.getmembers()

        for test_vcf_info in test_vcf.getmembers():
            test_vcf_fname = os.path.basename(test_vcf_info.name)
            #print("Test VCF filename: {}".format(test_vcf_fname))

            if test_vcf_info.isfile() and \
                    os.path.basename(
                        test_vcf_info.name).startswith("chr") and \
                    os.path.basename(
                        test_vcf_info.name).endswith("vcf.gz"):

                print("Checking to see if truth vcf file {} is present in {}".format(
                        test_vcf_info.name, test_fn))
                # If a VCF file is missing in the test output then
                # the VCFs are not the same and return error
                if test_vcf_info.name not in test_vcf_fnames:
                    print(
                        "VCF file {} is missing from variant caller output".format(
                            test_vcf_info.name), file=sys.stderr)
                    sys.exit(1)

                # Get file like objects for the gzipped vcf files
                test_vcf_file_info = test_vcf.getmember(test_vcf_info.name)
                test_vcf_file = test_vcf.extractfile(test_vcf_file_info)

                truth_vcf_file_info = truth_vcf.getmember(test_vcf_info.name)
                truth_vcf_file = truth_vcf.extractfile(truth_vcf_file_info)
                print("Test VCF: {}".format(test_vcf_file))

                p = Popen(['java', '-cp',
                           'java -jar /home/michael/bin/GenomeAnalysisTK.jar',
                           '-T', 'GenotypeConcordance',
                           '-R', reference,
                           '-eval', test_vcf_file,
                           '-comp', truth_vcf_file,
                           '-o', output],
                          stdout=PIPE, stderr=STDOUT)

                print("p output arg from Java command: {}".format(p))


test_fn = '/home/ubuntu/vcf_test/test.vcf'
truth_fn = '/home/ubuntu/vcf_test/truth.vcf'
reference = '/home/ubuntu/hg38/hs38DH.fa'
output = '/home/ubuntu/vcf_test/out.grp'
run_concordance(test_fn, truth_fn, reference, output)

#if __name__=='__main__':
#    run_concordance()

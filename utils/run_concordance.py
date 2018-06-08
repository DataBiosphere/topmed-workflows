#!/usr/bin/env python

import os
import tarfile
import sys
import tempfile
import shutil
from subprocess import Popen, PIPE, STDOUT

uname = home = os.path.expanduser('~')

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

        for truth_vcf_member in truth_vcf.getmembers():
            #test_vcf_fname = os.path.basename(test_vcf_info.name)
            #print("Test VCF filename: {}".format(test_vcf_fname))

            if truth_vcf_member.isfile() and \
                    os.path.basename(
                        truth_vcf_member.name).startswith("chr") and \
                    os.path.basename(
                        truth_vcf_member.name).endswith("vcf.gz"):

                print("Checking to see if truth vcf file {} is present in {}".format(
                        truth_vcf_member.name, test_fn))
                # If a VCF file is missing in the test output then
                # the VCFs are not the same and return error
                if truth_vcf_member.name not in test_vcf_fnames:
                    # (checks whether string is in list)
                    print("VCF file {} is missing from test variant caller output".
                          format(truth_vcf_member.name))
                    #print("VCF file {} is missing from variant caller output".format(test_vcf_info.name), file=sys.stderr)
                    sys.exit(1)

                with tempfile.TemporaryDirectory() as tmpdirname:
                    os.chdir(tmpdirname)
                    truth_vcf_file = truth_vcf.getmember(truth_vcf_member.name)
                    truth_vcf.extract(truth_vcf_file, path=tmpdirname)
                    # Need to rename directory as we need to extract another one.
                    os.rename(os.listdir(tmpdirname)[0], 'truth')
                    truth_path = os.path.join(tmpdirname,
                                              os.listdir(tmpdirname)[0])
                    truth_vcf_name = get_vcf_filename(truth_path)
                    print(truth_vcf_name)

                    truth_vcf.extract(truth_vcf_file, path=tmpdirname)
                    # Need to rename directory as we need to extract another one.
                    os.rename(os.listdir(tmpdirname)[0], 'test')
                    test_path = os.path.join(tmpdirname,
                                              os.listdir(tmpdirname)[0])
                    test_vcf_name = get_vcf_filename(test_path)
                    print(test_vcf_name)

                    jar = uname + '/bin/GenomeAnalysisTK.jar'

                    p = Popen(['java', '-jar', str(jar),
                               '-T', 'GenotypeConcordance',
                               '-R', str(reference),
                               '-eval', str(truth_vcf_name),
                               '-comp', str(test_vcf_name),
                               '-o', str(output)],
                              stdout=PIPE, stderr=STDOUT)

                    print("Java output: {}".format(p.communicate()))

def get_vcf_filename(dirname):
    for root, dirs, files in os.walk(dirname):
        for name in files:
            return os.path.join(root, name)


test_fn = uname + '/dev/topmed-workflows-local/vcf_files/test.varcall.tar.gz'
truth_fn = uname + '/dev/topmed-workflows-local/vcf_files/truth.varcall.tar.gz'
reference = uname + '/dev/hg38/hs38DH.fa'
output = uname + '/dev/topmed-workflows-local/vcf_test/out.grp'

main(test_fn, truth_fn, reference, output)

if __name__=='__main__':
    main(test_fn, truth_fn, reference, output)














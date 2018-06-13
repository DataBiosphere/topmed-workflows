#!/usr/bin/env python

import os
import csv
import math
import tarfile
import sys
import tempfile
import shutil
import gzip
from subprocess import Popen, PIPE, STDOUT


#def main(test_fn, truth_fn, reference, output):
def main():
    """Open a terminal shell to run a command in a Docker 
    image with Genotype Concordance installed.
     :return: none 
     """
    user_name = os.path.expanduser('~')
    output = user_name + '/dev/topmed-workflows/utils/concordance_output.tsv'

    cmd = ['docker', 'run', '-i', '-t', '-v',
           '/home/michael:/home/michael',
           'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135',
           '/usr/gitc/gatk4/gatk-launch',
           'Concordance',
           '-R', '/home/michael/dev/hg38/hs38DH.fa',
           '--eval', '/home/michael/dev/topmed-workflows-local/vcf_files/root/topmed_freeze3_calling/out/paste/chr17_1_83257441_paste.sites.vcf.gz',
           '--truth', '/home/michael/dev/topmed-workflows-local/vcf_files/root/topmed_freeze3_calling/out/paste/chr17_1_83257441_paste.sites.vcf.gz',
           '--summary', str(output)]

    p = Popen(cmd, stdout=PIPE, stderr=STDOUT)
    print("GenotypeConcordance out: {}".format(p.communicate()))
    #print(os.getcwd())
    d = process_output_tsv(output_tsv=output)
    print(d)

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

    # Convert all values in dict to floats.
    snp_prec = float(D['type']['SNP']['precision'])
    snp_sens = float(D['type']['SNP']['sensitivity'])
    ind_prec = float(D['type']['INDEL']['precision'])
    ind_sens = float(D['type']['SNP']['sensitivity'])

    vals = [snp_prec, snp_sens, ind_prec, ind_sens]
    # The next line is needed as I (MK) encounted NaNs in the output
    # after I run Concordance with two identical inputs for truth and test.
    vals = [x for x in vals if not math.isnan(x)]

    # Test whether all values pass the threshold test:
    if all(val >= 0.95 for val in vals):
        message = 'The VCFs can be considered identical.'
        print(message)
        return 0
    else:
        message = 'The VCFs do not have enough overlap.'
        print(message)
        return 1


def list2dict(L):
    """Returns a dictionary from input list, originating from the
    Concordance TSV file."""

    dd = {i: L[i].split('\t') for i in range(len(L))}  # auxiliary dict
    D = {}
    # Construct output dictionary of key-value pairs:
    D[dd[0][0]] = {dd[1][0]: dict(zip(dd[0][1:], dd[1][1:])),
                   dd[2][0]: dict(zip(dd[0][1:], dd[2][1:]))}
    return D


if __name__=='__main__':
    main()














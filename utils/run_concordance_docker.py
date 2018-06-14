#!/usr/bin/env python

import os
import csv
import math
import tempfile
from subprocess import Popen, PIPE, STDOUT


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
            D['type']['INDEL']['sensitivity']]

    vals = [float(val) for val in vals]

    # The next line is needed as we encountered NaNs in the output
    # after we ran Concordance with two identical inputs for truth and test.
    # The following lines removes NaNs from the list.
    vals = [x for x in vals if not math.isnan(x)]

    # Test whether all values pass the threshold test:
    if all(val >= threshold for val in vals):
        message = 'The VCFs can be considered identical.'
        print(message)
        return 0  # this line is sys.exit(0) in the WDL

    else:
        message = 'The VCFs do not have enough overlap.'
        print(message)
        return 1  # this line is sys.exit(0) in WDL


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
    user_name = os.path.expanduser('~')
    thisdir = os.getcwd()
    os.chdir(thisdir)
    newdir = os.path.split(thisdir)
    print(os.getcwd())
    ref = user_name + '/dev/hg38/hs38DH.fa'
    eval_fn = '../../test_data/chr17_1_83257441_paste.sites.vcf.gz'
    truth = eval_fn
    run_concordance(ref, eval_fn, truth, output_file=outfile)


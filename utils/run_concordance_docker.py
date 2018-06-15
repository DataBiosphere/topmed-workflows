#!/usr/bin/env python

import os
import csv
import math
from pprint import pprint
from subprocess import Popen, PIPE, STDOUT


def run_concordance(reference, eval_file, truth_file, output_file):
    """Run Genotype Concordance installed in a Docker image.
    
     :return: none 
     """
    user_name = os.path.expanduser('~')
    # Give Docker access to local system.
    docker_permission = user_name + ':' + user_name
    docker_url = 'us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135'

    cmd = ['docker', 'run', '-i', '-v',
           str(docker_permission),
           str(docker_url),
           '/usr/gitc/gatk4/gatk-launch',
           'Concordance',
           '-R', str(reference),
           '--eval', str(eval_file),
           '--truth', str(truth_file),
           '--summary', str(output_file)]

    p = Popen(cmd, stdout=PIPE, stderr=STDOUT)
    p.wait()
    print("GenotypeConcordance out: {}".format(p.communicate()))


def process_output_tsv(output_tsv, threshold=None, print_dict=False):
    """Process TSV file written to the current directory.
    
    :parameter: output_tsv: (string) path to a TSV file from Concordance VCF
    :parameter: threshold: (float) 0 < thresh < 1, sensitivity and precision
                default: 0.95
    :return: 0 (if VCFs overlap within 95%), otherwise 1.
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
    if print_dict:
        pprint('Concordance output: {}'.format(D))

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
    Concordance TSV file.
    
    :parameter L: a list from a Concordance summary output file.
    :return D: a dict with key / value pairs from input list
    """

    dd = {i: L[i].split('\t') for i in range(len(L))}  # auxiliary dict
    D = {}
    # Construct output dictionary of key-value pairs:
    D[dd[0][0]] = {dd[1][0]: dict(zip(dd[0][1:], dd[1][1:])),
                   dd[2][0]: dict(zip(dd[0][1:], dd[2][1:]))}
    return D


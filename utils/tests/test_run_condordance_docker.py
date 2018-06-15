#!/usr/bin/env python

import unittest
from pprint import pprint
import csv
import os
from utils.run_concordance_docker import process_output_tsv, \
    list2dict, run_concordance


class TestRunConcordanceDocker(unittest.TestCase):
    """Test the functionality of the Python wrapper for Broad's GATK
    GenomeAnalysis to perform a concordance test on a truth compared 
    to a test VCF from the Variant Calling workflow.
    
    NOTE:
        Most tests depend on the hg38 reference genome. For the tests to run
        hg38 needs to be residing on the local system, and you need to 
        specify the relative path for file hs38DH.fa in the variable 
        "rel_path_to_hg38_ref" in the set up. But Concordance needs the 
         complete hg38 directory. The code expects to find it in 
        "/home/$USER". The Unix user name is added automatically.
    """

    def setUp(self):
        current_path = os.getcwd()
        # Move up two levels.
        test_data_path = os.path.split(os.path.split(current_path)[0])[0]
        user_name = os.path.expanduser('~')

        rel_path_to_hg38_ref = 'dev/hg38/hs38DH.fa'

        self.reference = os.path.join(os.sep, user_name, rel_path_to_hg38_ref)
        self.tsv_file = os.path.join(
            os.sep, test_data_path,
            'test_data/concordance_output_sample.tsv')
        self.vcf_chr17 = os.path.join(
            os.sep, test_data_path,
            'test_data/chr17_1_83257441_paste.sites.vcf.gz')
        self.vcf_chr04 = os.path.join(
            os.sep, test_data_path,
            'test_data/chr4_1_190214555_paste.sites.vcf.gz')
        self.vcf_chr04_compromised = os.path.join(
            os.sep, test_data_path,
            'test_data/chr4_1_190214555_paste_SNPs_del.sites.vcf.gz')
        self.output_path = os.path.join(os.sep, test_data_path,
                                        'test_data/test_output.tsv')

    def tearDown(self):
        pass

    def test_process_output_tsv(self):
        d = process_output_tsv(self.tsv_file)
        pprint(d)

    def test_list2dict(self):
        L = []  # list to capture results
        with open(self.tsv_file, newline='') as csvfile:
            file_reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
            for row in file_reader:
                L.append(row[0])
        d = list2dict(L)

        self.assertEqual(d['type']['SNP']['true-positive'], '12')
        self.assertEqual(d['type']['INDEL']['sensitivity'], 'NaN')

    def test_run_concordance_same_file(self):
        """Testing perfect overlap by using the same file as truth and test."""
        run_concordance(reference=self.reference,
                        eval_file=self.vcf_chr17,
                        truth_file=self.vcf_chr17,
                        output_file=self.output_path)
        tf = process_output_tsv(self.output_path)

        self.assertEqual(tf, 0)

    def test_run_concordance_different_file(self):
        # Challenging the method by using files from two different chromosomes,
        # which likely have zero overlap. This test should fail.
        run_concordance(reference=self.reference,
                        eval_file=self.vcf_chr17,
                        truth_file=self.vcf_chr04,
                        output_file=self.output_path)
        tf = process_output_tsv(self.output_path)

        self.assertEqual(tf, 1)

    def test_run_concordance_SNPs_del(self):
        # Challenging the method by using two files from the same chromosome,
        # but in one of them 20% of SNPs have been deleted. Should fail.
        run_concordance(reference=self.reference,
                        eval_file=self.vcf_chr04_compromised,
                        truth_file=self.vcf_chr04,
                        output_file=self.output_path)
        tf = process_output_tsv(self.output_path)

        self.assertEqual(tf, 1)

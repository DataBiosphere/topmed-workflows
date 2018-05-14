#!/usr/bin/env python

import unittest
import json
import tempfile
import shutil
import os
import argparse
from utils.update_json_file import main, parse_args


def create_parser():
    parser = argparse.ArgumentParser()
    #parser.add_argument...
    # ...Create your parser as you like...
    return parser


class TestUpdateJsonFile(unittest.TestCase):

    def setUp(self):
        self.json_path = tempfile.mkdtemp()

        # Original data in JSON file.
        d = {'key1': '/home/user/old_path/file1',
             'key2': '/home/user/old_path/file2',
             'key3': '/home/user/old_path/file3'}

        with open(os.path.join(self.json_path, 'json_data.json'), 'w') as fp:
            json.dump(d, fp)

        self.parser = create_parser()

    def tearDown(self):
        shutil.rmtree(self.json_path)

    def test_parser(self):
        parser = parse_args(['-l', '-m', '-s'])
        #            ['input.json', 'old_path', 'new_path'])
        self.assertTrue(parser.long)

    def test_something(self):
        parsed = self.parser.parse_args(['--something', 'test'])
        self.assertEqual(parsed.something, 'test')

    def test_check_test_data(self):
        with open(os.path.join(self.json_path, 'json_data.json'), 'r') as fp:
            data = json.load(fp)

        self.assertEqual(data, {'key1': '/home/user/old_path/file1',
                                'key2': '/home/user/old_path/file2',
                                'key3': '/home/user/old_path/file3'}
                         )

    def test_update_json_file(self):


        main(
            os.path.join(os.path.join(self.json_path, 'json_data.json')),
            '/home/user/old_path/',
            '/home/user/new_path')
        with open(os.path.join(self.json_path, 'json_data.json'), 'r') as fp:
            data = json.load(fp)

        print(data)


if __name__=='__main__':
    unittest.main()


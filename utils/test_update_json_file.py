#!/usr/bin/env python

import unittest
import json
import tempfile
import shutil
import os
from utils.update_json_file import update_json_file

class TestUpdateJsonFile(unittest.TestCase):

    def setUp(self):
        self.json_path = tempfile.mkdtemp()

        # Original data in JSON file.
        d = {'key1': '/home/user/old_path/file1',
             'key2': '/home/user/old_path/file2',
             'key3': '/home/user/old_path/file3'}

        with open(os.path.join(self.json_path, 'json_data.json'), 'w') as fp:
            json.dump(d, fp)

    def tearDown(self):
        shutil.rmtree(self.json_path)

    def test_check_test_data(self):
        with open(os.path.join(self.json_path, 'json_data.json'), 'r') as fp:
            data = json.load(fp)

        self.assertEqual(data, {'key1': '/home/user/old_path/file1',
                                'key2': '/home/user/old_path/file2',
                                'key3': '/home/user/old_path/file3'}
)

    def test_update_json_file(self):

        new_path = '/home/user/new_path/'

        update_json_file(
            os.path.join(os.path.join(self.json_path, 'json_data.json')),
                                      '/home/user/old_path/',
                                      new_path)
        with open(os.path.join(self.json_path, 'json_data.json'), 'r') as fp:
            data = json.load(fp)

        print(data)




if __name__=='__main__':
    unittest.main()


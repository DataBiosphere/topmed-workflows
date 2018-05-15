#!/usr/bin/env python

import json

def update_json_file(json_file_name, old_base_path, new_base_path):
    """Substitutes a specific path (=value) in a JSON file with a new one.
    :param json_file_name: (str) absolute path
    :param new_base_path: (list) array of strings with new paths
    :return: writes JSON file back to path"""

    json_file = open(json_file_name, 'r')
    data = json.load(json_file)
    for key, val in data.items():
        if old_base_path in val:
            data[key] = val.replace(old_base_path, new_base_path)
    with open(json_file_name, 'w') as fp:
        json.dump(data, fp)

if __name__=='__main__':
    update_json_file()
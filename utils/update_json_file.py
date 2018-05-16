#!/usr/bin/env python

import json
import argparse
import sys
import os

def parse_args(args):
    """Helper method to parse command line arguments."""
    parser = argparse.ArgumentParser(
        prog=os.path.basename(__file__),
        description='Replace base path in JSON file values with new path.')
    parser.add_argument(
        "--input_JSON", type=str, dest='json_file_name', action='store',
        help="absolute path of JSON input file")
    parser.add_argument(
        "--old_path", type=str, dest='old_base_path', action='store',
        help="base path this method replaces in input file value")
    parser.add_argument(
        "--new_path", type=str, dest='new_base_path', action='store',
        help="new path that replaces the old base path")
    return parser.parse_args(args)


def main():
#def update_json_file(json_file_name, old_base_path, new_base_path):
    """Substitutes a specific path (=value) in a JSON file with a new one.
    :param json_file_name: (str) as absolute path
    :param old_base_path: (str) 
    :param new_base_path: (str)
    :return: writes JSON file back to path"""

    parser = parse_args(sys.argv[1:])

    # Open JSON input file, and update its keys.
    json_file = open(parser.json_file_name, 'r')
    data = json.load(json_file)
    for key, val in data.items():
        if parser.old_base_path in val:
            data[key] = val.replace(parser.old_base_path,
                                    parser.new_base_path)
    with open(parser.json_file_name, 'w') as fp:
        json.dump(data, fp, indent=4, sort_keys=True)

if __name__=='__main__':
    main()
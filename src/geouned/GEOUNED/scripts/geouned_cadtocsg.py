#!/usr/bin/env python3

"""
Convert Cad geometry to CSG format for use in Monte Carlo Codes
"""

import argparse

import geouned

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "-i",
    "--input",
    type=str,
    default="config.json",
    help="The path to the config JSON file",
)
args = parser.parse_args()


def main():
    geouned.CadToCsg.from_json(args.input)


if __name__ == "__main__":
    main()

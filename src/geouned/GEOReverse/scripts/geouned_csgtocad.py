#!/usr/bin/env python3

"""
Convert Cad geometry to CSG format for use in Monte Carlo Codes
"""

import argparse
import json

from pathlib import Path

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

    if not Path(args.input).exists():
        raise FileNotFoundError(f"config file {args.input} not found")

    with open(args.input) as f:
        config = json.load(f)

    geo = geouned.CsgToCad()
    geo.read_csg_file(**config["read_csg"])

    if "cell_filter" in config.keys():
        geo.cell_filter(**config["cell_filter"])

    if "mat_filter" in config.keys():
        geo.material_filter(**config["mat_filter"])

    geo.build_universe()

    if "export_cad" in config.keys():
        geo.export_cad(**config["export_cad"])
    else:
        geo.export_cad()


if __name__ == "__main__":
    main()

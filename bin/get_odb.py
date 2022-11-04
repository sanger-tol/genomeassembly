#!/usr/bin/env python3

# ----
# Version from Sanger Genomenote pipeline by @priyanka-surana
# https://github.com/sanger-tol/genomenote/blob/383f23e6b7a89f9aad6b85c8f7320b5c5825de73/bin/get_odb.py
# ----

import argparse
import os
import json
import sys
import requests
import re

def parse_args(args=None):
    Description = "Get ODB database value using GOAT API"

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument("TOLID", help="like icAdaBipu29")
    parser.add_argument("FILE_OUT", help="Output CSV file.")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    return parser.parse_args(args)

def make_dir(path):
    if len(path) > 0:
        os.makedirs(path, exist_ok=True)

def get_odb(tolid, file_out):
    # Using API, get JSON file with all associated ODB10 databases
    #response = requests.get("https://goat.genomehubs.org/api/v2/search?query=tax_lineage%28queryA.taxon_id%29&queryA=assembly--assembly_id%3D"+asm+"&result=taxon&fields=odb10_lineage").json()["results"]
    tolid = ''.join(filter(str.isalpha, tolid))
    response = requests.get("https://goat.genomehubs.org/api/v2/search?result=taxon&taxonomy=ncbi&includeEstimates=false&offset=0&summaryValues=count&query=tax_name%28"+tolid+"%29&fields=busco_lineage").json()["results"]
    # Extract ODB10 lineage values
#    odb_arr = [ r["result"]["fields"]["odb10_lineage"]["value"] for r in response ]
    odb_arr = [ r['result']['fields']['busco_lineage']['value'] for r in response ]
    # The closest [0] OBD10 lineage is selected, unless one of the lineage values is "eutheria", then choose "eutheria"
    odb_val = odb_arr[0] if odb_arr else "eukaryota_odb10"

    out_dir = os.path.dirname(file_out)
    make_dir(out_dir)

    with open(file_out, "w") as fout:
        print("busco_lineage", odb_val, sep=",", file=fout)

def main(args=None):
    args = parse_args(args)
    get_odb(args.TOLID, args.FILE_OUT)

if __name__ == "__main__":
    sys.exit(main())

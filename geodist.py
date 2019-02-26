#!/usr/bin/env python3

### This script was written by Bradley T. Martin, PhD candidate,
### University of Arkansas, Dept. of Biological Sciences,
### btm002@email.uark.edu
###
### It is intended to calculate pairwise geographic distances for
### input into the conStruct R package:
### https://cran.r-project.org/web/packages/conStruct/index.html

import argparse
import sys

import numpy as np
import pandas as pd

from geopy import distance
from pyproj import Proj



def main():

    arguments = Get_Arguments()

    lat = arguments.lat.strip().lower()
    long = arguments.long.strip().lower()

    df = pd.read_csv(arguments.coords, header=0)
    df.columns = [x.strip().lower() for x in df.columns]
    #print(df.columns)



    return 0

def Get_Arguments():
# Parse command-line arguments.
    parser = argparse.ArgumentParser(description="Calculates geodesic distance "
                                    "between samples given XY coordinates",
                                    add_help=False)

    required_args = parser.add_argument_group("Required Arguments")
    optional_args = parser.add_argument_group("Optional Arguments")

    ## Required Arguments
    required_args.add_argument("-c", "--coords",
                                type=str,
                                required=True,
                                help="String; Comma separated input file containing coordinates")
    required_args.add_argument("--lat",
                                type=str,
                                required=True,
                                help="String; Column name for latitude coordinates")
    required_args.add_argument("--long",
                                type=str,
                                required=True,
                                help="String; Column name for longitude coordinates")

    ## Optional Arguments
    optional_args.add_argument("-f", "--format",
                                type=str,
                                required=False,
                                default="dd",
                                nargs="?",
                                help="String; Specify input coordinate format: "
                                    "[ dd || utm ]; default=dd (decimal degrees)")
    optional_args.add_argument("-z", "--zone",
                                type=int,
                                required=False,
                                default=None,
                                nargs="?",
                                help="Int; Specify UTM zone if '--format utm' option "
                                "is specified; default=None")
    optional_args.add_argument("--zone_column",
                                type=str,
                                required=False,
                                default=None,
                                nargs="?",
                                help="String; Specify column name for UTM zones; "
                                "either --zone_column or --zone can be used, but not both")
    optional_args.add_argument("-h", "--help",
                                action="help",
                                help="Displays this help menu")


    args = parser.parse_args()

    if args.format == "utm" and (args.zone or args.zone_column) is None:
        parser.error("'--format utm' option requires --zone or --zone_column to be specified")

    if args.zone and args.zone_column:
        parser.error("--zone and --zone_column options cannot be used together")

    return args

if __name__ == "__main__":

    rtrn_code = main()
    print("Program finished with exit status " + str(rtrn_code) + "\n")
    sys.exit(rtrn_code)

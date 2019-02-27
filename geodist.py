#!/usr/bin/env python3

### This script was written by Bradley T. Martin, PhD candidate,
### University of Arkansas, Dept. of Biological Sciences,
### btm002@email.uark.edu
###
### It is intended to calculate pairwise geographic distances for
### input into the conStruct R package:
### https://cran.r-project.org/web/packages/conStruct/index.html

import argparse
import itertools
import operator
import sys

import numpy as np
import pandas as pd

from collections import defaultdict
from geopy import distance
from pyproj import Proj



def main():

    arguments = Get_Arguments()

    validate_file_exists(arguments.coords)
    validate_file_exists(arguments.popmap)

    popmap = read_popmap(arguments.popmap)
    popcounts = get_popcounts(popmap)

    popd = get_popdict(popmap)
    print(popd)

    lat = arguments.lat.strip().lower()
    long = arguments.long.strip().lower()

    df = pd.read_csv(arguments.coords, header=0)
    df.columns = [x.strip().lower() for x in df.columns]
    #print(df.columns)

    return 0

def get_popdict(list_of_tuples):
    # Input is list of tuples, returns dictionary of lists
    d = defaultdict(list)
    for v, k in list_of_tuples:
        d[k].append(v)
    return d

def get_popcounts(list_of_tuples):

    ind_list = list()
    popinfo = list()
    popcounts = list()
    d = defaultdict(str)

    for ind, pop in list_of_tuples:
        if not pop in popinfo:
            popinfo.append(pop)
            d[pop]=1
            ind_list.append(ind)
        else:
            d[pop]+=1

    print("\nDetected {} populations with sample counts as shown below:".format(len(popinfo)))

    for pop in popinfo:
        print("{}: {}".format(pop, d[pop]))
        my_tuple = (pop, d[pop])
        popcounts.append(my_tuple)

    return popcounts

def read_popmap(file):

    list_of_tuples = list()
    with open(file, "r") as fin:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            cols = line.split()
            my_tuple = (cols[0], cols[1])
            list_of_tuples.append(my_tuple)
    return list_of_tuples


def validate_file_exists(filename):

    try:
        file = open(filename, "r")
        file.close()
    except IOError:
        print("\nError: The file " + filename + " does not exist or could not be read.\n")
        sys.exit(1)

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
    required_args.add_argument("-p", "--popmap",
                                type=str,
                                required=True,
                                help="String; Filename for population map")
    required_args.add_argument("-i", "--id",
                                type=str,
                                required=True,
                                help="String; Column name for individual IDs")

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

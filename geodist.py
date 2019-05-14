#!/usr/bin/env python3

### This script was written by Bradley T. Martin, PhD candidate,
### University of Arkansas, Dept. of Biological Sciences
### Please submit any issues or bug reports to: btm002@email.uark.edu

### The script is intended to calculate pairwise geographic distances between
### populations for input into the conStruct R package:
### https://cran.r-project.org/web/packages/conStruct/index.html

import argparse
import sys
import utm

import numpy as np
import pandas as pd
import geopandas as gp
import matplotlib as mpl
import matplotlib.pyplot as plt

from collections import defaultdict, OrderedDict
from geopy import distance
from scipy.cluster.vq import kmeans, vq
from scipy.spatial.distance import cdist
from shapely.geometry import Point, Polygon

mpl.use("Agg")

def main():

    arguments = Get_Arguments()

    validate_file_exists(arguments.coords)
    validate_file_exists(arguments.popmap)
    popmap = read_popmap(arguments.popmap)
    popcounts = get_popcounts(popmap)

    id_col = arguments.id.lower().strip()

    popdf = pd.DataFrame.from_records(popmap, columns=[id_col, "POPDF"])
    popd = get_popdict(popmap)

    lat_col = arguments.lat.strip().lower()
    lon_col = arguments.long.strip().lower()

    df = pd.read_csv(arguments.coords, header=0)
    df.columns = [x.strip().lower() for x in df.columns]
    coord_format = arguments.format.lower().strip()

    if coord_format == "utm" and arguments.zone:
        zone = int(arguments.zone)
        latlon_list = convert_utm_to_latlon(df, lat_col, lon_col, arguments.zone, arguments.hemisphere)
        if latlon_list == 1:
            return 1

        df = add_latlon_to_df(df, latlon_list)

    elif coord_format == "utm" and arguments.zone_column:
        zone_col = arguments.zone_column.strip().lower()
        latlon_list = convert_utm_to_latlon(df, lat_col, lon_col, zone_col, arguments.hemisphere)
        if latlon_list == 1:
            return 1

        df = add_latlon_to_df(df, latlon_list)

    else:
        latlon_list = zip(df[lat_col].tolist(), df[lon_col].tolist())
        df = add_latlon_to_df(df, latlon_list)

    # Merge df and popdf on arguments.id. Will drop individuals not in popmap.
    df2 = df.merge(popdf, how="inner", on=id_col)
    df2.dropna(subset=["DD_LATCOL", "DD_LONCOL"], inplace=True)

    # Initialize coordinate reference system
    # Can be specified by user with --epsg option.
    e = "epsg:" + str(arguments.epsg)
    crs = {"init": e}
    geodf = coordinates2polygons(df2, crs)
    # Write polygons to shapefile if option toggled on.
    if arguments.polygons:
        write_shapefile(geodf, arguments.polygons)

    centroids_df = get_centroids(geodf, crs)

    # Write centroids to shapefile if option toggled on.
    if arguments.centroids:
        write_shapefile(centroids_df, arguments.centroids)

    # Calculate great circle distance using geopy.distance
    distances = calculate_greatcircledist(centroids_df)

    # Replace distances index with population IDs
    distances.set_index(centroids_df["POPDF"], inplace=True)

    # Write coordinates to CSV
    write_coords(centroids_df, arguments.coord_outfile)

    # Write matrix to CSV
    write_matrix(distances, arguments.geodist)

    return 0

def plot_kmeans(K, avgWithinSS):
    # Under development.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(K, avgWithinSS, "b*-")
    plt.grid(True)
    plt.xlabel("Number of Cluster")
    plt.ylabel("Average Within-cluster Sum of Squares")
    tt = plt.title("Elbow Plot for K-Means Clustering")
    fig.savefig("Kmeans.png")

def get_kmeans(maxK, df, latName, lonName):
    # Under development.
    # Get pd.DataFrame for lat/lon coordinates
    coords = df[[latName, lonName]].values

    # Specify all K values
    K = range(1, maxK+1)

    # Apply kmeans for K=1 to K=maxK
    # scipy.cluster.vq.kmeans
    KM = [kmeans(coords, k) for k in K]

    # Cluster centroids
    centroids = [cent for (cent,var) in KM]
    D_k = [cdist(coords, cent, "euclidean") for cent in centroids]

    cIdx = [np.argmin(D, axis=1) for D in D_k]
    dist = [np.min(D,axis=1) for D in D_k]

    # Get average sum of squares within K groups
    avgWithinSS = [sum(d)/coords.shape[0] for d in dist]

    return K, avgWithinSS, centroids

def msg(name=None):
    return """
    geodist.py
    [-c COORDS_FILE]
    [-p POPMAP_FILE]
    [-i IND_COLUMN]
    [--lat LAT_COLUMN]
    [--long LONG_COLUMN]
    [Any Optional Arguments]

    This program was written by Bradley T. Martin, PhD Candidate, University of Arkansas
    Please submit any issues or bug reports to btm002@email.uark.edu"""

def write_coords(gdf, outfile):
    """
    Writes comma-separated lon/lat coordinates to file
    Input:
        geopandas DataFrame with point objects (geopandas.GeoDataFrame)
    Returns:
        None
    """
    gdf["COORDS"] = gdf["geometry"].apply(lambda x: x.coords[0])
    gdf["Lon"] = gdf["COORDS"].apply(lambda x: x[0])
    gdf["Lat"] = gdf["COORDS"].apply(lambda x: x[1])
    columns = ["Lon", "Lat"]
    gdf.to_csv(outfile, sep="\t", header = True, index = False, columns = columns)


def write_matrix(dist, outfile):
    """
    Writes geographic distance matrix to output file as csv
    Input:
        Pandas DataFrame (pandas.DataFrame)
        Output filename obtained from argparse (String)
    Returns:
        None
    """
    dist.to_csv(outfile, sep="\t", header=False, index=True)

def calculate_greatcircledist(gdf):
    """
    Calculates Great Circle Distance from coordinates in GeoPandas GeoDataFrame
    Input:
        Geopandas DataFrame (GeoDataFrame)
    Returns:
        New pandas.DataFrame with Geographic Distance matrix (pandas.DataFrame)
    """
    list_of_tuples = list(map(lambda p: (p.y, p.x), gdf["geometry"]))
    result = [[distance.great_circle(p1, p2).km for p2 in list_of_tuples] for p1 in list_of_tuples]
    res_df = pd.DataFrame(result)
    return res_df

def get_centroids(gdf, crs):
    """
    Gets the centroids from shapely.geometry.Polygon objects
    Inputs:
        geopandas.GeoDataFrame object
        coordinate reference system object (set via command-line)
    Returns:
        New geopandas.GeoDataFrame with centroids as geometry object.
    """
    centroids = gdf.copy()
    centroids["geometry"] = centroids["geometry"].centroid
    centroids.crs = gdf.crs
    return centroids

def write_shapefile(gdf, shapefile):
    """
    Writes a shapefile from geopandas.GeoDataFrame objects
    Input:
        geopandas.GeoDataFrame (geopandas.GeoDataFrame)
        Shapefile name as defined by command-line argument
    Returns:
        Nothing.
    """
    gdf.to_file(driver = "ESRI Shapefile", filename = shapefile)

def coordinates2polygons(df, crs):
    """
    Converts lat/long columns from pandas dataframe to Polygon objects
    Grouped by population ID from popmap

    Input:
        pandas.DataFrame containing lat/lon coordinates columns (pandas.DataFrame)
        epsg geodetic parameter dataset ID to specify coordinate reference system (String)
            Default = 4269 (for NAD83)
    Returns:
        GeoPandas dataframe with polygon objects (geopandas.GeoDataFrame)
    """
    # Aggregates lat/lon coordinates into pandas.DataFrame column
    df["coordinates"] = list(zip(df.DD_LONCOL, df.DD_LATCOL))

    # Construct shapely.geometry Point objects from lat/lon
    df["geometry"] = df["coordinates"].apply(Point)

    # Construct geopandas.GeoDataFrame using shapely.geometry points
    gdf = gp.GeoDataFrame(df, geometry="geometry", crs=crs)

    # Extract the coordinates from the shapely.geometry.Point object
    gdf["geometry"] = gdf["geometry"].apply(lambda x: x.coords[0])

    # Group by population
    #   1. Get all coordinates for each population as a list
    #   2. Convert that list to shapely.geometry.Polygon objects
    gdf = gdf.groupby("POPDF", sort = False)["geometry"].apply(lambda x: Polygon(x.tolist())).reset_index()

    # Save result to new GeoDataFrame
    gdf = gp.GeoDataFrame(gdf, geometry="geometry", crs=crs)

    # This fixes invalid polygons
    gdf["geometry"] = gdf["geometry"].buffer(0.1)

    return gdf

def add_latlon_to_df(df, list_of_tuples):
    """
    Adds lat and lon columns to pandas dataframe from list of tuples
    Input:
        Pandas DataFrame (pandas.DataFrame)
    Returns:
        Modified pandas DataFrame (pandas.DataFrame)
    """
    lat = list()
    lon = list()
    for coord in list_of_tuples:
        lat.append(coord[0])
        lon.append(coord[1])
    df["DD_LATCOL"] = lat
    df["DD_LONCOL"] = lon
    return df

def convert_utm_to_latlon(df, latcol, loncol, zone, hemi):
    """
    Uses utm library to convert UTM coordinates to decimal degrees
    Input:
        pandas dataframe (pandas.DataFrame)
        latitude column name (string)
        longitude column name (string)
        UTM zone (int)
        UTM hemishpere (N or S) (string)
    Returns:
        list of tuples with converted coordinates (lat, lon) (list)
    """

    list_of_tuples = list()

    if isinstance(zone, int):
        for i in range(len(df[latcol])):
            if not np.isnan(df[loncol][i]) and not np.isnan(df[latcol][i]):
                try:
                    my_tuple = utm.to_latlon(df[loncol][i], df[latcol][i], zone, hemi)
                    list_of_tuples.append(my_tuple)
                except ValueError:
                    print("Error: The coordinates for individual {} do not fall within the specified UTM zone".format(i+1))
                    return 1
            else:
                print("Warning: Individual {} is missing coordinates; it be removed".format(str(i+1)))
                my_tuple = (df[latcol][i], df[loncol][i])
                list_of_tuples.append(my_tuple)

    else:
        if df[zone].isnull().values.any():
            print("Error: Zone value must be specified in coordinates file if using --zone_column option; one or more values is missing.")
            return 1
        for i in range(len(df[latcol])):
            if not np.isnan(df[loncol][i]) and not np.isnan(df[latcol][i]):
                try:
                    my_tuple = utm.to_latlon(df[loncol][i], df[latcol][i], df[zone][i], hemi)
                    list_of_tuples.append(my_tuple)
                except ValueError:
                    print("Error: The coordinates for individual {} do not fall within the specified UTM zone".format(i+1))
                    return 1
            else:
                print("Warning: Individual {} is missing coordinates".format(str(i+1)))
                my_tuple = (df[latcol][i], df[loncol][i])
                list_of_tuples.append(my_tuple)
    return list_of_tuples

def get_popdict(list_of_tuples):
    """
    This function is to use popID as dictionary keys, with lists of indIDs
    for each population as the values.
    Input:
        List of tuples containing (indID, popID) (list)
    Returns:
        Dictionary of lists with tuple[1] as key
    """
    d = defaultdict(list)
    for v, k in list_of_tuples:
        d[k].append(v)
    return d

def get_popcounts(list_of_tuples):
    """
    Get counts for each populations and return as list of tuples
    Input: list of tuples (list)
    Returns: list of tuples with (popID, indcounts) (list)
    """
    popinfo = list()
    popcounts = list()
    d = OrderedDict()

    for ind, pop in list_of_tuples:
        if not pop in popinfo:
            popinfo.append(pop)
            d[pop]=1
        else:
            d[pop]+=1

    print("\nDetected {} populations with sample counts as shown below:".format(len(popinfo)))

    for pop in popinfo:
        print("{}: {}".format(pop, d[pop]))
        my_tuple = (pop, d[pop])
        popcounts.append(my_tuple)
    print("\n")

    for pop in popinfo:
        assert int(d[pop]) >= 3, \
        "At least three individuals per population are required to create the polygons"


    return popcounts

def read_popmap(file):
    """
    Function to read a population map file in the format: indID\tpopID
    Input:
        filename (string)
    Returns:
        list of tuples containing (indID, popID) (list)
    """
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
    """
    Function to validate that an input file exists.
    Input:
        filename (string)
    Returns:
        None
    """
    try:
        file = open(filename, "r")
        file.close()
    except IOError:
        print("\nError: The file " + filename + " does not exist or could not be read.\n")
        sys.exit(1)

def Get_Arguments():
    """
    Parse command-line arguments. Imported with argparse.
    Returns: object of command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Calculates geodesic distance "
                                    "between samples given XY coordinates",
                                    add_help=False, usage=msg())

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
    optional_args.add_argument("--hemisphere",
                                type=str,
                                required=False,
                                default="N",
                                nargs="?",
                                help="String; Specify hemisphere (N || S) if --zone or "
                                "--zone_column is used; default=N")
    optional_args.add_argument("-e", "--epsg",
                                type=str,
                                required=False,
                                default="4269",
                                nargs="?",
                                help="String; Define epsg for spatial projection; "
                                "default=4269 (NAD83)")
    optional_args.add_argument("--polygons",
                                type=str,
                                required=False,
                                default=None,
                                nargs="?",
                                const="polygons",
                                help="String; Write population polygons to shapefile; "
                                "if toggled, specify shapefile name as string; default = off")
    optional_args.add_argument("--centroids",
                                type=str,
                                required=False,
                                default=None,
                                nargs="?",
                                const="centroids",
                                help="String; Write centroids of each population to shapefile; "
                                "if toggled, specify shapefile name as string; default = off")
    optional_args.add_argument("-g", "--geodist",
                                type=str,
                                required=False,
                                default="out.geodist.txt",
                                nargs="?",
                                help="String; Specify output filename for "
                                "geodist matrix; default = out.geodist.txt")
    optional_args.add_argument("--coord_outfile",
                                type=str,
                                required=False,
                                default="coords.out.txt",
                                nargs="?",
                                help="String; Specify output file for lat/lon "
                                "coordinates (Tab-delimited; "
                                "default = coords.out.txt")
    optional_args.add_argument("-h", "--help",
                                action="help",
                                help="Displays this help menu")


    if len(sys.argv)==1:
        print("\nExiting because no command-line options were called.\n")
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.format == "utm" and (args.zone or args.zone_column) is None:
        parser.error("'--format utm' option requires --zone or --zone_column to be specified")

    if args.zone and args.zone_column:
        parser.error("--zone and --zone_column options cannot be used together")

    return args

if __name__ == "__main__":

    rtrn_code = main()
    print("\nProgram finished with exit status " + str(rtrn_code) + "\n")
    sys.exit(rtrn_code)

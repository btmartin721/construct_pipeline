# construct_pipeline  
Pipeline for running conStruct on a distributed HPC system  

TO-DO:  
- Add clustering to geodist, for automatic pop assignment  

## Step 1 - geodist.py  

geodist.py calculates geographic distances between populations for input into the conStruct R package. It does so by first grouping by population IDs specified in a popmap file, then collapsing the individuals into distinct populations. It then calculates the great circle distances between each population and outputs two files: A distances file and a two-column file with lon/lat coordinates.  

### geodist.py dependencies  

```
python3 (tested with python 3.6)
geopandas
pandas
numpy
utm
geopy
shapely
```

The easiest way to install all the dependencies is with Anaconda or Miniconda. I have provided a conda environment file called geodist.yml that will allow you to create a conda environment containing all the dependencies. E.g.,  

`conda env create -f geodist.yml`  

Alternatively, you can install the dependencies yourself. E.g.,  

```
conda create --name geodist python=3.6
conda install -n geodist -c conda-forge geopandas 
conda install -n geodist -c conda-forge utm
conda install -n geodist -c conda-forge geopy
```

Installing geopandas using conda will install pandas, numpy, and shapely as well.  

I recommend using only getting the packages from conda-forge. I had compatibility issues when I tried getting packages from mixed repositories.  

### geodist.py options  

Use the help option to see full list of options:  

`geodist.py -h`  

__Required Arguments:__
```
-c [STRING], --coords [STRING]: Name of comma-separated file containing coordinates for each sample  
-i [STRING], --id [STRING]: Specify name of column containing sample IDs in --coords file.  
-p [STRING], --popmap [STRING]: Specify filename for tab-separated population map  
--lat [STRING]: Specify name of latitude column in --coords file  
--long [STRING]: Specify name of longitude column in --coords file  
```

__Optional Arguments:__
```
  -f [FORMAT], --format [STRING]:
                        Specify input coordinate format: [ dd || utm
                        ]; default=dd (decimal degrees); can be decimal degrees or utm
  -z [INT], --zone [INT]:
                        Specify UTM zone if '--format utm' option is
                        specified; default=None
  --zone_column [STRING]:
                        Specify column name for UTM zones in --coords file; either
                        --zone_column or --zone can be used, but not both
  --hemisphere [STRING]:
                        Specify hemisphere (N || S) if --zone or
                        --zone_column is used; default=N
  -e [EPSG], --epsg [STRING]:
                        Define epsg for spatial projection;
                        default=4269 (NAD83)
  --polygons [STRING]:
                        Write population polygons to shapefile;
                        can be name of shapefile or name of directory (created if !exists)
                        default = off; 
  --centroids [STRING]
                        Write centroids of each population to shapefile; 
                        can be name of shapefile or name of directory (created if !exists)
                        default = off
  -g [STRING], --geodist [STRING]
                        Specify output filename for geodist matrix between populations;
                        default = out.geodist.txt
  --coord_outfile [STRING]
                        Specify output file for lon/lat coordinates of populations
                        (Tab-delimited; default = coords.out.txt
  -h, --help            Displays this help menu
  ```
  
  ### Running geodist.py    
  
geodist.py takes as input:  

1. A CSV file containing coordinates. They coordinates can be in decimal degrees or UTM format with an easting or longitude and a northing or latitude column. If the UTM format is specified, you need to specify either the name of a column in the --coords file containing the UTM zones, or use the --zone option to specify a single UTM zone from the command-line. The script will then convert the UTM coordinates to decimal degrees.    

2. A two-column, tab-separated population map file with sampleID\tpopulationID as the columns. There should not be a header. E.g.,:  

```
sample1\tpop1
sample2\tpop1
sample3\tpop1
sample4\tpop2
sample5\tpop2
sample6\tpop2
```

The script will group by populationID and calculate great circle distances between the centroid points of each population. You can have it save shapefiles containing the polygons and centroids if you want to visualize them on a map.  

Finally, the script will output the centroid point coordinates for each population as a tab-separated file (lon\tlat), and it will output a geographic distance matrix that can be input into the runConstruct.R script.  

## Step 2 - runConstruct.R  

The next step in the pipeline is to run conStruct using the runConstruct.R script.  This script is intended to be run in batch mode from the command-line.  The script runs conStruct on a single K value and outputs the data to an RData file. Multiple K values can be run in parallel using GNU parallel, as shown in the included .pbs script (coming soon).  

### runConstruct.R dependencies  

The following R packages must be installed to use runConstruct.R and the following runXvalidation.R scripts:  

```
R version >3.4.0 (tested on R 3.5.1)
optparse
conStruct
doParallel
foreach
parallel
```

If you are working on a distributed HPC enviroment without sudo (root) privelidges, I highly recommend installing R using conda. You will also need to install r-essentials and r-base. E.g.:   

`conda create -n r_env r-essentials r-base`

Then you can install the packages that are available through conda:  

```
conda install -n r_env -c r r-doparallel
conda install -n r_env -c r r-foreach
```

conStruct and optparse are not available through conda, so you will have to install them within R:

```
R
install.packages("optparse", repos="https://cran.r-project.org")
install.packages("conStruct", repos="https://cran.r-project.org")
```

I have also included an R.yml file that you can use to install the R packages available on conda. But you'll still need to install conStruct and optparse using install.packages()

### runConstruct.R options  

Call the help menu with: 

`Rscript runConstruct.R -h`  

Options:  

```
        -s CHARACTER, --str=CHARACTER
                required; Input structure filename; default = NULL

        -p CHARACTER, --popmap=CHARACTER
                required; Input popmap filename; default = NULL

        -g CHARACTER, --geodist=CHARACTER
                required; Input geoDist matrix file (tab-delimited)

        -c COORDS, --coords=COORDS
                required; Input coordinates file (tab-delimited)

        -w CHARACTER, --wd=CHARACTER
                optional; Set working directory; default = ./

        --prefix=CHARACTER
                required; Specify prefix for output files in conStruct analysis

        --nchains=INTEGER
                optional; Specify number of independend MCMC chains to run

        -i INTEGER, --niter=INTEGER
                optional; Specify length of MCMC chain to run

        -K INTEGER, --K=INTEGER
                required; Specify K value (Number of popualations) to run

        -a CHARACTER, --afreq=CHARACTER
                optional; Specify allele frequency file to output

        -o CHARACTER, --outdir=CHARACTER
                optional; Specify directory to write output files

        -h, --help
                Show this help message and exit
```

### Running runConstruct.R  

Required inputs are:  

1. A STRUCTURE file. This file must have sampleIDs as the first column, population IDs as the second column, and the rest of the columns must contain the SNPs. There cannot be a header line.  

2. A tab-separated geographic distance matrix file. You can use the geodist.py script (see above) to calculate distances and output the matrix.  

3. A coordinates file.  This is also outputted from geodist.py.  Must be tab-separated file with longitude as the first column and latitude as the second.  

4. Prefix for all the output files.  

5. Population map (popmap) file. Should be the same one input into geodist.py.  This is to collapse all the individuals into populations. The resulting allele frequencies object created in the script will contain one line per population.  

I also strongly suggest changing the --niter option to a more appropriate (much longer) chain length. You can also tell it to run additional independent chains with the --nchains option.  

Finally, you can change the output directory if you want. I recommend it because a lot of files will be produced.  

I also strongly recommend that you assess the trace output files produced in this step before going on to the next one. This will help you to determine whether you ran the MCMC chains long enough by assessing convergence.  

Run the script as follows:

`Rscript runConstruct.R [options]`

## Step 3 - runXvalidation.R  

The runXvalidation.R script will run cross-validation for conStruct to evaluate whether 1) the spatial versus non-spatial model is better supported and 2) determine the "best" K value.  

### runXvalidation dependencies  

The dependencies for runConstruct.R are all required for runXvalidation.R as well. However, no additional dependences are required.  

### runXvalidation options  

Options:  
```
        -k INTEGER, --minK=INTEGER
                required; Specify minimum K value (lower case k)

        -K INTEGER, --maxK=INTEGER
                required; Specify maximum K value (upper case K)

        -r INTEGER, --nreps=INTEGER
                required; Specify number of cross-validation replicates

        -t NUMERIC, --trainProp=NUMERIC
                optional; Specify training proportion for cross-validation

        -n INTEGER, --nodes=INTEGER
                required; Specify number of CPU cores for parallelization

        -f, --saveFiles
                optional; Boolean; Toggle save Robj files from cross-validation; will save lots of files; default=FALSE

        -F, --saveFigs
                optional; Boolean; Don't save figures from cross-validation; default=TRUE

        -o CHARACTER, --outdir=CHARACTER
                optional; Specify directory for output files; will be created if doesn't exist

        -w CHARACTER, --wd=CHARACTER
                required; Specify directory containing output from runConstruct.R

        -h, --help
                Show this help message and exit
```

### Running runXvalidation.R

Required Input:  

1. You need to have run runConstruct.R first, which will write the .RData files to the output directory specified in the runConstruct.R options.  

2. You need to specify the minimum and maximum K values that you ran runConstruct.R on. Currently, minK must be 1.  It's kind of a deprecated option.  

3. You need to specify the directory where the runConstruct.R files are stored using the --wd option.  

4. You need to indicate the number of CPU threads to use with the --nodes option.  

5. Specify the number of cross-validation replicates with the --nreps option.  

The script will output lots of files, so it is recommended that you specify an output directory with the --outdir option. It is also  highly recommended that you inspect the cross-validation replicate traces for convergence. 

Finally, run the script as follows:  

`Rscript runXvalidation.R [options...]`  


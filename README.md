# construct_pipeline  
Pipeline for running conStruct on a distributed HPC system  

TO-DO:  
- Add clustering to geodist, for automatic pop assignment  

## Step 1 - geodist.py  

geodist.py calculates geographic distances between populations for input into the conStruct R package. It does so by first grouping by population IDs specified in a popmap file, then collapsing the individuals into distinct populations. It then calculates the distances between each population and outputs two files: A distances file and a two-column file with lon/lat coordinates.  

### geodist.py dependencies  

utm


#!/bin/bash

# Load required modules
module load PYTHON/2.7.5
module load lims/1.2

# Get information for each library (flow cell, lane, sample id, etc.)
# $1  needs to be the name of the project

project="CLARIAJOA"
/scratch/project/production/DAT/apps/LIMSQ/limsq -p $project | sed 's/;/\t/g' > ../data/${project}_info.tsv

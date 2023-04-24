#! /usr/bin/env python

"""
Author: "Laura Jiménez Gracia"
Date: 2020-10-12

This script merges performance metrics of cellranger counts [with Feature Barcoded]
for all libraries (GEM ids) in this project into 
a single file cellranger_mapping_metrics_counts.csv located in /results folder.
"""

# Load packages
import os
import pandas as pd
import config as cfg

# Get data paths
project = cfg.project
project_path = cfg.project_path
metadata_path = cfg.metadata_path

# Define list of libraries to merge 'metrics'
metadata_df = pd.read_csv(metadata_path, sep=",", header=0)
mask = (metadata_df["project"] == project)
libraries = metadata_df.loc[mask, "gem_id"]
libraries = list(set(libraries))

# Create results directory
results_directory = "{}/results".format(project_path)
if not os.path.exists(results_directory):
    os.mkdir(results_directory)


## FOR COUNTS
metrics_files = ["{}/jobs/{}/{}/outs/metrics_summary.csv".format(project_path, gem_id, gem_id) for gem_id in libraries]
# Combine metrics
combined_csv = pd.concat([pd.read_csv(f) for f in metrics_files], sort=False)
# Add gem_ids to metrics
combined_csv.insert(loc=0, column='gem_id', value=libraries)
# Save combined metrics
combined_csv.to_csv("{}/cellranger_mapping_metrics_count.csv".format(results_directory), index=False)
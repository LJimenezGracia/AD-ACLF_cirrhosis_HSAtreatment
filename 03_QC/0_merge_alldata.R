# ---
# author: "Laura Jiménez Gracia"
# date: 2021-24-03
# ---
# This R script merges the hashed and non-hashed Seurat objects into a single one.


# Load packages
library(Seurat)
library(tidyverse)

# Paths
path_r_objects <- here::here("02_demultiplexing/results/R_objects")
path_save_seuratobj <- here::here("03_QC/results/R_objects")

# Load data
seurat_obj_hashed <- readRDS(paste0(path_r_objects, "/CLARIAJOA_demultiplexed_margin1_merged.rds"))
seurat_obj_not_hashed <- readRDS(paste0(path_r_objects, "/CLARIAJOA_nothashed_merged.rds"))

# Homogenize metadata
seurat_obj_not_hashed$hash.ID <- NA
seurat_obj_hashed$treatment <- str_split_fixed(seurat_obj_hashed$hash.ID, "-", 2)[,2]

# Merge Seurat objects
seurat_obj <- merge(x = seurat_obj_hashed, y = seurat_obj_not_hashed)
Idents(seurat_obj) <- "library_name"
seurat_obj
rm(seurat_obj_hashed)
rm(seurat_obj_not_hashed)

# Remove unused metadata variables
seurat_obj$hash.ID <- NULL
seurat_obj$hashing_snr <- NULL

# Save Seurat object
saveRDS(seurat_obj, 
        file = paste0(path_save_seuratobj, "/CLARIAJOA_alldata_merged.rds"))

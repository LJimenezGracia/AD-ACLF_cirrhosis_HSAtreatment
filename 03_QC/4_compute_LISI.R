# ---
# author: Laura Jiménez Gracia
# date: 2021-02-24
# ---

# Compute Local Inverse Simpson Index (LISI)

# To evaluate the success of the different integration methods, we compute the LISI score
# running the LISI algorithm (https://www.nature.com/articles/s41592-019-0619-0)

# Pre-processing
## Load libraries
library(tidyverse)
library(lisi)

## Paths
path_qc_res_obj <- here::here("03_QC/results/R_objects")

## Functions
source(here::here("bin/utils.R"))

## Load data
unintegrated_coords <- readRDS(paste0(path_qc_res_obj, "/CLARIAJOA_pca_uncorrected.rds"))
harmony_coords <- readRDS(paste0(path_qc_res_obj, "/CLARIAJOA_pca_integrated_harmony.rds"))
confounders_df <- readRDS(paste0(path_qc_res_obj, "/CLARIAJOA_confounders.rds"))


# Compute LISI score
## Prepare list of PCA embeddings
pca_coords_list <- list(unintegrated_coords, harmony_coords)
names(pca_coords_list) <- list("Uncorrected", "Harmony")

## Compute LISI score to assess multiple confounder variables
lisi_scores <- purrr::map(pca_coords_list, function(pca_coords) {
  scores <- lisi::compute_lisi(
    X = pca_coords,
    meta_data = confounders_df,
    label_colnames = colnames(confounders_df),
    perplexity = 60
  )
})

## Merge LISI scores into a single table
lisi_scores_df <- bind_rows(lisi_scores, .id = "correction")

# Save LISI results
saveRDS(lisi_scores_df, 
        file = paste0(path_qc_res_obj, "/CLARIAJOA_lisi_scores.rds"))

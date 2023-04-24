# ---
# author: "Laura Jiménez Gracia"
# date: 2020-12-02
# ---
# This R script is used to run `1_demultiplex_hashtags.Rmd` Rmarkdown document.

# Load packages
library(knitr)
library(markdown)
library(rmarkdown)

# Paths
path_project_metadata <- here::here("01_cellranger_mapping/data/CLARIAJOA_metadata.csv")
path_demultiplexing_rmd <- here::here("02_demultiplexing/1_demultiplex_hashtags.Rmd")
path_demultiplexing_report <- here::here("02_demultiplexing/reports/")

# Load data
metadata <-  read.csv(path_project_metadata)
libraries_list <- as.vector(unique(metadata$gem_id[metadata$hashing == "hashed"]))

# Creating output directories
dir.create(file.path("results"), showWarnings = FALSE)
dir.create(file.path("results/R_objects"), showWarnings = FALSE)
dir.create(file.path("reports"), showWarnings = FALSE)

# Run demultiplexing .Rmd file
for (lib in libraries_list) {
  rmarkdown::render(path_demultiplexing_rmd,
                    output_file =  paste0(lib, "_demultiplexed_hashtags_margin1.html"),
                    output_dir = path_demultiplexing_report
  )
}

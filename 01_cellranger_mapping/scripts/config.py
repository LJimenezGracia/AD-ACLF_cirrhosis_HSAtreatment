#! /usr/bin/env python

"""
Python configuration file
"""

# SUBPROJECT INFO
project = "CLARIAJOA"
chemistry = "SC3Pv3"
reference = "human"

# Path to working project folder
project_path = f"/scratch/devel/ljimenez/projects/{project}/01_cellranger_mapping"

# Path to metadata file
metadata_path = f"/scratch/devel/ljimenez/projects/{project}/01_cellranger_mapping/data/{project}_metadata.csv"
infofile_path = f"/scratch/devel/ljimenez/projects/{project}/01_cellranger_mapping/data/{project}_info.tsv"
feat_ref_path = f"/scratch/devel/ljimenez/projects/{project}/01_cellranger_mapping/data/{project}_feature_reference.csv"

# Path to fastq files
fastq_path = "/scratch/project/production/fastq"

# Path to cellranger v5 - patch
cellranger_path = "/scratch/groups/singlecell/software/cellranger/6.1.1/cellranger"
slurmtemplate_path = "/scratch/groups/singlecell/software/cellranger/6.1.1/external/martian/jobmanagers/slurm.template"

# Path to reference genomes
Mmus_path = "/scratch/groups/singlecell/data/reference/refdata-gex-mm10-2020-A"
Mmus_vdj_path = "/scratch/groups/singlecell/data/reference/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0"

Hsapiens_path = "/scratch/groups/singlecell/data/reference/refdata-gex-GRCh38-2020-A"
Hsapiens_vdj_path = "/scratch/groups/singlecell/data/reference/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0"
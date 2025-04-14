# AD-ACLF_cirrhosis_HSAtreatment

This is the GitHub repository for the **Albumin reprograms the B cell transcriptional landscape and improves neutrophil antimicrobial function in patients with decompensated cirrhosis** manuscript.

The manuscript is out now at [JHEP Reports - Innovation in Hepatology EASL](https://doi.org/10.1016/j.jhepr.2024.101184).


## Abstract
**Background & Aims**: Patients with acutely decompensated (AD) cirrhosis are immunocompromised and particularly susceptible to infections. This study investigated the immunomodulatory actions of albumin by which this protein may lower the incidence of infections.

**Methods**: Blood immunophenotyping was performed in 11 patients with AD cirrhosis and 10 healthy volunteers (HV). Bulk and single-cell RNA sequencing (scRNA-seq) and flow cytometry were performed in peripheral blood mononuclear cells (PBMCs) from 20 patients with AD cirrhosis and 34 HV exposed to albumin. Albumin’s effects on degranulation, phagocytosis, chemotaxis, and swarming of neutrophils from six patients with AD cirrhosis and nine HV were assessed by measuring myeloperoxidase enzymatic activity, the engulfment of fluorescent-labeled Escherichia coli and zymosan, and interactions of neutrophils with Candida albicans at single-cell resolution in microfluidic chambers, respectively. Whole blood RNA sequencing (RNA-seq) analyses were performed in 49 patients admitted for severe AD cirrhosis, of whom 30 received albumin during hospitalization.

**Results**: Compared with HV, patients with AD cirrhosis showed severe lymphopenia and defective neutrophil antimicrobial function. Bulk and scRNA-seq analyses revealed significantly (false discovery rate [FDR] <0.05) increased signatures related to B cells, myeloid cells, and CD4+ T cells in PBMCs incubated with albumin. Changes in the B cell population were confirmed by flow cytometry. Neutrophils exposed to albumin also exhibited augmented chemotactic and degranulation responses, enhanced phagocytosis, and increased pathogen-restrictive swarming. RNA-seq data analysis in patients who had received albumin revealed specific upregulation of signatures related to B cells and neutrophils together with transcriptional changes in CD4+ T cells (FDR <0.05).

**Conclusions**: The finding that albumin promotes the transcriptional reprogramming and expansion of the B cell compartment and improves neutrophil antimicrobial functions indicates mechanisms that may lower the incidence of infections in patients with severe AD cirrhosis receiving albumin therapy.

**Impact and implications**: Patients with acutely decompensated cirrhosis receiving albumin as treatment have a lower incidence of infections. The reason for this protection is currently unknown, but the present study provides data that support the ability of albumin to boost the antimicrobial functions of immune cells in these patients. Moreover, these findings encourage the design of controlled clinical studies specifically aimed at investigating the effects of albumin administration on the immune system.


## Code implementation

The repository is organized into the following folder tree, which contains all the necessary files and scripts to perform the detailed tasks and reproduce all our results.

* **01_cellranger_mapping** --> It includes an overview of the project data information, including which samples and 10X libraries were generated. Also, it contains the scripts needed to create a folder directory to perform the sequencing read mapping to the reference genome. Finally, it includes R markdown notebooks to perform a general quality control on the raw sequencing reads from each library, considering different organism and tissue types independently.

* **02_demultiplexing** --> All scripts to demultiplex samples and to predict doublets for each 10x generated library. 
 
* **03_QC** --> R markdown notebooks to perform the quality control and the first pre-processing, including data normalization, scaling, dimensionality reduction and integration (and its evaluation).

* **04_clustering_annotation** --> All R markdown notebooks to perform a top-down clustering annotation approach, as well as scripts to find differential expressed markers for each clustering, and to assign a biological-relevant identity to each cluster.

* **05_GEX_analysis** --> All the code used to performed further downstream analysis on the processed data.
  
* **ACLF_figures** --> Scripts to generate the figures shown in the manuscript.


## Data accessibility

All data associated with this work is presented in the main manuscript or Supplementary material. 
For availability of any other type of data, contact the corresponding authors.

## Code accessibility

You can easily download a copy of all the files contained in this repository by:

* Cloning the git repository using the following command in the terminal:

`git clone https://github.com/LJimenezGracia/AD-ACLF_cirrhosis_HSAtreatment.git`

* Downloading a .ZIP archive [HERE](https://github.com/LJimenezGracia/AD-ACLF_cirrhosis_HSAtreatment/archive/refs/heads/main.zip).

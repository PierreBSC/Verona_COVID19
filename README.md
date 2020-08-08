# Analysis of the Verona COVID-19 cohort

This GitHub contains all the scripts used to analyze scRNA-seq and functional data from the paper 'Deciphering the state of immune silence in COVID-19 patients'.
Our single-cell analysis relies on the [**Pagoda2 pipeline**](https://github.com/kharchenkolab/pagoda2/) while the viral landscape of patients' lung was assessed using [**Viral-Track**](https://github.com/PierreBSC/Viral-Track/).

The ScRNA_seq_Verona_cohort_script.R file contains all the code required to generate the analysis and figures linked to scRNA-seq data. Many packages are required and have to be installed before running the script.

The Immune_suppression_analysis.R script details how to analyze the Immune-suppression data and the associated flow-cytometry and cytokine secretion data.

The Clinical_data_processing.R script is required for the analysis of the various biological and clinical data collected, including cytokine and blood cell count data.

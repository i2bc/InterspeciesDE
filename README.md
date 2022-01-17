# Repository to reproduce the analyses in Bastide et al. 2022.

## Content

* `InterspeciesDE.Rproj`: `R` project that can be opened with `R Studio` for convenience. 
All the `R` scripts assume that the working directory is the root directory of the project.

* `R_scripts`: `R` scripts to reproduce the simulations and plots.
The scripts are numbered in the order that they should be run.
Please see the header of each file for more a mode detailed description of each script.
  * `00_construct_rawTable_Sern2018.R`: format the raw counts table from Stern et al. 2018. DO NOT RUN.
  * `01_data_format_Stern2018.R`: format the data into a usable `PhyloCompData` object.
  * `02_simulations_Stern2018.R`: simulations using parameters drawn from the data.
  * `03_simulations_data_analysis.R`: differential expression analysis of the simulated datasets.
  * `04_simulations_data_analysis_results.R`: comparisons of the results from various methods on the same simulated datasets.
  * `05_plot_results.R`: plot of the results as presented in the manuscript.
  * `06_simulations_data_check.R`: check of the data using `countsimQC`.

* `data`: data issued from Stern et al. 2018.
  * `crayfish.nodelabels.tre`: crayfish tree from Stern et al 2018.
  * `states.csv`: states of each species (sight / no sight), from Stern et al 2018.
  * `rawCounts.txt`: raw counts table, produced by `R_scripts/00_construct_rawTable_Sern2018.R` (DO NOT RUN).
  * `rawLengths.txt`: raw length table, produced by `R_scripts/00_construct_rawTable_Sern2018.R` (DO NOT RUN).
  * `stern2018_cpd.rds`: formatted `PhyloCompData` dataset, produced by `R_scripts/01_data_format_Stern2018.R`.

* `2021-12-01_simulations_stern2018`: some datasets simulated using `R_scripts/02_simulations_Stern2018.R`, with base parameters (see text).
  * `simulation_parameters.RData`: all parameters used for the simulations.
  * `all_with_lengths_3_0.8_1_countsim_report.html`: reports obtained using the package `countsimQC` on selected simulated datasets.
  * `no_tree_alt_uni_NB_with_lengths_3_1_1_1.rds`: dataset simulated without the tree, using the NB framework.
  * `real_tree_alt_BM_with_lengths_3_0.8_1_1.rds`: dataset simulated using the BM on the tree, and the "alt" design.
  * `real_tree_bloc_BM_with_lengths_3_0.8_1_1.rds`: dataset simulated using the BM on the tree, and the "bloc" design.
  * `real_tree_sights_BM_with_lengths_3_0.8_1_1.rds`: dataset simulated using the BM on the tree, and the "sight" design.
  * `us_star_tree_alt_uni_BM_with_lengths_3_1_1_1.rds`: dataset simulated using the BM on a star tree.

* `2021-12-01_2021-12-01_simulations_stern2018_results`: selected results from differential analysis on the simulated datasets.
  * `full_result_table_1_50.rds`: scores of the various differential expression analysis methods applied on the various simulated datasets.
    Produced by the `R` script `R_scripts/04_simulations_data_analysis_results.R`.

## Requirements

* `R`: https://cran.r-project.org/index.html. Version 4.0 or later.
* `R` package `compcodeR`: https://doi.org/10.18129/B9.bioc.compcodeR.
* `R` package `countsimQC`: www.doi.org/10.18129/B9.bioc.countsimQC

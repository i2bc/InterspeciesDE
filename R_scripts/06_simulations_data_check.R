# Copyright (C) {2021} {PB, CS, MG}
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# Program to compute scores of the differential expression analysis methods.
# Use input from script 02_simulations_data_Stern2018.R

## Install correct version of compcodeR
# devtools::install_github("csoneson/compcodeR", ref = "phylocomp")

library(here)
library(compcodeR)
library(countsimQC)
library(DESeq2)

################################################################################
## File management
################################################################################
datestamp_day_simus <- "2021-12-01" ## Change here for the simulation date
name_data <- "stern2018"
simus_directory <- paste0(datestamp_day_simus, "_simulations_", name_data)

################################################################################
## Original Stern Data
################################################################################

data_stern <- readRDS(here("data", "stern2018_cpd.rds"))

data_stern <- DESeqDataSetFromMatrix(countData = data_stern@count.matrix,
                                     colData = data_stern@sample.annotations,
                                     design = ~condition)

################################################################################
## Simulation Parameters
################################################################################
load(here(simus_directory, "simulation_parameters.RData"))

## Base parameters
base_effect_size <- 3
base_prop_var_tree <- 0.8
base_fact_disp <- 1
base_cond_type <- ".*"
base_tree_type <- "real_tree"
base_model_process <- "BM"
base_use_lengths <- "with_lengths"

################################################################################
## Find base dataset simulation files
################################################################################

# List first replicate of all datasets
base_dataset <- paste(base_tree_type, base_cond_type, base_model_process, base_use_lengths, base_effect_size, base_prop_var_tree, base_fact_disp, sep = "_")
all_datasets <- list.files(path = here(simus_directory), pattern = paste0(base_dataset, "_1.rds"), full.names = FALSE)

# Add unstructured
us_datase <- paste("us_star_tree", base_cond_type, base_model_process, base_use_lengths, base_effect_size, 1, base_fact_disp, sep = "_")
all_datasets <- c(all_datasets, list.files(path = here(simus_directory), pattern = paste0(us_datase, "_1.rds"), full.names = FALSE))
no_dataset <- paste("no_tree", base_cond_type, "NB", base_use_lengths, base_effect_size, 1, base_fact_disp, sep = "_")
all_datasets <- c(all_datasets, list.files(path = here(simus_directory), pattern = paste0(no_dataset, "_1.rds"), full.names = FALSE))

################################################################################
## Analyse with countsimQC
################################################################################

base_dataset <- paste("all", base_use_lengths, base_effect_size, base_prop_var_tree, base_fact_disp, sep = "_")

get_results <- function(dataset){
  data_sim <- readRDS(here(simus_directory, dataset))
  
  return(DESeqDataSetFromMatrix(countData = data_sim@count.matrix,
                                colData = data_sim@sample.annotations,
                                design = ~condition))
}

data_list <- lapply(all_datasets, get_results)
names(data_list) <- all_datasets

data_list$Original <- data_stern

countsimQCReport(ddsList = data_list,
                 outputFile = here(simus_directory, paste0(base_dataset, "_countsim_report.html")),
                 outputDir = here(simus_directory),
                 outputFormat = "html_document", 
                 showCode = TRUE, forceOverwrite = TRUE,
                 savePlots = FALSE, description = NULL, 
                 maxNForCorr = 500, maxNForDisp = Inf, 
                 calculateStatistics = TRUE, subsampleSize = 50,
                 kfrac = 0.01, kmin = 5, 
                 permutationPvalues = FALSE, nPermutations = NULL)

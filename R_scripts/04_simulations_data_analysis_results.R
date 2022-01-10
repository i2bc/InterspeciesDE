# Copyright (C) {2021} {PB, MG}
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
# Use input from script 03_simulations_data_analysis.R

# NOTE : The full result table from our own run (result of 04_simulations_data_analysis.R)
#        is provided in the repository for convenience

## Install correct version of compcodeR
# devtools::install_github("csoneson/compcodeR", ref = "phylocomp")

library(here)
library(compcodeR)

################################################################################
## File management
################################################################################
# here_dir <- sub("save", "work", here())
here_dir <- here()

datestamp_day_simus <- "2021-11-25" ## Change here for the simulation date
datestamp_day_anaysis <- "2021-11-25" ## Change here for the analysis date
name_data <- "stern2018"
simus_directory <- paste0(datestamp_day_simus, "_simulations_", name_data)
simus_directory <- file.path(here_dir, simus_directory)
results_directory <- paste0(datestamp_day_anaysis, "_", datestamp_day_simus, "_simulations_", name_data, "_results")
results_directory <- file.path(here_dir, results_directory)

################################################################################
## Simulation Parameters
################################################################################
load(file.path(simus_directory, "simulation_parameters.RData"))

################################################################################
## Find all dataset simulation files
################################################################################
## List all files
all_datasets <- list.files(path = file.path(simus_directory), pattern = ".*.rds" )
## remove rds
all_datasets <- sub(".rds", "", all_datasets)
## Remove range
all_datasets <- sub("_[0-9]+$", "", all_datasets)
all_datasets <- unique(all_datasets)

################################################################################
## Create result table
################################################################################
## result table
result.table <- NULL

## Loop on all simulations
for (dataset in all_datasets) {
  
  ## Simulation data file
  data_file <- file.path(simus_directory, dataset)
  
  ## Find all analysis files that deal with the simulation dataset
  all_method_files <- list.files(path = results_directory,
                                 pattern = paste0(dataset, ".*.rds"),
                                 full.names = TRUE)
  
  ## If no analysis file, skip this dataset
  if (length(all_method_files) == 0) {
    warning(paste0("For dataset ", dataset, " there were no analysis file."))
    next
  }
  
  ## Extract range of simulation replicates
  range_N <- regmatches(all_method_files, regexpr(paste0(dataset, "_[0-9]+_"), all_method_files))
  range_N <- sub("_$", "", range_N)
  range_N <- as.numeric(regmatches(range_N, regexpr("[0-9]+$", range_N)))
  Nmin <- min(range_N)
  Nmax <- max(range_N)
  
  ##################################################################
  ## Comparison scores
  ##################################################################
  ## File table with all analyses for a given dataset
  file.table <- data.frame(input.files = all_method_files, stringsAsFactors = FALSE)
  
  ## Parameters and computed scores
  th <- 5e-2
  parameters <- list(incl.nbr.samples = NULL,
                     incl.replicates = NULL,
                     incl.dataset = dataset, 
                     incl.de.methods = NULL, 
                     fdr.threshold = th, tpr.threshold = th, typeI.threshold = th,
                     ma.threshold = th, fdc.maxvar = 1500, overlap.threshold = th,
                     fracsign.threshold = th, mcc.threshold = th, 
                     nbrtpfp.threshold = th, 
                     comparisons = c("auc",
                                     "mcc",
                                     "fdr", "tpr",
                                     "fdrvsexpr",
                                     "typeIerror", "fracsign", "nbrsign", "nbrtpfp",
                                     "maplot",
                                     "fdcurvesone",
                                     "rocone",
                                     "overlap", "sorensen",
                                     "scorevsexpr",
                                     "scorevssignal",
                                     "correlation"
                     ))
  
  ## Compute scores and save results in a table
  runComparison(file.table = file.table,
                parameters = parameters,
                output.directory = file.path(results_directory, dataset),
                save.result.table = TRUE,
                knit.results = TRUE)
  
  ## Read and format result table
  res.file <- list.files(file.path(results_directory, dataset),
                         pattern = "compcodeR_result_table_.*", full.names = TRUE)
  res.table <- readRDS(res.file[length(res.file)])
  
  ## Extract simu params
  find_tree <- function(dataset) {
    for (tt in all_tree_types) {
      if (grepl(tt, dataset)) return(tt)
    }
  }
  find_cond <- function(dataset) {
    for (cc in all_cond_types_tree[[find_tree(dataset)]]) {
      if (grepl(cc, dataset)) return(cc)
    }
  }
  find_model <- function(dataset) {
    for (mm in all_model_process_tree[[find_tree(dataset)]]) {
      if (grepl(mm, dataset)) return(mm)
    }
  }
  find_length <- function(dataset) {
    for (tt in all_use_lengths) {
      if (grepl(tt, dataset)) return(tt)
    }
  }
  find_effect_size <- function(dataset) {
    for (tt in all_effect_size) {
      if (grepl(paste(find_length(dataset), tt, sep = "_"), dataset)) return(tt)
    }
  }
  find_prop_tree <- function(dataset) {
    for (tt in rev(all_prop_var_tree)) {
      if (grepl(paste(find_length(dataset), find_effect_size(dataset), tt, sep = "_"), dataset)) return(tt)
    }
  }
  find_disp_factor <- function(dataset) {
    for (tt in all_fact_disp) {
      if (grepl(paste(find_length(dataset), find_effect_size(dataset), find_prop_tree(dataset), tt, sep = "_"), dataset)) return(tt)
    }
  }
  res.table$tree_type <- find_tree(dataset)
  res.table$cond_type <- find_cond(dataset)
  res.table$model_process <- find_model(dataset)
  res.table$use_lengths <- find_length(dataset)
  res.table$effect_size <- find_effect_size(dataset)
  res.table$prop_var_tree <- find_prop_tree(dataset)
  res.table$fact_disp <- find_disp_factor(dataset)
  
  ## Extract method parameters
  # methods
  res.table$method <- NA
  for (nn in c("DESeq2", "limma", "phylolimma", "phylolm")) {
    res.table$method[grepl(nn, res.table$de.methods)] <- nn
  }
  res.table$methodMod <- res.table$method
  # length Normalization
  res.table$lengthNormalization <- NA
  for (nn in c("none", "length", "TPM", "RPKM")) {
    res.table$lengthNormalization[grepl(paste0("\\.", nn, "\\."), res.table$de.methods)] <- nn
  }
  res.table$lengthNormalization[is.na(res.table$lengthNormalization)] <- "none"
  # transformation
  res.table$transformation <- NA
  for (nn in c("none", "log2", "sqrt", "asin(sqrt)")) {
    res.table$transformation[grepl(nn, res.table$de.methods)] <- nn
  }
  res.table$transformation[is.na(res.table$transformation)] <- "none"
  # Model process
  res.table$methodProcess <- NA
  for (nn in c("BM", "OUfixedRoot")) {
    res.table$methodProcess[grepl(nn, res.table$de.methods)] <- nn
    res.table$methodMod[grepl(nn, res.table$de.methods)] <- paste0(res.table$methodMod[grepl(nn, res.table$de.methods)], "_", nn)
  }
  # trend eBayes
  res.table$trendeBayes <- NA
  for (nn in c("no_trend", "with_trend")) {
    res.table$trendeBayes[grepl(nn, res.table$de.methods)] <- nn
    if (nn == "with_trend") res.table$methodMod[grepl(nn, res.table$de.methods)] <- paste0(res.table$methodMod[grepl(nn, res.table$de.methods)], "_trend")
  }
  # eBayes correction
  res.table$phyloeBayes <- NA
  res.table$phyloeBayes[grepl("moderation.eBayes", res.table$de.methods)] <- "eBayes"
  res.table$phyloeBayes[grepl("moderation.none", res.table$de.methods)] <- "vanilla"
  res.table$methodMod[grepl("moderation.none", res.table$de.methods)] <- paste0(res.table$methodMod[grepl("moderation.none", res.table$de.methods)], "_vanilla")
  # corBlock limma
  res.table$blockCor <- NA
  for (nn in c("id\\.species")) {
    res.table$blockCor[grepl(nn, res.table$de.methods)] <- "with_species_blocks"
    res.table$methodMod[grepl(nn, res.table$de.methods)] <- paste0(res.table$methodMod[grepl(nn, res.table$de.methods)], "_cor")
  }
  res.table$blockCor[is.na(res.table$blockCor)] <- "none"
  
  ## Effective Sample Size
  oneres <- readRDS(file.table$input.files[1])
  res.table$nEff <- oneres@info.parameters$nEff
  res.table$nEffRatio <- oneres@info.parameters$nEffRatio
  
  ## Merge
  result.table <- rbind(result.table, res.table)
  
  ## Temp file in case it crashes
  saveRDS(result.table, file = file.path(results_directory, paste0("tmp_full_result_table_", Nmin, "_", Nmax, ".rds")))
  
  ## Delete compcodeR file (not used)
  unlink(file.path(results_directory, dataset), recursive = TRUE)
}

## Save global table
saveRDS(result.table, file = file.path(results_directory, paste0("full_result_table_", Nmin, "_", Nmax, ".rds")))
##  Delete temp file
file.remove(file.path(results_directory, paste0("tmp_full_result_table_", Nmin, "_", Nmax, ".rds")))


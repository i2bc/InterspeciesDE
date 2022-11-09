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

# Program to run differential expression analysis on simulated data.
# Use input from script 02_simulations_Stern2018.R

# NOTE: This is a computation intensive script
# NOTE: The full result table from our own run (result of 04_simulations_data_analysis_results.R)
#       is provided in the repository for convenience

## Install correct version of compcodeR
# devtools::install_github("csoneson/compcodeR", ref = "phylocomp")

library(here)
library(compcodeR)
library(doParallel)

################################################################################
## Load Simulation Parameters
################################################################################
# here_dir <- sub("save", "work", here())
here_dir <- here()

name_data <- "stern2018"
datestamp_day_simus <- "2021-12-01" ## Change here to match simulation date
simus_directory <- paste0(datestamp_day_simus, "_simulations_", name_data)
simus_directory <- file.path(here_dir, simus_directory)
load(file.path(simus_directory, "simulation_parameters.RData"))

################################################################################
## File management
################################################################################
datestamp_day_anaysis <- format(Sys.time(), "%Y-%m-%d")
results_directory <- paste0(datestamp_day_anaysis, "_", datestamp_day_simus, "_simulations_", name_data, "_results")
results_directory <- file.path(here_dir, results_directory)
dir.create(results_directory)


################################################################################
## Parallel Computations settings
################################################################################
## Iterations to analyze
Nmin <- 1
Nmax <- 50

## Required packages to pass to all nodes
reqpckg <- c("compcodeR", "here", "foreach")

## Keep only some prop var tree
all_prop_var_tree <- list(0.6, 0.8, 1.0)

## Number of cores
Ncores <- 3

## Register nodes
cl <- makeCluster(Ncores, outfile = "")
registerDoParallel(cl)

################################################################################
## Loop on all settings
################################################################################

foreach (tree_type = all_tree_types) %:%
  foreach (cond_type = all_cond_types_tree[[tree_type]]) %:%
  foreach (use_lengths = all_use_lengths) %:%
  foreach (model_process = all_model_process_tree[[tree_type]]) %:%
  foreach (effect_size = all_effect_size) %:%
  foreach(prop_var_tree = all_prop_var_tree) %:%
  foreach(fact_disp = all_fact_disp, .packages = reqpckg, .verbose = TRUE) %dopar% {
    
    ## Parameters compatibility
    if (tree_type %in% c("us_star_tree", "no_tree") && prop_var_tree != 1.0) return(NULL)
    
    ## Remove some combinations of parameters
    if(tree_type == "ss_star_tree") return(NULL)
    if(use_lengths == "no_lengths") return(NULL)
    
    ## Dataset names
    dataset <- paste(tree_type, cond_type, model_process, use_lengths, effect_size, prop_var_tree, fact_disp, sep = "_")
    
    ## Data file name
    data_file <- file.path(simus_directory, dataset)
  
    ######################################################################
    ## DESeq2
    ######################################################################
    method <- "DESeq2"

    ## Length normalization when there are lengths
    all_length_norm <- c(FALSE)
    if (use_lengths == "with_lengths") all_length_norm <- c(FALSE, TRUE)

    foreach(lnorm = all_length_norm) %:%
      foreach(i = Nmin:Nmax, .packages = reqpckg, .verbose = TRUE) %do% {
        
        ## Remove some combinations of parameters
        if(!lnorm) return(NULL)

        ## Find the correct method to use (with or without length normalisation)
        method_name <- method
        if (lnorm) method_name <- paste0(method, ".length")

        ## Sanity check : if result file already exists, do not run the analysis
        if (file.exists(file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds")))) return(NULL)

        ## Run analysis
        runDiffExp(data.file = paste0(data_file, "_", i, ".rds"),
                   result.extent = method_name,
                   Rmdfunction = paste0(method_name, ".createRmd"),
                   output.directory = results_directory,
                   fit.type = "parametric",
                   test = "Wald",
                   nas.as.ones = TRUE)

        ## For the first iteration only, generate the analysis code
        if (i == 1) generateCodeHTMLs(file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds")), results_directory)
      }
    
    ####################################################################
    ## Limma
    ####################################################################
    method <- "lengthNorm.limma"
    
    ## Length normalization when there are lengths only
    all_length_norm <- c("none")
    if (use_lengths == "with_lengths") all_length_norm <- c("none", "RPKM", "TPM")
    ## Replicates correlation when there are replicates only
    all_blocks <- c("no_blocks")
    if (!tree_type %in% c("no_tree", "us_star_tree")) all_blocks <- c("no_blocks", "with_blocks")
    
    foreach (lnorm =  all_length_norm) %:%
      foreach (ltrans = c("log2", "sqrt")) %:%
      foreach (trend = c("no_trend", "with_trend")) %:%
      foreach (block = all_blocks) %:%
      foreach(i = Nmin:Nmax, .packages = reqpckg, .verbose = TRUE) %do% {
        
        ## Method id
        method_name <- paste(method, lnorm, ltrans, trend, block, sep = ".")

        ## If replicate correlation, find the right blocks
        is_block <- NULL
        if (block == "with_blocks") is_block <- "id.species"
        
        ## Sanity check : if result file already exists, do not run the analysis
        if (file.exists(file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds")))) return(NULL)
        
        ## Remove some combinations of parameters
        if(cond_type != "sights" && (lnorm != "TPM" || ltrans != "log2")) return(NULL)
        if(block != "with_blocks" && (lnorm != "TPM" || ltrans != "log2")) return(NULL)
        if(model_process == "OU" && (lnorm != "TPM" || ltrans != "log2" || block != "with_blocks")) return(NULL)
        if(trend != "no_trend" && (lnorm != "TPM" || ltrans != "log2" || cond_type != "sights")) return(NULL)
        
        ## Run analysis
        runDiffExp(data.file = paste0(data_file, "_", i, ".rds"),
                   result.extent = method_name,
                   Rmdfunction = paste0(method, ".createRmd"),
                   output.directory = results_directory,
                   norm.method = "TMM",
                   length.normalization = lnorm,
                   data.transformation = ltrans,
                   trend = (trend == "with_trend"),
                   block.factor = is_block)
        
        ## For the first iteration only, generate the analysis code
        if (i == 1) generateCodeHTMLs(file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds")), results_directory)
      }
    
    ##################################################################
    ## phylolm - with lengths
    ##################################################################
    method <- "phylolm"
    
    if (tree_type != "no_tree") { ## Run phylolm only when there is a tree
      
      ## Length normalization when there are lengths only
      all_length_norm <- c("none")
      if (use_lengths == "with_lengths") all_length_norm <- c("none", "RPKM", "TPM")
      
      foreach (sproc = c("BM", "OUfixedRoot")) %:%
        foreach (lnorm = all_length_norm) %:%
        foreach (ltrans = c("log2", "sqrt")) %:%
        foreach(i = Nmin:Nmax, .packages = reqpckg, .verbose = TRUE) %do% {
          
          ## Method id
          method_name <- paste(method, lnorm, ltrans, sproc, sep = ".")

          ## Sanity check : if result file already exists, do not run the analysis
          if (file.exists(file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds")))) return(NULL)
          
          ## Remove some combinations of parameters
          if(cond_type != "sights" && (lnorm != "TPM" || ltrans != "log2")) return(NULL)
          if(model_process == "OU" && (lnorm != "TPM" || ltrans != "log2")) return(NULL)

          ## Run analysis
          runDiffExp(data.file = paste0(data_file, "_", i, ".rds"),
                     result.extent = method_name,
                     Rmdfunction = paste0(method, ".createRmd"),
                     output.directory = results_directory,
                     norm.method = "TMM",
                     model = sproc,
                     measurement_error = TRUE,
                     extra.design.covariates = NULL,
                     length.normalization = lnorm,
                     data.transformation = ltrans)

          ## For the first iteration only, generate the analysis code
          if (i == 1) generateCodeHTMLs(file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds")), results_directory)
        }
      
      ####################################################################
      ## SVA + limma
      ####################################################################
      method <- "lengthNorm.sva.limma"
      
      ## Length normalization when there are lengths only
      all_length_norm <- c("none")
      if (use_lengths == "with_lengths") all_length_norm <- c("none", "RPKM", "TPM")
      ## Replicates correlation when there are replicates only
      all_n.sv <- c("auto", "one")
      if (!tree_type %in% c("no_tree", "us_star_tree")) all_blocks <- c("no_blocks", "with_blocks")
      
      foreach (lnorm =  all_length_norm) %:%
        foreach (ltrans = c("log2", "sqrt")) %:%
        foreach (trend = c("no_trend", "with_trend")) %:%
        foreach (n.sv = all_n.sv) %:%
        foreach(i = Nmin:Nmax, .packages = reqpckg, .verbose = TRUE) %do% {
          
          ## Method id
          method_name <- paste(method, lnorm, ltrans, trend, n.sv, sep = ".")
          
          ## If replicate correlation, find the right blocks
          # is_block <- NULL
          # if (block == "with_blocks") is_block <- "id.species"
          
          ## Sanity check : if result file already exists, do not run the analysis
          if (file.exists(file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds")))) return(NULL)
          
          ## Remove some combinations of parameters
          if(cond_type != "sights" && (lnorm != "TPM" || ltrans != "log2")) return(NULL)
          if(model_process == "OU" && (lnorm != "TPM" || ltrans != "log2")) return(NULL)
          if(trend != "no_trend") return(NULL)
          
          ## Run analysis
          runDiffExp(data.file = paste0(data_file, "_", i, ".rds"),
                     result.extent = method_name,
                     Rmdfunction = paste0(method, ".createRmd"),
                     output.directory = results_directory,
                     norm.method = "TMM",
                     length.normalization = lnorm,
                     data.transformation = ltrans,
                     trend = (trend == "with_trend"),
                     n.sv = ifelse(n.sv == "auto", "auto", 1))
          
          ## For the first iteration only, generate the analysis code
          if (i == 1) generateCodeHTMLs(file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds")), results_directory)
        }
    }
  }


## Stop the parallel cluster
stopCluster(cl)


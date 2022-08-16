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

# Program to re-analyse the empirical dataset from:
# Stern D. B., Crandall K. A. 2018. The Evolution of Gene Expression Underlying Vision Loss in Cave Animals. Molecular Biology and Evolution. 35:2005–2014.
# Use input from script 01_data_format_Stern2018.R

## Install correct version of compcodeR
# devtools::install_github("csoneson/compcodeR", ref = "phylocomp")

# NOTE: To reproduce the results from the paper, lines between "DO NOT RUN" tags can be ignored.
#       Required output analysis results are provided in the repository for convenience.
#       The html output of the comparison of all methods is also provided in the repository.

library(here)
library(compcodeR)
library(knitr)
library(ggplot2)
library(ape)
library(doParallel)

################################################################################
## Load Data
################################################################################
## read data
dataset <- "stern2018_cpd"
data_file <- here("data", dataset)

################################################################################
## Files management
################################################################################
# here_dir <- sub("save", "work", here())
here_dir <- here()

# datestamp_day <- format(Sys.time(), "%Y-%m-%d")
datestamp_day <- "2022-02-25"
results_directory <- paste0(datestamp_day, "_analysis_stern2018")
results_directory <- file.path(here_dir, results_directory)
dir.create(results_directory)

# BEGIN DO NOT RUN 

################################################################################
## Parallel Computations settings
################################################################################
## Required packages to pass to all nodes
reqpckg <- c("compcodeR", "here", "foreach")

## Number of cores
Ncores <- 3

## Register nodes
cl <- makeCluster(Ncores, outfile = "")
registerDoParallel(cl)

######################################################################
## DESeq2
######################################################################
method <- "DESeq2"

## Method name
method_name <- method
method_name <- paste0(method, ".length")

## Sanity check : if result file already exists, do not run the analysis
if (!file.exists(file.path(results_directory, paste0(dataset, "_", method_name, ".rds")))) {
  ## Run analysis
  runDiffExp(data.file = paste0(data_file, ".rds"),
             result.extent = method_name,
             Rmdfunction = paste0(method_name, ".createRmd"),
             output.directory = results_directory,
             fit.type = "parametric",
             test = "Wald",
             nas.as.ones = TRUE)
  
  ## Generate the analysis code
  generateCodeHTMLs(file.path(results_directory, paste0(dataset, "_", method_name, ".rds")), results_directory)
}

####################################################################
## Limma
####################################################################
method <- "lengthNorm.limma"

## Length normalization when there are lengths only
all_length_norm <- c("none", "RPKM", "TPM")
## Replicates correlation when there are replicates only
all_blocks <- c("no_blocks", "with_blocks")

foreach (lnorm =  all_length_norm) %:%
  foreach (ltrans = c("log2", "sqrt")) %:%
  foreach (trend = c("no_trend", "with_trend")) %:%
  foreach (block = all_blocks, .packages = reqpckg, .verbose = TRUE) %dopar% {
    
    ## Method id
    method_name <- paste(method, lnorm, ltrans, trend, block, sep = ".")
    
    ## If replicate correlation, find the right blocks
    is_block <- NULL
    if (block == "with_blocks") is_block <- "id.species"
    
    ## Sanity check : if result file already exists, do not run the analysis
    if (file.exists(file.path(results_directory, paste0(dataset, "_", method_name, ".rds")))) return(NULL)
    
    ## Run analysis
    runDiffExp(data.file = paste0(data_file, ".rds"),
               result.extent = method_name,
               Rmdfunction = paste0(method, ".createRmd"),
               output.directory = results_directory,
               norm.method = "TMM",
               length.normalization = lnorm,
               data.transformation = ltrans,
               trend = (trend == "with_trend"),
               block.factor = is_block)
    
    ## Generate the analysis code
    generateCodeHTMLs(file.path(results_directory, paste0(dataset, "_", method_name, ".rds")), results_directory)
  }

##################################################################
## phylolm - with lengths
##################################################################
method <- "phylolm"

## Length normalization when there are lengths only
all_length_norm <- c("none", "RPKM", "TPM")

foreach (sproc = c("BM", "OUfixedRoot")) %:%
  foreach (lnorm = all_length_norm) %:%
  foreach (ltrans = c("log2", "sqrt"), .packages = reqpckg, .verbose = TRUE) %dopar% {
    
    ## Method id
    method_name <- paste(method, lnorm, ltrans, sproc, sep = ".")
    
    ## Sanity check : if result file already exists, do not run the analysis
    if (file.exists(file.path(results_directory, paste0(dataset, "_", method_name, ".rds")))) return(NULL)
    
    ## Run analysis
    runDiffExp(data.file = paste0(data_file, ".rds"),
               result.extent = method_name,
               Rmdfunction = paste0(method, ".createRmd"),
               output.directory = results_directory,
               norm.method = "TMM",
               model = sproc,
               measurement_error = TRUE,
               extra.design.covariates = NULL,
               length.normalization = lnorm,
               data.transformation = ltrans)
    
    ## Generate the analysis code
    generateCodeHTMLs(file.path(results_directory, paste0(dataset, "_", method_name, ".rds")), results_directory)
  }

##################################################################
## Stop the parallel cluster
##################################################################
stopCluster(cl)

##################################################################
## Comparison scores
##################################################################
## Find all analysis files that deal with the simulation dataset
all_method_files <- list.files(path = results_directory,
                               pattern = paste0(dataset, ".*.TPM.log2.*.rds"),
                               full.names = TRUE)
all_method_files <- c(list.files(path = results_directory,
                                 pattern = paste0(dataset, ".*DESeq2.length.rds"),
                                 full.names = TRUE),
                      all_method_files)

## File table with all analyses for a given dataset
file.table <- data.frame(input.files = all_method_files, stringsAsFactors = FALSE)

## Parameters and computed scores
th <- 5e-2
parameters <- list(incl.nbr.samples = NULL,
                   incl.replicates = NULL,
                   incl.dataset = dataset,
                   incl.de.methods = NULL,
                   comparisons = c("maplot",
                                   "scorevsexpr",
                                   "fracsign",
                                   "nbrsign",
                                   "overlap",
                                   "sorensen",
                                   "correlation"
                   ))

## Compute scores and save results in a table
runComparison(file.table = file.table,
              parameters = parameters,
              output.directory = file.path(results_directory, dataset),
              save.result.table = TRUE,
              knit.results = TRUE)


# END DO NOT RUN 

##################################################################
## List of genes
##################################################################

## Results from Stern et al 2018
results_stern <- read.table(here("data", "Table_S2.txt"), sep = "\t", header = TRUE)

## All OGs
all_ogs <- read.table(here("data", "uniprot_sprot_blastout.tsv"), sep = "\t", header = FALSE)
all_ogs <- all_ogs[, c(1, 2)]
colnames(all_ogs) <- c("Orthogroup", "Uniprot.top.hit")

## Limma Cor
th <- 0.05
results_limma_cor <- readRDS(here(results_directory, "stern2018_cpd_lengthNorm.limma.TPM.log2.with_trend.with_blocks.rds"))
results_limma_cor <- subset(results_limma_cor@result.table, adjpvalue <= th)
results_limma_cor$Orthogroup <- rownames(results_limma_cor)
results_limma_cor$Uniprot.top.hit <- all_ogs[match(results_limma_cor$Orthogroup, all_ogs$Orthogroup), "Uniprot.top.hit"]
limma_OG <- results_limma_cor[, c("Orthogroup", "adjpvalue", "Uniprot.top.hit")]
limma_OG$Protein.name <- results_stern[match(limma_OG$Orthogroup, results_stern$Orthogroup), 4]
limma_OG$Protein.name[limma_OG$Uniprot.top.hit == "CSK2B_RAT"] <- "Casein kinase II subunit beta"
limma_OG
## Match Stern - Limma Cor
limma_stern <- results_stern[match(results_limma_cor$Orthogroup, results_stern$Orthogroup), -2]
limma_stern

## OU
th <- 0.05
results_OU <- readRDS(here(results_directory, "stern2018_cpd_phylolm.TPM.log2.OUfixedRoot.rds"))
results_OU <- subset(results_OU@result.table, adjpvalue <= th)
results_OU$Orthogroup <- rownames(results_OU)
results_OU
## Match Stern - OU
OU_stern <- results_stern[match(results_OU$Orthogroup, results_stern$Orthogroup), -2]
OU_stern

## OU vs limma cor
limma_OU <- results_OU[match(results_limma_cor$Orthogroup, results_OU$Orthogroup), -2]
limma_OU

## OG in all methods
limma_OU_stern <- results_stern[match(limma_OU$Orthogroup, results_stern$Orthogroup), -2]
limma_OU_stern
n_match_all <- nrow(na.omit(limma_OU_stern))
n_match_all

##################################################################
## List of genes - Other thresholds
##################################################################

## 0.1
th <- 0.1
results_limma_cor_001 <- readRDS(here(results_directory, "stern2018_cpd_lengthNorm.limma.TPM.log2.with_trend.with_blocks.rds"))
results_limma_cor_001 <- subset(results_limma_cor_001@result.table, adjpvalue <= th)
results_limma_cor_001$Orthogroup <- rownames(results_limma_cor_001)
results_limma_cor_001
# Match Stern - Limma Cor
limma_stern_001 <- results_stern[match(results_limma_cor_001$Orthogroup, results_stern$Orthogroup), -2]
limma_stern_001
diff_001 <- nrow(results_limma_cor_001) - nrow(results_limma_cor)
diff_001
match_001 <- nrow(na.omit(limma_stern_001)) - nrow(na.omit(limma_stern))
match_001

## 0.2
th <- 0.2
results_limma_cor_002 <- readRDS(here(results_directory, "stern2018_cpd_lengthNorm.limma.TPM.log2.with_trend.with_blocks.rds"))
results_limma_cor_002 <- subset(results_limma_cor_002@result.table, adjpvalue <= th)
results_limma_cor_002$Orthogroup <- rownames(results_limma_cor_002)
results_limma_cor_002
# Match Stern - Limma Cor
limma_stern_002 <- results_stern[match(results_limma_cor_002$Orthogroup, results_stern$Orthogroup), -2]
limma_stern_002
diff_002 <- nrow(results_limma_cor_002) - nrow(results_limma_cor)
diff_002
match_002 <- nrow(na.omit(limma_stern_002)) - nrow(na.omit(limma_stern))
match_002

# BEGIN DO NOT RUN 
####################################################################
## Limma by clades - Analysis
####################################################################
require(limma)
require(edgeR)
cdata <- readRDS(paste0(data_file, ".rds"))
is.valid <- check_phyloCompData(cdata)
if (!(is.valid == TRUE)) stop('Not a valid phyloCompData object.')

# Design
design_data <- cdata@sample.annotations
design_data$sight <- design_data$condition == 2
design_data$blind_Procambarus <- design_data$species %in% c("PLUCI", "PHORS", "PPALL")
design_data$blind_Cambarus <- design_data$species %in% c("CCRYP", "CHAMU", "CSETO")
design_data$blind_Orconectes <- design_data$species %in% c("OAUST", "OINCO")

# check
all_blind <- design_data$blind_Procambarus | design_data$blind_Cambarus | design_data$blind_Orconectes
all.equal(all_blind, design_data$condition == 1)

# clade design
design_formula <- as.formula(~ blind_Procambarus + blind_Cambarus + blind_Orconectes)
design_data <- design_data[, c("blind_Procambarus", "blind_Cambarus", "blind_Orconectes"), drop = FALSE]
design_data <- as.data.frame(lapply(design_data, function(x) as.factor(x + 0)))
design <- model.matrix(design_formula, design_data)

# Normalisation
nf <- edgeR::calcNormFactors(cdata@count.matrix / cdata@length.matrix, method = 'TMM')
lib.size <- colSums(cdata@count.matrix / cdata@length.matrix) * nf
data.norm <- sweep((cdata@count.matrix + 0.5) / cdata@length.matrix, 2, lib.size + 1, '/')
data.norm <- data.norm * 1e6

# Transformation
data.trans <- log2(data.norm)
rownames(data.trans) <- rownames(cdata@count.matrix)

# ANOVA fit
fit_anova <- lm(data.trans[1,] ~ ., data = design_data)
summary(fit_anova)

# Fitting Block correlations
block <- cdata@sample.annotations[['id.species']]
corfit <- duplicateCorrelation(data.trans, design = design, block = block, ndups = 1)

# Fit
length.fitlimma <- limma::lmFit(data.trans, design = design, correlation = corfit$consensus, block = block)
length.fitbayes <- limma::eBayes(length.fitlimma, trend = TRUE)
length.pvalues <- length.fitbayes$p.value[, -1]
length.adjpvalues <- as.data.frame(apply(length.pvalues, 2, function(x) p.adjust(x, method = 'BH')))
length.adjpvalues$Orthogroup <- rownames(length.adjpvalues)
length.logFC <- length.fitbayes$coefficients[, ncol(length.fitbayes$coefficients)]
length.score <- 1 - length.pvalues
result.table <- data.frame('pvalue' = length.pvalues, 'adjpvalue' = length.adjpvalues, 'logFC' = length.logFC, 'score' = length.score)
rownames(result.table) <- rownames(cdata@count.matrix)
cdata@result.table <- result.table
cdata@package.version <- paste('limma,', packageVersion('limma'), ';', 'edgeR,', packageVersion('edgeR'))
cdata@analysis.date <- date()
cdata@method.names <- list('short.name' = 'sqrtTPM', 'full.name' = 'length.3.48.3.limma.TMM.lengthNorm.TPM.dataTrans.log2.no_trend.id.species')
is.valid <- check_compData_results(cdata)
saveRDS(cdata, here(results_directory, "stern2018_cpd_lengthNorm.limma.TPM.log2.no_trend.with_blocks.byClade.rds"))

# END DO NOT RUN 

####################################################################
## Limma by clades - Genes
####################################################################

## Genes by clade
th <- 0.05
results_limma_cor <- readRDS(here(results_directory, "stern2018_cpd_lengthNorm.limma.TPM.log2.no_trend.with_blocks.byClade.rds"))
result.table <- results_limma_cor@result.table
result.table$Orthogroup <- rownames(result.table)

OG_Procambarus <- subset(result.table, adjpvalue.blind_Procambarus1 <= th)
OG_Cambarus <- subset(result.table, adjpvalue.blind_Cambarus1 <= th)
OG_Orconectes <- subset(result.table, adjpvalue.blind_Orconectes1 <= th)
OG_all <- subset(result.table,
                 adjpvalue.blind_Procambarus1 <= th & 
                   adjpvalue.blind_Cambarus1 <= th & 
                   adjpvalue.blind_Orconectes1 <= th)
OG_ProCam <- subset(result.table,
                    adjpvalue.blind_Procambarus1 <= th & 
                      adjpvalue.blind_Cambarus1 <= th)
OG_ProOrc <- subset(result.table,
                    adjpvalue.blind_Procambarus1 <= th & 
                      adjpvalue.blind_Orconectes1 <= th)
OG_CamOrc <- subset(result.table,
                    adjpvalue.blind_Cambarus1 <= th & 
                      adjpvalue.blind_Orconectes1 <= th)

OG_Procambarus <- OG_Procambarus[, c("Orthogroup", "adjpvalue.blind_Procambarus1"), drop = F]
OG_Cambarus <- OG_Cambarus[, c("Orthogroup", "adjpvalue.blind_Cambarus1"), drop = F]
OG_Orconectes <- OG_Orconectes[, c("Orthogroup", "adjpvalue.blind_Orconectes1"), drop = F]
OG_all <- OG_all[, c("Orthogroup"), drop = F]
OG_ProCam <- OG_ProCam[, c("Orthogroup"), drop = F]
OG_ProOrc <- OG_ProOrc[, c("Orthogroup"), drop = F]
OG_CamOrc <- OG_CamOrc[, c("Orthogroup"), drop = F]

## Match OG Uniprot top hit
OG_Procambarus$Uniprot.top.hit <- all_ogs[match(OG_Procambarus$Orthogroup, all_ogs$Orthogroup), "Uniprot.top.hit"]
OG_Cambarus$Uniprot.top.hit <- all_ogs[match(OG_Cambarus$Orthogroup, all_ogs$Orthogroup), "Uniprot.top.hit"]
OG_Orconectes$Uniprot.top.hit <- all_ogs[match(OG_Orconectes$Orthogroup, all_ogs$Orthogroup), "Uniprot.top.hit"]
OG_all$Uniprot.top.hit <- all_ogs[match(OG_all$Orthogroup, all_ogs$Orthogroup), "Uniprot.top.hit"]
OG_ProCam$Uniprot.top.hit <- all_ogs[match(OG_ProCam$Orthogroup, all_ogs$Orthogroup), "Uniprot.top.hit"]
OG_ProOrc$Uniprot.top.hit <- all_ogs[match(OG_ProOrc$Orthogroup, all_ogs$Orthogroup), "Uniprot.top.hit"]
OG_CamOrc$Uniprot.top.hit <- all_ogs[match(OG_CamOrc$Orthogroup, all_ogs$Orthogroup), "Uniprot.top.hit"]

## Get Uniprot names
cambarus_names <- read.table(here("data", "cambarus_list.protein_names.txt"), sep = "\t", header = TRUE)
orconectes_names <- read.table(here("data", "orconectes_list.protein_names.txt"), sep = "\t", header = TRUE)
procambarus_names <- read.table(here("data", "procambarus_list.protein_names.txt"), sep = "\t", header = TRUE)
colnames(cambarus_names) <- colnames(orconectes_names) <- colnames(procambarus_names) <- c("Uniprot.top.hit", "Protein.name")
stern_names <- results_stern[, c("Uniprot.top.hit", "Protein.name")]
all_names <- rbind(cambarus_names, orconectes_names, procambarus_names, stern_names)
all_names <- unique(all_names)

## Match Names
OG_Procambarus$Protein.name <- all_names[match(OG_Procambarus$Uniprot.top.hit, all_names$Uniprot.top.hit), "Protein.name"]
OG_Cambarus$Protein.name <- all_names[match(OG_Cambarus$Uniprot.top.hit, all_names$Uniprot.top.hit), "Protein.name"]
OG_Orconectes$Protein.name <- all_names[match(OG_Orconectes$Uniprot.top.hit, all_names$Uniprot.top.hit), "Protein.name"]
OG_all$Protein.name <- all_names[match(OG_all$Uniprot.top.hit, all_names$Uniprot.top.hit), "Protein.name"]
OG_ProCam$Protein.name <- all_names[match(OG_ProCam$Uniprot.top.hit, all_names$Uniprot.top.hit), "Protein.name"]
OG_ProOrc$Protein.name <- all_names[match(OG_ProOrc$Uniprot.top.hit, all_names$Uniprot.top.hit), "Protein.name"]
OG_CamOrc$Protein.name <- all_names[match(OG_CamOrc$Uniprot.top.hit, all_names$Uniprot.top.hit), "Protein.name"]

## Print result
OG_Procambarus
OG_Cambarus
OG_Orconectes
OG_all
OG_ProCam
OG_ProOrc
OG_CamOrc

## NAs
n_na <- sum(is.na(OG_Procambarus$Uniprot.top.hit)) + sum(is.na(OG_Cambarus$Uniprot.top.hit)) + sum(is.na(OG_Orconectes$Uniprot.top.hit))
n_na

# BEGIN DO NOT RUN 

####################################################################
## Separated setotus (extra analysis not used in the ms) - Analysis
####################################################################
require(limma)
require(edgeR)
cdata <- readRDS(paste0(data_file, ".rds"))
is.valid <- check_phyloCompData(cdata)
if (!(is.valid == TRUE)) stop('Not a valid phyloCompData object.')

# Design
design_data <- cdata@sample.annotations
design_data$sight <- design_data$condition == 2
design_data$blind_Procambarus <- design_data$species %in% c("PLUCI", "PHORS", "PPALL")
design_data$blind_Cambarus_Setosus <- design_data$species %in% c("CSETO")
design_data$blind_Cambarus <- design_data$species %in% c("CCRYP", "CHAMU")
design_data$blind_Orconectes <- design_data$species %in% c("OAUST", "OINCO")

# check
all_blind <- design_data$blind_Procambarus | design_data$blind_Cambarus_Setosus | design_data$blind_Cambarus | design_data$blind_Orconectes
all.equal(all_blind, design_data$condition == 1)

# clade design
design_formula <- as.formula(~ blind_Procambarus + blind_Cambarus_Setosus + blind_Cambarus + blind_Orconectes)
design_data <- design_data[, c("blind_Procambarus", "blind_Cambarus_Setosus", "blind_Cambarus", "blind_Orconectes"), drop = FALSE]
design_data <- as.data.frame(lapply(design_data, function(x) as.factor(x + 0)))
design <- model.matrix(design_formula, design_data)

# Normalisation
nf <- edgeR::calcNormFactors(cdata@count.matrix / cdata@length.matrix, method = 'TMM')
lib.size <- colSums(cdata@count.matrix / cdata@length.matrix) * nf
data.norm <- sweep((cdata@count.matrix + 0.5) / cdata@length.matrix, 2, lib.size + 1, '/')
data.norm <- data.norm * 1e6

# Transformation
data.trans <- log2(data.norm)
rownames(data.trans) <- rownames(cdata@count.matrix)

# ANOVA fit
fit_anova <- lm(data.trans[1,] ~ ., data = design_data)
summary(fit_anova)

# Fitting Block correlations
block <- cdata@sample.annotations[['id.species']]
corfit <- duplicateCorrelation(data.trans, design = design, block = block, ndups = 1)

# Fit
length.fitlimma <- limma::lmFit(data.trans, design = design, correlation = corfit$consensus, block = block)
length.fitbayes <- limma::eBayes(length.fitlimma, trend = TRUE)
length.pvalues <- length.fitbayes$p.value[, -1]
length.adjpvalues <- as.data.frame(apply(length.pvalues, 2, function(x) p.adjust(x, method = 'BH')))
length.adjpvalues$Orthogroup <- rownames(length.adjpvalues)
length.logFC <- length.fitbayes$coefficients[, ncol(length.fitbayes$coefficients)]
length.score <- 1 - length.pvalues
result.table <- data.frame('pvalue' = length.pvalues, 'adjpvalue' = length.adjpvalues, 'logFC' = length.logFC, 'score' = length.score)
rownames(result.table) <- rownames(cdata@count.matrix)
cdata@result.table <- result.table
cdata@package.version <- paste('limma,', packageVersion('limma'), ';', 'edgeR,', packageVersion('edgeR'))
cdata@analysis.date <- date()
cdata@method.names <- list('short.name' = 'sqrtTPM', 'full.name' = 'length.3.48.3.limma.TMM.lengthNorm.TPM.dataTrans.log2.no_trend.id.species')
is.valid <- check_compData_results(cdata)
saveRDS(cdata, here(results_directory, "stern2018_cpd_lengthNorm.limma.TPM.log2.no_trend.with_blocks.byCladeSetosus.rds"))

# END DO NOT RUN 

####################################################################
## Separated setotus (extra analysis not used in the ms) - Genes
####################################################################

## Genes by clade
th <- 0.05
results_limma_cor <- readRDS(here(results_directory, "stern2018_cpd_lengthNorm.limma.TPM.log2.no_trend.with_blocks.byCladeSetosus.rds"))
result.table <- results_limma_cor@result.table
result.table$Orthogroup <- rownames(result.table)

OG_Procambarus <- subset(result.table, adjpvalue.blind_Procambarus1 <= th)
OG_Cambarus <- subset(result.table, adjpvalue.blind_Cambarus1 <= th)
OG_Orconectes <- subset(result.table, adjpvalue.blind_Orconectes1 <= th)
OG_all <- subset(result.table,
                 adjpvalue.blind_Procambarus1 <= th & 
                   adjpvalue.blind_Cambarus1 <= th & 
                   adjpvalue.blind_Orconectes1 <= th)
OG_ProCam <- subset(result.table,
                    adjpvalue.blind_Procambarus1 <= th & 
                      adjpvalue.blind_Cambarus1 <= th)
OG_ProOrc <- subset(result.table,
                    adjpvalue.blind_Procambarus1 <= th & 
                      adjpvalue.blind_Orconectes1 <= th)
OG_CamOrc <- subset(result.table,
                    adjpvalue.blind_Cambarus1 <= th & 
                      adjpvalue.blind_Orconectes1 <= th)

OG_Procambarus <- OG_Procambarus[, c("Orthogroup", "adjpvalue.blind_Procambarus1"), drop = F]
OG_Cambarus <- OG_Cambarus[, c("Orthogroup", "adjpvalue.blind_Cambarus1"), drop = F]
OG_Orconectes <- OG_Orconectes[, c("Orthogroup", "adjpvalue.blind_Orconectes1"), drop = F]
OG_all <- OG_all[, c("Orthogroup"), drop = F]
OG_ProCam <- OG_ProCam[, c("Orthogroup"), drop = F]
OG_ProOrc <- OG_ProOrc[, c("Orthogroup"), drop = F]
OG_CamOrc <- OG_CamOrc[, c("Orthogroup"), drop = F]

## Match OG Uniprot top hit
OG_Procambarus$Uniprot.top.hit <- all_ogs[match(OG_Procambarus$Orthogroup, all_ogs$Orthogroup), "Uniprot.top.hit"]
OG_Cambarus$Uniprot.top.hit <- all_ogs[match(OG_Cambarus$Orthogroup, all_ogs$Orthogroup), "Uniprot.top.hit"]
OG_Orconectes$Uniprot.top.hit <- all_ogs[match(OG_Orconectes$Orthogroup, all_ogs$Orthogroup), "Uniprot.top.hit"]
OG_all$Uniprot.top.hit <- all_ogs[match(OG_all$Orthogroup, all_ogs$Orthogroup), "Uniprot.top.hit"]
OG_ProCam$Uniprot.top.hit <- all_ogs[match(OG_ProCam$Orthogroup, all_ogs$Orthogroup), "Uniprot.top.hit"]
OG_ProOrc$Uniprot.top.hit <- all_ogs[match(OG_ProOrc$Orthogroup, all_ogs$Orthogroup), "Uniprot.top.hit"]
OG_CamOrc$Uniprot.top.hit <- all_ogs[match(OG_CamOrc$Orthogroup, all_ogs$Orthogroup), "Uniprot.top.hit"]

## Get Uniprot names
cambarus_names <- read.table(here("data", "cambarus_list.protein_names.txt"), sep = "\t", header = TRUE)
orconectes_names <- read.table(here("data", "orconectes_list.protein_names.txt"), sep = "\t", header = TRUE)
procambarus_names <- read.table(here("data", "procambarus_list.protein_names.txt"), sep = "\t", header = TRUE)
colnames(cambarus_names) <- colnames(orconectes_names) <- colnames(procambarus_names) <- c("Uniprot.top.hit", "Protein.name")
stern_names <- results_stern[, c("Uniprot.top.hit", "Protein.name")]
all_names <- rbind(cambarus_names, orconectes_names, procambarus_names, stern_names)
all_names <- unique(all_names)

## Match Names
OG_Procambarus$Protein.name <- all_names[match(OG_Procambarus$Uniprot.top.hit, all_names$Uniprot.top.hit), "Protein.name"]
OG_Cambarus$Protein.name <- all_names[match(OG_Cambarus$Uniprot.top.hit, all_names$Uniprot.top.hit), "Protein.name"]
OG_Orconectes$Protein.name <- all_names[match(OG_Orconectes$Uniprot.top.hit, all_names$Uniprot.top.hit), "Protein.name"]
OG_all$Protein.name <- all_names[match(OG_all$Uniprot.top.hit, all_names$Uniprot.top.hit), "Protein.name"]
OG_ProCam$Protein.name <- all_names[match(OG_ProCam$Uniprot.top.hit, all_names$Uniprot.top.hit), "Protein.name"]
OG_ProOrc$Protein.name <- all_names[match(OG_ProOrc$Uniprot.top.hit, all_names$Uniprot.top.hit), "Protein.name"]
OG_CamOrc$Protein.name <- all_names[match(OG_CamOrc$Uniprot.top.hit, all_names$Uniprot.top.hit), "Protein.name"]

## Print result
OG_Procambarus
OG_Cambarus
OG_Orconectes
OG_all
OG_ProCam
OG_ProOrc
OG_CamOrc

## NAs
n_na <- sum(is.na(OG_Procambarus$Uniprot.top.hit)) + sum(is.na(OG_Cambarus$Uniprot.top.hit)) + sum(is.na(OG_Orconectes$Uniprot.top.hit))
n_na

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

# Program to simulate data using empirical dataset from:
# Stern D. B., Crandall K. A. 2018. The Evolution of Gene Expression Underlying Vision Loss in Cave Animals. Molecular Biology and Evolution. 35:2005–2014.
# Use input from script 01_data_format_Stern2018.R

# NOTE : This is a computation intensive script
# NOTE : The full result table from our own run (result of 04_simulations_data_analysis.R)
#        is provided in the repository for convenience

## Install correct version of compcodeR
# devtools::install_github("csoneson/compcodeR", ref = "phylocomp")

library(here)
library(compcodeR)
library(knitr)
library(ggplot2)
library(ape)

################################################################################
## Load Data
################################################################################
## read data
data_stern <- readRDS(here("data", "stern2018_cpd.rds"))

counts_noNA <- data_stern@count.matrix
length_noNA <- data_stern@length.matrix
colData <- data_stern@sample.annotations
colnames(colData)[colnames(colData) == "condition"] <- "sights"
tree <- data_stern@tree

################################################################################
## Files management
################################################################################
# here_dir <- sub("save", "work", here())
here_dir <- here()

datestamp_day <- format(Sys.time(), "%Y-%m-%d")
simus_directory <- paste0(datestamp_day, "_simulations_stern2018")
simus_directory <- file.path(here_dir, simus_directory)
dir.create(simus_directory)

################################################################################
## Analyse Data to get parameters
################################################################################
## run DESeq2 to get parameters estimation
library(DESeq2)
# Correct by families
dds <- DESeqDataSetFromMatrix(counts_noNA, colData, ~ sights )
# normalisation factors for library size
dds <-  estimateSizeFactors(dds)
size_fac <- sizeFactors(dds)
# compute normalization with length
mat_size_fac <- matrix(size_fac, ncol=length(size_fac), nrow=length(counts_noNA[,1]), byrow=T)
normFactors <- (mat_size_fac*length_noNA) / exp(rowMeans(log(mat_size_fac*length_noNA)))
normalizationFactors(dds) <- as.matrix(normFactors)
# analysis
dds <- DESeq2::DESeq(dds, fitType = 'parametric', test = 'Wald', betaPrior = TRUE)
dds <- DESeq2::replaceOutliersWithTrimmedMean(dds)
dds <- DESeq2::DESeq(dds, fitType = 'parametric', test = 'Wald', betaPrior = TRUE)
res <- DESeq2::results(dds, independentFiltering = TRUE, cooksCutoff = TRUE)

################################################################################
## Get empirical moments for simulation
################################################################################
dispersions_stern <- DESeq2::dispersions(dds)
seqdepth_stern <- mean(colSums(counts_noNA))
minfact_stern <- min(colSums(counts_noNA)/median(colSums(counts_noNA)))
maxfact_stern <- max(colSums(counts_noNA)/median(colSums(counts_noNA)))
relmeans_stern <- rowMeans(counts_noNA[,colData$sights == 2])
n.vars <- dim(counts_noNA)[1]

plot(log(relmeans_stern), log(dispersions_stern))

################################################################################
## Format Tree
################################################################################
n <- length(unique(sub("_.*", "", tree$tip.label))) # number of species
N <- length(tree$tip.label) # total number of observations (with replicates)

tree$edge.length <- tree$edge.length / max(diag(vcv(tree))[1:N]) # normalize tree to unit height

id_species <- sub("_.*", "", tree$tip.label)
id_species <- factor(id_species)
names(id_species) <- tree$tip.label

species_names <- sapply(1:length(colnames(length_noNA)), function(i) strsplit(colnames(length_noNA),"_")[[i]][1])
species_names <- sapply(1:length(species_names), function(i) strsplit(species_names,"[[:upper:]]$")[[i]][1])

################################################################################
## Format Lengths and compute empirical stats
################################################################################
unique_length <- length_noNA[, !duplicated(species_names)]

mean_lengths <- rowMeans(unique_length)
var_lengths <- apply(unique_length,1,var)
disp_lengths <- (var_lengths - mean_lengths) / mean_lengths^2 # size = mu^2/(var-mu), disp = 1 / size
disp_lengths[disp_lengths < 0] <- 1

plot(log(mean_lengths), log(var_lengths))
abline(0,1)

################################################################################
## Empirical fraction of data explained by the tree
################################################################################
# Normalisation TPM
nf <- edgeR::calcNormFactors(counts_noNA / length_noNA, method = 'TMM')
lib.size <- colSums(counts_noNA / length_noNA) * nf
data.norm <- sweep((counts_noNA + 0.5) / length_noNA, 2, lib.size + 1, '/')
data.norm <- data.norm * 1e6
# Transformation log2
data.trans <- log2(data.norm)
rownames(data.trans) <- rownames(counts_noNA)
# phylolm pagel lambda analysis
lambdas <- apply(data.trans, 1, function(x) phylolm::phylolm(x ~ colData$sights, phy = tree, model = "lambda")$optpar)
hist(lambdas)

################################################################################
## Parameters for the simulation
################################################################################
samples.per.cond <- N / 2
n.diffexp <- 150

Nrep <- 20 # number of replicates
all_effect_size <- 3.0

selection.strength <- log(2) / 0.5 # half life is 50% of tree height

all_prop_var_tree <- list("emp", 0.6, 0.8, 0.9, 1.0)

all_fact_disp <- c(1)

all_use_lengths <- c("with_lengths", "no_lengths")

################################################################################
## Conditions design
################################################################################
all_cond_types <- NULL
all_conds <- NULL

## sights
id_cond <- colData$sights
levels(id_cond) <- c(1, 2)
names(id_cond) <- rownames(colData)
# reorder
id_cond <- id_cond[match(tree$tip.label, names(id_cond))]

plot(tree)
tiplabels(pch = 21, col = id_cond, bg = id_cond)

all_cond_types <- c(all_cond_types, "sights")
all_conds[["sights"]] <- id_cond

## alternate tips
species_names <- unique(sub("_.*", "", tree$tip.label))
species_names[c(n-1, n)] <- species_names[c(n, n-1)]
cond_species <- rep(c(1, 2), length(species_names) / 2)
names(cond_species) <- species_names

id_cond <- id_species
id_cond <- cond_species[as.vector(id_cond)]
id_cond <- as.factor(id_cond)
names(id_cond) <- tree$tip.label

plot(tree)
tiplabels(pch = 21, col = id_cond, bg = id_cond)

all_cond_types <- c(all_cond_types, "alt")
all_conds[["alt"]] <- id_cond

## Blocs
species_names <- unique(sub("_.*", "", tree$tip.label))
cond_species <- rep(c(1, 2), each = length(species_names) / 2)
names(cond_species) <- species_names
cond_species["CHAMU"] <- 1

id_cond <- id_species
id_cond <- cond_species[as.vector(id_cond)]
id_cond <- as.factor(id_cond)
names(id_cond) <- tree$tip.label

plot(tree)
tiplabels(pch = 21, col = id_cond, bg = id_cond)

all_cond_types <- c(all_cond_types, "bloc")
all_conds[["bloc"]] <- id_cond

################################################################################
## Simulations 
################################################################################
all_tree_types <- NULL

tree_type <- "real_tree"
all_tree_types <- c(all_tree_types, tree_type)

all_cond_types_tree <- NULL
all_cond_types_tree[[tree_type]] <- all_cond_types
all_model_process_tree <- NULL
all_model_process_tree[[tree_type]] <- c("BM", "OU")

for (cond_type in all_cond_types_tree[[tree_type]]) {
  for (use_lengths in all_use_lengths) {
    for (model_process in all_model_process_tree[[tree_type]]) {
      for (effect_size in all_effect_size) {
        for(prop_var_tree in all_prop_var_tree) {
          for(fact_disp in all_fact_disp) {
            
            dataset <- paste(tree_type, cond_type, model_process, use_lengths, effect_size, prop_var_tree, fact_disp, sep = "_")
            
            set.seed(18570823)
            for (i in 1:Nrep) {
              simus_obj <- generateSyntheticData(dataset = dataset,
                                                 n.vars = n.vars,
                                                 samples.per.cond = samples.per.cond,
                                                 n.diffexp = n.diffexp,
                                                 repl.id = i,
                                                 seqdepth = seqdepth_stern,
                                                 minfact = minfact_stern,
                                                 maxfact = maxfact_stern,
                                                 relmeans = relmeans_stern,
                                                 dispersions = dispersions_stern * fact_disp,
                                                 fraction.upregulated = 0.5,
                                                 between.group.diffdisp = FALSE,
                                                 filter.threshold.total = 1,
                                                 filter.threshold.mediancpm = 0,
                                                 fraction.non.overdispersed = 0,
                                                 random.outlier.high.prob = 0,
                                                 random.outlier.low.prob = 0,
                                                 single.outlier.high.prob = 0,
                                                 single.outlier.low.prob = 0,
                                                 effect.size = effect_size,
                                                 output.file = file.path(simus_directory, paste0(dataset, "_", i, ".rds")),
                                                 tree = tree,
                                                 prop.var.tree = ifelse(prop_var_tree == "emp", lambdas, prop_var_tree),
                                                 id.condition = all_conds[[cond_type]],
                                                 id.species = id_species,
                                                 model.process = model_process,
                                                 selection.strength = ifelse(model_process == "OU", selection.strength, 0),
                                                 lengths.relmeans = if (use_lengths == "with_lengths") mean_lengths else NULL,
                                                 lengths.dispersions = if (use_lengths == "with_lengths") disp_lengths else NULL
              )
              
              if (i == 1) { ## Show summary only for the first iteration
                summarizeSyntheticDataSet(data.set = file.path(simus_directory, paste0(dataset, "_", i, ".rds")),
                                          output.filename = file.path(simus_directory, paste0(dataset, "_", i, "_check.html")))
              }
            }
          }
        }
      }
    }
  }
}

################################################################################
## Star tree - species structured
################################################################################

## Star tree - species structured
ss_star_tree <- tree
# Temp tips for virtual node
for (tip_label in c("PPALL_KC8786", "CCRYP_KC8782")) {
  ss_star_tree <- phytools::bind.tip(ss_star_tree, tip.label = paste0("TEMP_", tip_label),
                                     where = which(ss_star_tree$tip.label == tip_label))
}
ss_star_tree <- di2multi(ss_star_tree)
ss_star_tree$edge.length <- rep(0, length(ss_star_tree$edge.length)) # all branches to 0
nodes_before_tip <- unique(ss_star_tree$edge[ss_star_tree$edge[, 2] %in% 1:N, 1]) # nodes above a species
ss_star_tree$edge.length[ss_star_tree$edge[, 2] %in% nodes_before_tip] <- 1.0
ss_star_tree <- di2multi(ss_star_tree)
# remove temp tips
ss_star_tree <- ape::drop.tip(ss_star_tree, paste0("TEMP_", tree$tip.label))
plot(ss_star_tree)
is.ultrametric(ss_star_tree)
diag(vcv(ss_star_tree))[1:N]

## Simus
tree_type <- "ss_star_tree"
all_tree_types <- c(all_tree_types, tree_type)
cond_type <- "alt"
model_process <- "BM"

all_cond_types_tree[[tree_type]] <- cond_type
all_model_process_tree[[tree_type]] <- model_process

for (use_lengths in all_use_lengths) {
  for (effect_size in all_effect_size) {
    for(prop_var_tree in all_prop_var_tree) {
      for(fact_disp in all_fact_disp) {
        
        dataset <- paste(tree_type, cond_type, model_process, use_lengths, effect_size, prop_var_tree, fact_disp, sep = "_")
        
        set.seed(18570823)
        for (i in 1:Nrep) {
          simus_obj <- generateSyntheticData(dataset = dataset,
                                             n.vars = n.vars,
                                             samples.per.cond = samples.per.cond,
                                             n.diffexp = n.diffexp,
                                             repl.id = i,
                                             seqdepth = seqdepth_stern,
                                             minfact = minfact_stern,
                                             maxfact = maxfact_stern,
                                             relmeans = relmeans_stern,
                                             dispersions = dispersions_stern * fact_disp,
                                             fraction.upregulated = 0.5,
                                             between.group.diffdisp = FALSE,
                                             filter.threshold.total = 1,
                                             filter.threshold.mediancpm = 0,
                                             fraction.non.overdispersed = 0,
                                             random.outlier.high.prob = 0,
                                             random.outlier.low.prob = 0,
                                             single.outlier.high.prob = 0,
                                             single.outlier.low.prob = 0,
                                             effect.size = effect_size,
                                             output.file = file.path(simus_directory, paste0(dataset, "_", i, ".rds")),
                                             tree = ss_star_tree,
                                             prop.var.tree = ifelse(prop_var_tree == "emp", lambdas, prop_var_tree),
                                             id.condition = all_conds[[cond_type]],
                                             id.species = id_species,
                                             model.process = model_process,
                                             selection.strength = ifelse(model_process == "OU", selection.strength, 0),
                                             lengths.relmeans = if (use_lengths == "with_lengths") mean_lengths else NULL,
                                             lengths.dispersions = if (use_lengths == "with_lengths") disp_lengths else NULL
          )
          
          if (i == 1) { ## Show summary only for the first iteration
            summarizeSyntheticDataSet(data.set = file.path(simus_directory, paste0(dataset, "_", i, ".rds")),
                                      output.filename = file.path(simus_directory, paste0(dataset, "_", i, "_check.html")))
          }
        }
      }
    }
  }
}

################################################################################
## Star tree - unstructured
################################################################################

## Star tree - one ind per species
us_star_tree <- compcodeR:::prune_tree_one_obs(ss_star_tree)
plot(us_star_tree)
is.ultrametric(us_star_tree)
diag(vcv(us_star_tree))[1:n]

## alt_uni
all_cond_types <- c(all_cond_types, "alt_uni")
all_conds[["alt_uni"]] <- all_conds[["alt"]][us_star_tree$tip.label]

sum(all_conds[["alt_uni"]] == 1) == n / 2

## Simus
tree_type <- "us_star_tree"
all_tree_types <- c(all_tree_types, tree_type)
cond_type <- "alt_uni"
model_process <- "BM"

all_cond_types_tree[[tree_type]] <- cond_type
all_model_process_tree[[tree_type]] <- model_process

prop_var_tree <- 1.0

for (use_lengths in all_use_lengths) {
  for (effect_size in all_effect_size) {
    for(fact_disp in all_fact_disp) {
      
      dataset <- paste(tree_type, cond_type, model_process, use_lengths, effect_size, prop_var_tree, fact_disp, sep = "_")
      
      set.seed(18570823)
      for (i in 1:Nrep) {
        simus_obj <- generateSyntheticData(dataset = dataset,
                                           n.vars = n.vars,
                                           samples.per.cond = n / 2,
                                           n.diffexp = n.diffexp,
                                           repl.id = i,
                                           seqdepth = seqdepth_stern,
                                           minfact = minfact_stern,
                                           maxfact = maxfact_stern,
                                           relmeans = relmeans_stern,
                                           dispersions = dispersions_stern * fact_disp,
                                           fraction.upregulated = 0.5,
                                           between.group.diffdisp = FALSE,
                                           filter.threshold.total = 1,
                                           filter.threshold.mediancpm = 0,
                                           fraction.non.overdispersed = 0,
                                           random.outlier.high.prob = 0,
                                           random.outlier.low.prob = 0,
                                           single.outlier.high.prob = 0,
                                           single.outlier.low.prob = 0,
                                           effect.size = effect_size,
                                           output.file = file.path(simus_directory, paste0(dataset, "_", i, ".rds")),
                                           tree = us_star_tree,
                                           prop.var.tree = prop_var_tree,
                                           id.condition = all_conds[["alt_uni"]],
                                           id.species = id_species[us_star_tree$tip.label],
                                           check.id.species = FALSE,
                                           model.process = model_process,
                                           selection.strength = ifelse(model_process == "OU", selection.strength, 0),
                                           lengths.relmeans = if (use_lengths == "with_lengths") mean_lengths else NULL,
                                           lengths.dispersions = if (use_lengths == "with_lengths") disp_lengths else NULL
        )
        
        if (i == 1) { ## Show summary only for the first iteration
          summarizeSyntheticDataSet(data.set = file.path(simus_directory, paste0(dataset, "_", i, ".rds")),
                                    output.filename = file.path(simus_directory, paste0(dataset, "_", i, "_check.html")))
        }
      }
    }
  }
}

################################################################################
## Vanilla - unstructured
################################################################################

## Simus
tree_type <- "no_tree"
all_tree_types <- c(all_tree_types, tree_type)
cond_type <- "alt_uni"
model_process <- "NB"

all_cond_types_tree[[tree_type]] <- cond_type
all_model_process_tree[[tree_type]] <- model_process

prop_var_tree <- 1.0

for (use_lengths in all_use_lengths) {
  for (effect_size in all_effect_size) {
    for(fact_disp in all_fact_disp) {
      
      dataset <- paste(tree_type, cond_type, model_process, use_lengths, effect_size, prop_var_tree, fact_disp, sep = "_")
      
      set.seed(18570823)
      for (i in 1:Nrep) {
        simus_obj <- generateSyntheticData(dataset = dataset,
                                           n.vars = n.vars,
                                           samples.per.cond = n / 2,
                                           n.diffexp = n.diffexp,
                                           repl.id = i,
                                           seqdepth = seqdepth_stern,
                                           minfact = minfact_stern,
                                           maxfact = maxfact_stern,
                                           relmeans = relmeans_stern,
                                           dispersions = dispersions_stern * fact_disp,
                                           fraction.upregulated = 0.5,
                                           between.group.diffdisp = FALSE,
                                           filter.threshold.total = 1,
                                           filter.threshold.mediancpm = 0,
                                           fraction.non.overdispersed = 0,
                                           random.outlier.high.prob = 0,
                                           random.outlier.low.prob = 0,
                                           single.outlier.high.prob = 0,
                                           single.outlier.low.prob = 0,
                                           effect.size = effect_size,
                                           output.file = file.path(simus_directory, paste0(dataset, "_", i, ".rds")),
                                           id.condition = all_conds[[cond_type]],
                                           id.species = id_species[us_star_tree$tip.label],
                                           lengths.relmeans = if (use_lengths == "with_lengths") mean_lengths else NULL,
                                           lengths.dispersions = if (use_lengths == "with_lengths") disp_lengths else NULL,
                                           lengths.phylo = FALSE
        )
        
        if (i == 1) { ## Show summary only for the first iteration
          summarizeSyntheticDataSet(data.set = file.path(simus_directory, paste0(dataset, "_", i, ".rds")),
                                    output.filename = file.path(simus_directory, paste0(dataset, "_", i, "_check.html")))
        }
      }
    }
  }
}

save(Nrep,
     all_effect_size,
     all_use_lengths,
     all_tree_types,
     all_cond_types_tree,
     all_model_process_tree,
     all_conds,
     all_prop_var_tree,
     all_fact_disp,
     dispersions_stern,
     relmeans_stern,
     seqdepth_stern,
     minfact_stern,
     maxfact_stern,
     mean_lengths,
     disp_lengths,
     n.diffexp,
     id_species,
     tree,
     us_star_tree,
     ss_star_tree,
     file = file.path(simus_directory, "simulation_parameters.RData"))


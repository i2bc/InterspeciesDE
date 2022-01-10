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

# Program to format the data from:
# Stern D. B., Crandall K. A. 2018. The Evolution of Gene Expression Underlying Vision Loss in Cave Animals. Molecular Biology and Evolution. 35:2005–2014.
# Use input from script 00_construct_rawTable_Stern1018.R

## Install correct version of compcodeR
# devtools::install_github("csoneson/compcodeR", ref = "phylocomp")

library(here)
library(compcodeR)

################################################################################
## Read and format count and length data
################################################################################
# Read the raw counts
raw_counts <- read.table(here("data/rawCounts.txt"))
raw_counts <- round(raw_counts) # expected RSEM values -> counts
# Read length information
leng <- read.table(here("data/rawLengths.txt"))
# remove NA
raw_counts_noNA <- raw_counts[complete.cases(raw_counts),]
length_noNA <- leng[complete.cases(raw_counts),colnames(raw_counts_noNA)]

################################################################################
## Read and format condition data
################################################################################
condName <- "sights"
# Species condition (1 = blind, 2 = sighted)
states <- read.table(here("data/states.csv"),sep=",")
condstp <- states$V2
names(condstp) <- states$V1
# Sample Condition
idSpe <- sapply(1:dim(raw_counts)[2],function(i) strsplit(colnames(raw_counts)[i],"_")[[1]][1])
sights <- condstp[idSpe]  
names(sights) <- colnames(raw_counts)
# Format
sights <- factor(sights)
colData <- data.frame(sights)
colData$sights <- factor(colData$sights)
colData$species <- sub("_.*", "", rownames(colData))

################################################################################
## Read and format Tree
################################################################################
library(ape)
## Get the tree
tree <- read.tree(file = here("data/crayfish.nodelabels.tre"))
tree$node.label <- NULL
plot(tree)

## Species names
# Match
tree_data_cor <- match(tree$tip.label, colData$species)
data_tree_cor <- match(colData$species, tree$tip.label)
# Species in the tree NOT in data
tree$tip.label[is.na(tree_data_cor)]
# Species in data NOT in the tree
colData$species[is.na(data_tree_cor)]

## Format Tree
# Get rid of species not in data
tree <- drop.tip(tree, tip = tree$tip.label[is.na(tree_data_cor)])
plot(tree)
# function to add replicates
createTreeWithReplicates <- function(tree){
  tree_rep <- tree
  for (tip_label in tree$tip.label) {
    replis <- colnames(raw_counts_noNA)[grepl(tip_label, colnames(raw_counts_noNA))]
    for (rep in replis) {
      tree_rep <- phytools::bind.tip(tree_rep, tip.label = rep,
                                     where = which(tree_rep$tip.label == tip_label))
    }
  }
  # Remove original tips
  return(ape::drop.tip(tree_rep, tree$tip.label))
}
tree_rep <- createTreeWithReplicates(tree, tab)
# Plot
plot(tree_rep)

## Reorder data
corr <- match(tree_rep$tip.label, rownames(colData))
colData <- colData[corr, ]
raw_counts_noNA <- raw_counts_noNA[, corr]
length_noNA <- length_noNA[, corr]

################################################################################
# Format for compcodeR
################################################################################
## rename key column "condition" in colData
colnames(colData)[grep(condName, colnames(colData))] <- "condition"
## id species 
colData$id.species <- colData$species
## create object
cpd <- phyloCompData(tree = tree_rep,
                     count.matrix = raw_counts_noNA, 
                     sample.annotations = colData, 
                     info.parameters = list(dataset = "stern2018_cpd", uID = "1"),
                     length.matrix = as.matrix(length_noNA))
## check obj conformity
check_phyloCompData(cpd)

## save data to rds file
saveRDS(cpd, file = here("data", "stern2018_cpd.rds"))


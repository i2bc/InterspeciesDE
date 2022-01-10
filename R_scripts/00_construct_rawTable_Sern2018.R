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

# Program to download and format the data from:
# Stern D. B., Crandall K. A. 2018. The Evolution of Gene Expression Underlying Vision Loss in Cave Animals. Molecular Biology and Evolution. 35:2005–2014.
# https://github.com/TheDBStern/interspecific_rnaseq

# DO NOT RUN

library(here)
library(edgeR)

## Download data from Stern et al.
orthogroups <- "https://raw.githubusercontent.com/TheDBStern/interspecific_rnaseq/master/data/orthogroups.TMM.EXPR.matrix"
tab <- read.table(orthogroups)
write.table(tab, here("data/orthogroups.TMM.EXPR.matrix"))

tree <- "https://raw.githubusercontent.com/TheDBStern/interspecific_rnaseq/master/data/crayfish.nodelabels.tre"
tree <- ape::read.tree(tree)
write.table(tree, here("data/crayfish.nodelabels.tre"))

states <- "https://raw.githubusercontent.com/TheDBStern/interspecific_rnaseq/master/data/states.csv"
states <- read.table(states, sep = ",")
write.table(states, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE,
            here("data/states.csv"))

## raw expression matrix

# Raw data from Stern 2018 (personal comunication, available upon request)
filesToRead <- list.files(here("data/Stern_Expression_data_orthogroups"), pattern = "*.isoforms.results")

# Init tables
tab <- read.table(here("data/orthogroups.TMM.EXPR.matrix"))
tabRaw <- lengthRaw <- effLengthRaw <- TPMRaw <- tab
tabRaw[,] <- lengthRaw[,] <- effLengthRaw[,] <- TPMRaw[,]  <- NA

# Fill tables
for(i in 1:length(filesToRead)){
  print(i)
  current_sample <- strsplit(filesToRead[i],"[.]")[[1]][1]
  subTab <- read.table(here(paste0("data/Stern_Expression_data_orthogroups/",filesToRead[i])),header=TRUE)
  for(j in 1:length(rownames(tabRaw))){
    current_COG <- rownames(tabRaw)[j]
    ind_COG <- grep(current_COG,subTab$transcript_id)
    tabRaw[current_COG,current_sample] <- sum(subTab[ind_COG,"expected_count"])
    lengthRaw[current_COG,current_sample] <- sum(subTab[ind_COG,"length"])
    effLengthRaw[current_COG,current_sample] <- sum(subTab[ind_COG,"effective_length"])
    TPMRaw[current_COG,current_sample] <- sum(subTab[ind_COG,"TPM"])
  }
}

write.table(tabRaw,file=here("data/rawCounts.txt"))
write.table(lengthRaw,file=here("data/rawLengths.txt"))
write.table(effLengthRaw,file=here("data/rawEffLengths.txt"))
write.table(TPMRaw, file = here("data/rawTPM.txt"))
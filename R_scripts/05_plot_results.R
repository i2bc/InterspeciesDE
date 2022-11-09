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
# Use input from script 04_simulations_data_analysis_results.R
# The full result table from our own run is provided for convenience

## Install correct version of compcodeR
# devtools::install_github("csoneson/compcodeR", ref = "phylocomp")

library(compcodeR)
library(ggplot2)
library(dplyr)
library(here)

################################################################################
## File management
################################################################################
datestamp_day_simus <- "2021-12-01" ## Change here for the simulation date
datestamp_day_anaysis <- "2022-11-02" ## Change here for the analysis date
name_data <- "stern2018"
simus_directory <- paste0(datestamp_day_simus, "_simulations_", name_data)
simus_directory <- here(simus_directory)
results_directory <- paste0(datestamp_day_anaysis, "_", datestamp_day_simus, "_simulations_", name_data, "_results")
results_directory <- here(results_directory)

################################################################################
## Parameters
################################################################################
load(file.path(simus_directory, "simulation_parameters.RData"))
n.diffexp <- 150

## Iterations
Nmin <- 1
Nmax <- 50

## Base parameters
base_effect_size <- 3
base_prop_var_tree <- 0.8
base_fact_disp <- 1
base_cond_type <- "sights"

################################################################################
## Load and format result table for plotting
################################################################################
## Result table
result.table <- readRDS(file.path(results_directory, paste0("full_result_table_", Nmin, "_", Nmax, ".rds")))

## Replace NA in FDR when no discovery by 0
result.table[is.na(result.table$fdr) & !is.na(result.table$fp) & (result.table$fp == 0), "fdr"] <- 0.0

## Names
colnames(result.table)[colnames(result.table) == "fdr"] <- "FDR"
colnames(result.table)[colnames(result.table) == "tpr"] <- "TPR"
colnames(result.table)[colnames(result.table) == "auc"] <- "AUC"
colnames(result.table)[colnames(result.table) == "mcc"] <- "MCC"

result.table$methodMod <- sub("OUfixedRoot", "OU", result.table$methodMod)

result.table$Design <- result.table$cond_type
result.table$Design <- sub("sights", "sight", result.table$Design)
result.table$Design <- sub("bloc", "block", result.table$Design)
result.table$Design <- factor(result.table$Design, levels = c("block", "sight", "alt", "alt_uni"))

result.table$prop_var_tree <- as.factor(result.table$prop_var_tree)
result.table$effect_size <- as.factor(result.table$effect_size)

result.table$methodMod[grepl("n.sv.auto", result.table$de.methods)] <- "limma_sva_auto"
result.table$nsvmethod[grepl("n.sv.auto", result.table$de.methods)] <- "auto"

## Format for plot
result.table_plot <- reshape2::melt(result.table,
                                    measure.vars = c("AUC", "MCC", "FDR", "TPR"),
                                    variable.name = "score",
                                    value.name = "score_value")
result.table_plot$threshold <- NA
result.table_plot$threshold[result.table_plot$score == "FDR"] <- 0.05
result.table_plot$threshold[result.table_plot$score == "DR"] <- 1
result.table_plot$threshold[result.table_plot$score == "logDR"] <- 0
# result.table_plot$threshold[result.table_plot$score == "TD"] <- n.diffexp

pretty_names <- c("DESeq2" = "DESeq2",
                  "limma" = "limma",
                  "limma_cor" = "limma (cor)",
                  "limma_trend_cor" = "limma trend (cor)",
                  "limma_trend" = "limma trend",
                  "phylolm_BM" = "phylolm (BM)",
                  "phylolm_OU" = "phylolm (OU)",
                  "limma_sva_one" = "limma sva (one)",
                  "limma_sva_auto" = "limma sva (auto)")

pretty_names_tree <- c("no_tree" = "NB",
                       "us_star_tree" = "pPLN (star tree)",
                       "ss_star_tree" = "pPLN (species star tree)",
                       "real_tree" = "pPLN (real tree)")

plot_scores <- c("MCC", "FDR", "TPR")

################################################################################
## Figure : comparisons tree / no tree
################################################################################

df <- subset(result.table_plot,
             tree_type %in% c("real_tree",
                              "us_star_tree",
                              "no_tree")
             & lengthNormalization %in% c("TPM", "length")
             & use_lengths == "with_lengths"
             & transformation %in% c("log2", "none")
             & effect_size == base_effect_size
             & prop_var_tree %in% c(base_prop_var_tree, 1)
             & fact_disp == base_fact_disp
             & score %in% c("FDR")
             & model_process %in% c("NB", "BM")
             & methodMod %in% c("DESeq2", "limma")
)
df$methodMod <- factor(df$methodMod, levels = c("DESeq2", "limma"))
df$tree_type <- factor(df$tree_type, levels = c("no_tree",
                                                "us_star_tree",
                                                "real_tree"))
df <- subset(df, tree_type != "real_tree" | (tree_type == "real_tree" & cond_type == base_cond_type))
df <- subset(df, tree_type != "real_tree" | (tree_type == "real_tree" & prop_var_tree == base_prop_var_tree))
levels(df$tree_type) <- pretty_names_tree[levels(df$tree_type)]

if (nrow(df) / 3 / 2 / 1 != Nrep) message("wrong dimension")

ggplot(df, aes(x = tree_type, y = score_value, fill = methodMod, color = methodMod)) +
  geom_hline(aes(yintercept = threshold), linetype = 2) +
  facet_grid(score ~ ., scales = "free", switch = "y", labeller = labeller(Design = label_both)) +
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, show.legend = FALSE,
               position = "identity", width = 0.2*3/5,
               aes(group = interaction(effect_size, use_lengths, tree_type, methodMod))) +
  stat_summary(fun = median, geom = 'line', aes(group =  interaction(effect_size, use_lengths, methodMod))) +
  scale_fill_viridis_d(name = "Method", end = 0.8, option = "B", guide = "none") +
  scale_color_viridis_d(name = "Method",end = 0.8, option = "B", guide = guide_legend(override.aes = list(size = 2))) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size = 9),
        title = element_text(size = 7),
        strip.text.y = element_text(vjust = -3),
        strip.placement = "outside",
        strip.background = element_blank(),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position = c(0.01, 0.99),
        legend.justification = c("left", "top"),
        legend.key.size = unit(10, 'pt'),
        plot.margin = unit(c(0.1,0.1,-0.31,-0.5), "cm")
  )
# colorblindr::cvd_grid()


###############################################################################
## Figure : tree case, design and methods comparisons 
################################################################################

df <- subset(result.table_plot,
             lengthNormalization %in% c("TPM", "length")
             & transformation %in% c("log2", "none")
             & effect_size == base_effect_size
             & fact_disp == base_fact_disp
             & prop_var_tree == base_prop_var_tree
             & model_process == "BM"
             & score %in% plot_scores
             & !(methodMod %in% c(
               "limma_trend",
               "limma_trend_cor"))
)
df <- subset(df, use_lengths ==  "with_lengths" & tree_type == "real_tree")
df$methodMod <- factor(df$methodMod, levels = c("DESeq2", "limma", "limma_sva_one", "limma_sva_auto", "limma_cor", "phylolm_BM",
                                                "phylolm_OU"
))
levels(df$methodMod) <- pretty_names[levels(df$methodMod)]

if (nrow(df) / 3 / 3 / 7 != Nrep) message("wrong dimension")
if (anyNA(df$score_value)) message("Some NAs in the data")

dodge <- position_dodge(width = 0)

ggplot(df, aes(x = methodMod, y = score_value, fill = Design, color = Design)) +
  geom_hline(aes(yintercept = threshold), linetype = 2) +
  facet_grid(score ~ ., scales = "free", switch = "y", labeller = labeller(Design = label_both)) +
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, show.legend = FALSE,
               position = "identity", width = 0.2,
               aes(group = interaction(effect_size, use_lengths, tree_type, methodMod, prop_var_tree, Design))) +
  stat_summary(fun = median, geom = 'line', aes(group = interaction(effect_size, use_lengths, tree_type, prop_var_tree, Design), colour = Design)) +
  scale_fill_viridis_d(guide = "none", end = 0.8, option = "B") +
  scale_colour_viridis_d(end = 0.8, option = "B", guide = guide_legend(override.aes = list(size = 2))) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size = 9),
        title = element_text(size = 7),
        strip.text.y = element_text(vjust = -3),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.01, 0.01),
        legend.justification = c("left", "bottom"),
        legend.key.size = unit(10, 'pt'),
        plot.margin = unit(c(0.1,0.1,-0.31,-0.5), "cm")
  )
# colorblindr::cvd_grid()

################################################################################
## Figure : OU / BM comparisons
################################################################################

df <- subset(result.table_plot,
             lengthNormalization %in% c("TPM", "length")
             & transformation %in% c("log2", "none")
             & effect_size == base_effect_size
             & prop_var_tree == base_prop_var_tree
             & fact_disp == base_fact_disp
             # & cond_type == base_cond_type
             & score %in% plot_scores
             & methodMod %in% c(
               "phylolm_BM", "phylolm_OU", "limma_cor")
)
df <- subset(df, use_lengths ==  "with_lengths" & tree_type == "real_tree")
df$methodMod <- factor(df$methodMod,
                       levels = c("DESeq2", "limma",
                                  "limma_cor",
                                  "phylolm_BM",
                                  "phylolm_OU"
                       ))

levels(df$methodMod) <- pretty_names[levels(df$methodMod)]

if (nrow(df) / 3 / 3 / 3 / 2 != Nrep) message("wrong dimension")
if (anyNA(df$score_value)) message("Some NAs in the data")

ggplot(df, aes(x = methodMod, y = score_value, fill = model_process, color = model_process)) +
  facet_grid(score ~ cond_type, scales = "free_y", switch = "y", labeller = labeller(Design = label_both)) +
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, show.legend = FALSE,
               position = "identity", width = 0.2*3/5,
               aes(group = interaction(effect_size, methodMod, Design, model_process, factor(fact_disp), prop_var_tree))) +
  stat_summary(fun = median, geom = 'line', aes(group = interaction(effect_size, Design, model_process, factor(fact_disp), prop_var_tree))) +
  scale_fill_viridis_d(end = 0.8, begin = 0.4, option = "B", guide = "none") +
  scale_color_viridis_d(name = "", end = 0.8, begin = 0.4, option = "B", guide = guide_legend(override.aes = list(size = 2))) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size = 9),
        title = element_text(size = 7),
        strip.text.y = element_text(vjust = -3),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.01, 0.01),
        legend.justification = c("left", "bottom"),
        legend.direction = "horizontal",
        legend.key.size = unit(10, 'pt'),
        plot.margin = unit(c(0.1,0.1,-0.31,-0.5), "cm")
  )

# colorblindr::cvd_grid()

################################################################################
## Figure : variation of tree / noise variance proportion
################################################################################

df <- subset(result.table_plot,
             lengthNormalization %in% c("TPM", "length")
             & transformation %in% c("log2", "none")
             & effect_size == base_effect_size
             & prop_var_tree %in% c(0.6, 0.8, 1.0)
             & fact_disp == base_fact_disp
             & cond_type == base_cond_type
             & score %in% c("MCC")
             & model_process == "BM"
             & methodMod %in% c("limma_cor",
                                "phylolm_BM", "phylolm_OU")
)
df <- subset(df, use_lengths ==  "with_lengths" & tree_type == "real_tree")
df$methodMod <- factor(df$methodMod,
                       levels = c("limma_cor",
                                  "phylolm_BM", "phylolm_OU"
                       ))

levels(df$methodMod) <- pretty_names[levels(df$methodMod)]

if (nrow(df) / 3 / 1 / 3 != Nrep) message("wrong dimension")

df$prop_var_noise = df$prop_var_tree
levels(df$prop_var_noise) <- c(40, 20, 10, 0, "emp")
df$prop_var_noise <- relevel(df$prop_var_noise, "emp")

ggplot(df, aes(x = methodMod, y = score_value, fill = prop_var_noise, color = prop_var_noise)) +
  facet_grid(score ~ ., scales = "free_y", switch = "y", labeller = labeller(Design = label_both)) +
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, show.legend = FALSE,
               position = "identity", width = 0.2*3/5,
               aes(group = interaction(effect_size, methodMod, Design, model_process, factor(fact_disp), prop_var_noise))) +
  stat_summary(fun = median, geom = 'line', aes(group = interaction(effect_size, Design, model_process, factor(fact_disp), prop_var_noise))) +
  scale_fill_viridis_d(end = 0.8, option = "B", direction = -1, guide = "none") +
  scale_color_viridis_d(name = "Ind. Var. (%)", end = 0.8, option = "B", direction = -1,
                        guide = guide_legend(override.aes = list(size = 2))) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size = 9),
        title = element_text(size = 7),
        strip.text.y = element_text(vjust = -3),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.99, 0.99),
        legend.justification = c("right", "top"),
        legend.direction = "horizontal",
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.key.size = unit(10, 'pt'),
        legend.background = element_blank(),
        plot.margin = unit(c(0.1,0.1,-0.31,-0.5), "cm")
  )

# colorblindr::cvd_grid()


################################################################################
## Figure : normalisation and transformation
################################################################################

df <- subset(result.table_plot,
             effect_size == base_effect_size
             & prop_var_tree == base_prop_var_tree
             & cond_type == base_cond_type
             & model_process == "BM"
             & fact_disp == base_fact_disp
             & score %in% c("TPR")
             & methodMod %in% c("phylolm_BM"))
df <- subset(df, use_lengths ==  "with_lengths" & tree_type == "real_tree")

df$lengthNormalization <- factor(df$lengthNormalization, levels = c("none", "RPKM", "TPM"))

if (nrow(df) / 1 / 2 / 3 != Nrep) message("wrong dimension")
if (anyNA(df$score_value)) message("Some NAs in the data")

ggplot(df, aes(x = lengthNormalization, y = score_value, fill = transformation, color = transformation)) +
  facet_grid(score ~ ., scales = "free", switch = "both", labeller = labeller(Design = label_both)) +
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, show.legend = FALSE,
               position = "identity", width = 0.2*3/5,
               aes(group = interaction(effect_size, methodMod, Design, transformation, lengthNormalization, prop_var_tree))) +
  stat_summary(fun = median, geom = 'line', aes(group = interaction(effect_size, methodMod, Design, transformation, prop_var_tree))) +
  scale_fill_viridis_d(name = "Transformation", end =  0.8, begin = 0.4, option = "B", guide = "none") +
  scale_color_viridis_d(name = "Transformation", end = 0.8, begin = 0.4, option = "B", guide = guide_legend(override.aes = list(size = 2))) +
  guides(alpha = guide_legend(override.aes = list(fill = "black"))) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size = 9),
        title = element_text(size = 7),
        strip.text.y = element_text(vjust = -3),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.99, 0.01),
        legend.justification = c("right", "bottom"),
        legend.key.size = unit(10, 'pt'),
        plot.margin = unit(c(0.1,0.1,-0.31,-0.5), "cm")
  )

# colorblindr::cvd_grid()


################################################################################
## Sample Dataset
################################################################################

library(compcodeR)
tree_type <- "real_tree"
model_process <- "BM"
use_lengths <- "with_lengths"
effect_size <- base_effect_size
prop_var_tree <- base_prop_var_tree
fact_disp <- base_fact_disp
i <- 1

## Sight
cond_type <- "sights"
dataset <- paste(tree_type, cond_type, model_process, use_lengths, effect_size, prop_var_tree, fact_disp, sep = "_")
data_set_sights <- readRDS(file.path(simus_directory, paste0(dataset, "_", i, ".rds")))
# tree
tree <- data_set_sights@tree
# params
nGenes <- nrow(data_set_sights@count.matrix)
# Cond
conds_sights <- data.frame(label = rownames(data_set_sights@sample.annotations),
                           condition = as.factor(data_set_sights@sample.annotations$id.condition))
conds_sights$Design <- "sight"
# DE
data_set_sights <- data_set_sights@variable.annotations
data_set_sights$differential.expression <- as.factor(data_set_sights$differential.expression)
data_set_sights$Design <- "sight"

## bloc
cond_type <- "bloc"
dataset <- paste(tree_type, cond_type, model_process, use_lengths, effect_size, prop_var_tree, fact_disp, sep = "_")
data_set_bloc <- readRDS(file.path(simus_directory, paste0(dataset, "_", i, ".rds")))
# Cond
conds_bloc <- data.frame(label = rownames(data_set_bloc@sample.annotations),
                         condition = as.factor(data_set_bloc@sample.annotations$id.condition))
conds_bloc$Design <- "block"
# DE
data_set_bloc <- data_set_bloc@variable.annotations
data_set_bloc$differential.expression <- as.factor(data_set_bloc$differential.expression)
data_set_bloc$Design <- "block"

## alt
cond_type <- "alt"
dataset <- paste(tree_type, cond_type, model_process, use_lengths, effect_size, prop_var_tree, fact_disp, sep = "_")
data_set_alt <- readRDS(file.path(simus_directory, paste0(dataset, "_", i, ".rds")))
# Cond
conds_alt <- data.frame(label = rownames(data_set_alt@sample.annotations),
                        condition = as.factor(data_set_alt@sample.annotations$id.condition))
conds_alt$Design <- cond_type
# DE
data_set_alt <- data_set_alt@variable.annotations
data_set_alt$differential.expression <- as.factor(data_set_alt$differential.expression)
data_set_alt$Design <- cond_type

data_set <- rbind(#data_set_sights, 
  data_set_alt,
  data_set_bloc)
data_set$Design <- factor(data_set$Design, c("block", "sight", "alt"))

conds <- rbind(conds_bloc, conds_sights, conds_alt)

rm(tree_type, cond_type, model_process, use_lengths, effect_size, prop_var_tree, fact_disp, i)

################################################################################
## Plot Sample Dataset
################################################################################

ggplot(data_set,
       aes(x = A.value, y = M.value))+#, color = differential.expression)) + #, alpha = differential.expression)) + 
  stat_density2d(data = subset(data_set, differential.expression == 0),
                 aes(fill = ..count..), geom = "tile", contour = FALSE, n = 50) + 
  geom_point(data = subset(data_set, differential.expression == 1), color = "#D55E00", size = 0.5) +
  facet_grid(.~Design) + 
  scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256), limits = c(0, NA)) + 
  xlab("M value") +
  ylab("A value") +
  coord_cartesian(xlim = c(2, 14)) + 
  theme_bw() +
  theme(text = element_text(size = 9),
        title = element_text(size = 7),
        strip.text.x = element_text(vjust = -2),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(-0.1, -0.22),
        legend.justification = c("left", "top"),
        legend.key.height = unit(3, "pt"),
        legend.margin = margin(-0.1,0,0,0, "cm"),
        legend.direction = "horizontal",
        legend.key.size = unit(10, 'pt'),
        plot.margin = unit(c(-0.3,0,0.1,0), "cm")
  ) +
  guides(fill = guide_colorbar(label.position = "bottom",
                               title.position = "left",
                               title.vjust = 1.1,
                               label.vjust = 2.5))
# colorblindr::cvd_grid()

################################################################################
## Plot Tree Design
################################################################################
library(ggtree)
library(tidytree)
library(aplot)
library(cowplot)

treedata <- as_tibble(tree)
treedata <- full_join(treedata, conds, by = 'label')
treedata <- as.treedata(treedata)

## Tree
gt <- ggtree(tree)
gt <- flip(gt, 64, 62)
gt <- flip(gt, 60, 58)
gt <- flip(gt, 8, 42)
gt <- gt + geom_tippoint(size = 1)
gt <- gt + theme(legend.position = c(0.1, 0.8),
                 text = element_text(size = 9),
                 plot.margin = unit(c(0,0,0,0), "cm"))

## Data
d <- filter(gt, isTip) %>% select(c(label, y))
conds_plot <- left_join(conds, d, by='label')
conds_plot$label <- gsub("_.*", "", conds_plot$label)

aa <- subset(conds_plot, Design == "alt")
aa$ystart <- aa$y - 0.5
aa$yend <- aa$y + 0.5

conds_plot$Design <- factor(conds_plot$Design, c("block", "alt", "sight"))
conds_plot$real <- conds_plot$Design == "sight"
p1 <- ggplot(subset(conds_plot, real), aes(x = Design, y = y, col = condition)) +
  geom_point(size = 1.6, shape = 15) + 
  scale_color_viridis_d(end = 0.7, option = "D") + 
  theme_tree2() + ylim2(gt) + 
  theme(legend.position='none',
        text = element_text(size = 9),
        strip.placement = "outside",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))

ll <- subset(conds_plot, Design == "sight")
ll$ypos <- sapply(ll$label, function(x) sum(range(subset(ll, label == x)$y)) / 2)
names_species <- c(CTENE = "Cambarus \n tenebrosus",
                   CRUST = "Cambarus \n rusticiformis",
                   CGRAY = "Cambarus \n graysoni",
                   OINCO = "Orconectes \n incomptus",
                   OAUST = "Orconectes \n australis",
                   CNERT = "Cambarus \n nerterius",
                   CHAMU = "Cambarus \n hamulatus",
                   CCRYP = "C. cryptodytes",
                   CDUBI = "Cambarus \n dubius",
                   CSETO = "Cambarus \n setosus",
                   PPALL = "P. pallidus",
                   PHORS = "Procambarus \n horsti",
                   PLUCI = "Procambarus \n lucifugus",
                   PFALL = "Procambarus \n fallax"
)
ll$label_full <- names_species[ll$label]
ll$oneline <- sapply(ll$label, function(x) x %in% c("CCRYP", "PPALL"))
pLabel <- ggplot(ll, aes(x = 0, y = ypos, col = condition)) + 
  geom_text(data = subset(ll, !oneline), 
            aes(label = label_full, hjust = 0, vjust = 0.3), size = 1.8, fontface = "italic", lineheight = 0.5) +
  geom_text(data = subset(ll, oneline), 
            aes(label = label_full, hjust = 0, vjust = 0), size = 1.8, fontface = "italic", lineheight = 0.5) +
  xlim(0, 0.02) +
  theme_tree2() + ylim2(gt) + 
  theme(legend.position='none') +
  scale_color_viridis_d(end = 0.7, option = "D") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))

p2 <- ggplot(subset(conds_plot, !real), aes(x = Design, y = y, col = condition)) +
  geom_point(size = 1.6, shape = 15, alpha = 0.6) + 
  scale_color_viridis_d(end = 0.7, option = "A") + 
  theme_tree2() + ylim2(gt) + 
  theme(legend.position='none',
        # plot.title = element_text(size = 10),
        text = element_text(size = 9),
        strip.placement = "outside",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))

## Grid
pp <- plot_grid(gt, p2, p1, pLabel, ncol = 4, align = 'h', axis = "b", rel_widths = c(1, 0.2, 0.1, 0.3))
## background
pb <- ggplot(aa) + 
  geom_rect(aes(xmin = 0, xmax = 1,
                ymin = ystart, ymax = yend,
                fill = condition,
                col = NULL), alpha = 0.1) +
  ylim2(gt) + xlab("") + 
  scale_x_continuous(breaks = c(0)) + 
  scale_fill_grey(guide = "none") + 
  theme_tree2() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "white", size = 9),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
aligned_plots <- align_plots(pp, pb, align="hv", axis="tlr")
ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
# colorblindr::cvd_grid(ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]]))

###############################################################################
## Figure : SVA nsv
################################################################################
df <- subset(result.table,
             lengthNormalization %in% c("TPM", "length")
             & transformation %in% c("log2", "none")
             & effect_size == base_effect_size
             & fact_disp == base_fact_disp
             & prop_var_tree == base_prop_var_tree
             & model_process == "BM"
             # & score %in% c(plot_scores, "nsv")
             & methodMod %in% c(
               "limma_sva_auto")
)
df <- subset(df, use_lengths ==  "with_lengths" & tree_type == "real_tree")
df$methodMod <- factor(df$methodMod, levels = c("DESeq2", "limma", "limma_cor", "limma_sva_one", "limma_sva_auto", "phylolm_BM",
                                                "phylolm_OU"
))
levels(df$methodMod) <- pretty_names[levels(df$methodMod)]

if (nrow(df) / 3 / 1 / 1 != Nrep) message("wrong dimension")
if (anyNA(df$score_value)) message("Some NAs in the data")

dodge <- position_dodge(width = 0)

ggplot(df, aes(x = Design, y = nsv, fill = Design, color = Design)) +
  # geom_hline(aes(yintercept = threshold), linetype = 2) +
  # facet_grid(nsv ~ ., scales = "free", switch = "y", labeller = labeller(Design = label_both)) +
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, show.legend = FALSE,
               position = "identity", width = 0.2,
               aes(group = interaction(effect_size, use_lengths, tree_type, methodMod, prop_var_tree, Design))) +
  stat_summary(fun = median, geom = 'line', aes(group = interaction(effect_size, use_lengths, tree_type, prop_var_tree), colour = Design)) +
  scale_fill_viridis_d(guide = "none", end = 0.8, option = "B") +
  scale_colour_viridis_d(end = 0.8, option = "B", guide ="none") +
  ylab("Number of Surogate Variables") +
  # xlab("") +
  theme_bw() +
  theme(text = element_text(size = 9),
        title = element_text(size = 7),
        # strip.text.y = element_text(vjust = -3),
        # strip.placement = "outside",
        # strip.background = element_blank(),
        # panel.grid.minor = element_blank(),
        legend.position = c(0.01, 0.01),
        legend.justification = c("left", "bottom"),
        legend.key.size = unit(10, 'pt'),
        # plot.margin = unit(c(0.1,0.1,-0.31,-0.5), "cm")
  )
# colorblindr::cvd_grid

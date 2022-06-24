WDir="Your working directory"
setwd(WDir)


library(phyloseq)
library(dplyr)
library(ggplot2)
library(ggvenn)
library(patchwork)
library(microbiomeMarker)

# rm(list = ls())
# 
# 
# fdata16S.genus = readRDS("./RDSdata/fdata16S_genus.rds")
# fdata16S.genus
# 
# fdataITS.genus = readRDS("./RDSdata/fdataITS_genus.rds")
# fdataITS.genus
# 
# source("../tools/visual_tools.R")
# 
# set.seed(110)
# f16.lefse.genus = run_lefse(fdata16S.genus, group = "metamStatus", taxa_rank = "Genus")
# saveRDS(f16.lefse.genus, "./Field_analysis_output_v3/Lefse/f16_lfse_genus.rds")
# 
# 
# set.seed(102)
# fITS.lefse.genus = run_lefse(fdataITS.genus, group = "metamStatus", taxa_rank = "Genus")
# saveRDS(fITS.lefse.genus, "./Field_analysis_output_v3/Lefse/fITS_lfse_genus.rds")

#-----------------------------------------------------------------------------------------
rm(list = ls())

fdata16S.genus = readRDS("./RDSdata/fdata16S_genus.rds")
fdataITS.genus = readRDS("./RDSdata/fdataITS_genus.rds")
f16.lefse.genus = readRDS('./Field_analysis_output_v3/Lefse/f16_lfse_genus.rds')
fITS.lefse.genus = readRDS("./Field_analysis_output_v3/Lefse/fITS_lfse_genus.rds")

source("../tools/visual_tools.R")
# plot
(plot.16s = plot_lefse(f16.lefse.genus, phy = fdata16.genus, colorSet = ggsci::pal_jama()(2), label_print = 0, ylab = "Genus"))

#saveRDS(plot.16s, "./Field_analysis_output_v3/Lefse/plot_16S_metam_lefse_italic.rds")

(plot.ITS = plot_lefse(fITS.lefse.genus, phy = fdataITS.genus, colorSet = ggsci::pal_jama()(2), label_print = 0, ylab = "Genus"))

#saveRDS(plot.ITS, "./Field_analysis_output_v3/Lefse/plot_ITS_metam_lefse_italic.rds")


# combine plots figure 2
rm(list = ls())

plot.16s = readRDS("./Field_analysis_output_v3/Lefse/plot_16S_metam_lefse_italic.rds")
plot.ITS = readRDS("./Field_analysis_output_v3/Lefse/plot_ITS_metam_lefse_italic.rds")
comp.16S = readRDS("./Field_analysis_output_v3/comps/plot_f16_comp_metam_italic.rds")
comp.ITS = readRDS("./Field_analysis_output_v3/comps/plot_fITS_comp_metam_italic.rds")

comp.16S + plot.16s$graph + comp.ITS + plot.ITS$graph  + plot_layout(widths = c(1,1)) + plot_annotation(tag_levels = "A")


length(unique(plot.16s$df[plot.16s$df$enrich_group == "Not fumigated",]$feature)) #122
length(unique(plot.16s$df[plot.16s$df$enrich_group == "Fumigated",]$feature)) #110

length(unique(plot.ITS$df[plot.ITS$df$enrich_group == "Not fumigated",]$feature)) #53
length(unique(plot.ITS$df[plot.ITS$df$enrich_group == "Fumigated",]$feature)) #39

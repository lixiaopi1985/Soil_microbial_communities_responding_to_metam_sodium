WDir="Your working directory"
setwd(WDir)
library(phyloseq)
library(tidyverse)
library(tools)
library(vegan)
library(ampvis2)
library(microbiome)
library(microbiomeutilities)
library(ggrepel)
library(vegan)
library(ggplot2)
library(tidyverse)
library(ggvegan)
library(ggrepel)
library(ggthemes)
library(RColorBrewer)
library(stringr)
library(latex2exp)
library(pals)
library(lmerTest)
library(car)
library(indicspecies)
library(SpiecEasi)
library(igraph)
library(DESeq2)
library(emmeans)
library(pheatmap)
library(ggpubr)
library(FSA)
library(rcompanion)
library(patchwork)


rm(list = ls())

fdata_16s.prune = readRDS("../Data_prunning/RDS_Field_data_prunning/fdata16S_cleantaxa_bst2.rds")
fdata_its.prune = readRDS("../Data_prunning/RDS_Field_data_prunning/fdataITS_cleantaxa_bst2.rds")

# genus level data frame

# fdata16.genus = fdata_16s.prune %>%
#   tax_glom("Genus") %>%
#   psmelt()

# fdata16.species = fdata_16s.prune %>%
#   tax_glom("Species") %>%
#   psmelt()

# fdataITS.genus = fdata_its.prune %>%
#   tax_glom("Genus") %>%
#   psmelt()

# saveRDS(fdata16.genus, "./RDSdata/fdata16S_genus_df.rds")
# saveRDS(fdata16.species, "./RDSdata/fdata16S_species_df.rds")
# saveRDS(fdataITS.genus, "./RDSdata/fdataITS_genus_df.rds")


##############################################################################

rm(list = ls())

fdata_16s.prune = readRDS("../Data_prunning/RDS_Field_data_prunning/fdata16S_cleantaxa_bst2.rds")
fdata_its.prune = readRDS("../Data_prunning/RDS_Field_data_prunning/fdataITS_cleantaxa_bst2.rds")
fdata16.genus = readRDS("./RDSdata/fdata16S_genus_df.rds")
fdataITS.genus = readRDS("./RDSdata/fdataITS_genus_df.rds")

fdata16.genus$field_alias = factor(fdata16.genus$field_alias, levels = unique(fdata16.genus$field_alias))
fdataITS.genus$field_alias = factor(fdataITS.genus$field_alias, levels = unique(fdataITS.genus$field_alias))




source("../tools/visual_tools.R")


(f16.comp = plotBars(fdata16.genus, groupby = c("field_alias", "metamStatus"), plotx="field_alias", facet_group = "~ metamStatus", ylab = "Relative abundance (%)", nshow = 10, choosecolorset = ggthemes::calc_pal()(10), otherColor = "grey45", printlevels = 0, xlab = "Fields", plotMargin_r = 0, fontHjust = 0.5))


(fITS.comp=plotBars(fdataITS.genus, groupby = c("field_alias", "metamStatus"), plotx="field_alias", facet_group = "~ metamStatus", ylab = "Relative abundance (%)", nshow = 10, choosecolorset = ggsci::pal_igv()(10), otherColor = "grey45", printlevels = 0, xlab = "Fields", plotMargin_r = 20, fontHjust = 0.5))

f16.comp / fITS.comp + plot_annotation(tag_levels = "A")


(f16.comp.metam = plotBars(fdata16.genus, groupby = c("metamStatus"), plotx="metamStatus", facet_group = "", ylab = "Relative abundance (%)", nshow = 10, choosecolorset = ggthemes::calc_pal()(10), otherColor = "grey45", printlevels = 0, xlab = "Fumigation history", plotMargin_r = 0, fontHjust = 0.5))

(fITS.comp.metam=plotBars(fdataITS.genus, groupby = c("metamStatus"), plotx="metamStatus", facet_group = "", ylab = "Relative abundance (%)", nshow = 10, choosecolorset = ggsci::pal_igv()(10), otherColor = "grey45", printlevels = 0, xlab = "Fumigation history", plotMargin_r = 40, fontHjust = 0.5))


saveRDS(f16.comp.metam, "./Field_analysis_output_v3/comps/plot_f16_comp_metam_italic.rds")
saveRDS(fITS.comp.metam, "./Field_analysis_output_v3/comps/plot_fITS_comp_metam_italic.rds")

#-----------------------
# summary data
#-----------------------

fdata16.genus %>%
  mutate(total_abund = sum(Abundance)) %>%
  group_by(Genus) %>%
  mutate(rel_abund = 100*Abundance/total_abund) %>%
  summarise(acc_rel_abund = sum(rel_abund)) %>%
  arrange(desc(acc_rel_abund))

fdataITS.genus %>%
  mutate(total_abund = sum(Abundance)) %>%
  group_by(Genus) %>%
  mutate(rel_abund = 100*Abundance/total_abund) %>%
  summarise(acc_rel_abund = sum(rel_abund)) %>%
  arrange(desc(acc_rel_abund))


# by metam sodium
fdata16.genus %>%
  group_by(metamStatus) %>%
  mutate(total_abund = sum(Abundance)) %>%
  group_by(metamStatus, Genus) %>%
  mutate(rel_abund = 100*sum(Abundance)/total_abund) %>%
  summarize(sum_abund = mean(rel_abund)) %>%
  slice_max(order_by = sum_abund, n=10) 

fdataITS.genus %>%
  group_by(metamStatus) %>%
  mutate(total_abund = sum(Abundance)) %>%
  group_by(metamStatus, Genus) %>%
  mutate(rel_abund = 100*sum(Abundance)/total_abund) %>%
  summarize(sum_abund = round(mean(rel_abund), 2)) %>%
  slice_max(order_by = sum_abund, n=10)

Wdir="Your working directory"
setwd(Wdir)


library(phyloseq)
library(corncob)
library(dplyr)
library(ggpubr)
library(patchwork)


rm(list = ls())

# load data
m16.gen = readRDS("./RDSdata/mdata16S_bst_genus.rds")
mITS.gen = readRDS("./RDSdata/mdataITS_bst_genus.rds")

m16.gen.metam = subset_samples(m16.gen, metamStatus == "Fumigated")
m16.gen.NOmetam = subset_samples(m16.gen, metamStatus == "Not fumigated")

mITS.gen.metam =   subset_samples(mITS.gen, metamStatus == "Fumigated")
mITS.gen.NOmetam = subset_samples(mITS.gen, metamStatus == "Not fumigated")


m16.gen.df = readRDS("./RDSdata/mdata16S_bst_genus_df.rds")
mITS.gen.df = readRDS("./RDSdata/mdataITS_bst_genus_df.rds")

# model information
# corn16S.y = readRDS("./MCM_analysis_output_v3/corncob/metam16S_TOP10_period_nperm1000_DA_abundance_only.rds")
# corn16S.n = readRDS("./MCM_analysis_output_v3/corncob/NOmetam16S_TOP10_period_nperm1000_DA_abundance_only.rds")

corn16S.y = readRDS("./MCM_analysis_output_v3/corncob/metam16S_TOP10_period_nperm10_000_DA_abundance_only.rds")
corn16S.n = readRDS("./MCM_analysis_output_v3/corncob/NOmetam16S_TOP10_period_nperm10_000_DA_abundance_only.rds")

# cornITS.y = readRDS("./MCM_analysis_output_v3/corncob/metamITS_period_nperm1000_DA_abundance_only.rds")
# cornITS.n = readRDS("./MCM_analysis_output_v3/corncob/NOmetamITS_period_nperm1000_DA_abundance_only.rds")

cornITS.y = readRDS("./MCM_analysis_output_v3/corncob/metamITS_period_nperm10_000_DA_abundance_only.rds")
cornITS.n = readRDS("./MCM_analysis_output_v3/corncob/NOmetamITS_period_nperm10_000_DA_abundance_only.rds")

source("../tools/plotHeatMap2.R")
source("../tools/visual_tools.R")


top16S = plotHeat_matrix(m16.gen.df, groupby = c("period", "metamStatus"))
topITS = plotHeat_matrix(mITS.gen.df, groupby = c("period", "metamStatus"))

top16S$m
top16S$overal
top16S$dom



modelCorn.16S.y = plotCoef_mcm(corn16S.y, m16.gen.metam, orderList = top16S$dom)

modelCorn.16S.y$DF
modelCorn.16S.y$mod

modelCorn.16S.n = plotCoef_mcm(corn16S.n, m16.gen.NOmetam, orderList = top16S$dom)

modelCorn.16S.n$DF
modelCorn.16S.n$mod

modelCorn.ITS.y = plotCoef_mcm(cornITS.y, mITS.gen.metam, orderList = topITS$dom)
modelCorn.ITS.y$DF
modelCorn.ITS.y$mod


modelCorn.ITS.n = plotCoef_mcm(cornITS.n, mITS.gen.NOmetam, orderList = topITS$dom)

modelCorn.ITS.n$DF
modelCorn.ITS.n$mod

nbrk = seq(min(top16S$m), max(topITS$m), by=0.005)
nbrk

corn16.y.sigTaxa = as.data.frame(tax_table(corn16S.y$data))[corn16S.y$significant_taxa,"Genus"]
corn16.n.sigTaxa = as.data.frame(tax_table(corn16S.n$data))[corn16S.n$significant_taxa,"Genus"]
cornITS.y.sigTaxa = as.data.frame(tax_table(cornITS.y$data))[cornITS.y$significant_taxa,"Genus"]
cornITS.n.sigTaxa = as.data.frame(tax_table(cornITS.n$data))[cornITS.n$significant_taxa,"Genus"]

# add significant symbols
source("../tools/visual_tools.R")
new_row_label_16S = addSig(top16S$dom,  corn16.n.sigTaxa, corn16.y.sigTaxa)
new_row_label_16S
new_row_label_ITs = addSig(topITS$dom, cornITS.n.sigTaxa, cornITS.y.sigTaxa)
new_row_label_ITs
newcol_labels =rep(c("microcosm T0", "microcosm T1", "microcosm T2", "microcosm T3"),2)
newcol_labels

source("../tools/plotHeatMap2.R")
source("../tools/visual_tools.R")

newlabs16S = scientific_name_formatter(new_row_label_16S, plotdevice = "pheat")
newlabsITS = scientific_name_formatter(new_row_label_ITs, plotdevice = "pheat")

(p16S = plotHeat_plot(top16S, new_break = nbrk, out = T, outfile = "./MCM_analysis_output_v3/heatmap/hp_16S_metam_sig_italic.png", cellw = 32, cellh = 32, labels_col = newcol_labels, show_rownames=T, labels_row = as.expression(newlabs16S)))

(pITS = plotHeat_plot(topITS, new_break = nbrk, out = T,outfile = "./MCM_analysis_output_v3/heatmap/hp_ITS_metam_sig_italic.png", cellw = 32, cellh = 32, labels_col = newcol_labels, labels_row = as.expression(newlabsITS)))

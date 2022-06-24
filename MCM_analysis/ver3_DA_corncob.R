Wdir="Your working directory"
setwd(Wdir)

library(phyloseq)
library(corncob)
library(tidyr)
library(dplyr)

rm(list = ls())

m16.gen =  readRDS("./RDSdata/mdata16S_bst_genus.rds")
mITS.gen = readRDS("./RDSdata/mdataITS_bst_genus.rds")

# m16.gen.df = m16.gen %>%
#   psmelt()
# 
# mITS.gen.df = mITS.gen %>%
#   psmelt()

# saveRDS(m16.gen.df, "./RDSdata/mdata16S_bst_genus_df.rds")
# saveRDS(mITS.gen.df, "./RDSdata/mdataITS_bst_genus_df.rds")

m16.gen.df = readRDS("./RDSdata/mdata16S_bst_genus_df.rds")
mITS.gen.df = readRDS("./RDSdata/mdataITS_bst_genus_df.rds")


source("../tools/plotHeatMap2.R")


mat16S = plotHeat_matrix(m16.gen.df, groupby = c("period", "metamStatus"))

mat16S$dom
mat16S$overal

# # subset taxa in phyloseq
top16S = subset_taxa(m16.gen, Genus %in% mat16S$dom)
meta16S = data.frame(sample_data(m16.gen))

meta16S$metamStatus

fum16S.y = subset_samples(m16.gen, metamStatus == "Fumigated")
fum16S.n = subset_samples(m16.gen, metamStatus == "Not fumigated")


set.seed(111)
nboots = 10000
corn16S.y = differentialTest(formula = ~ period,
                             formula_null = ~ 1,
                             phi.formula = ~ 1, # test on variability
                             phi.formula_null = ~ 1,
                             test = "LRT",
                             boot = T,
                             B = nboots,
                             data = fum16S.y,
                             full_output = T,
                             try_only = taxa_names(top16S)
)

# corn16S.y$significant_taxa
# plot(corn16S.y, level="Genus")

# saveRDS(corn16S.y, "./MCM_analysis_output_v3/corncob/metam16S_TOP10_period_nperm1000_DA_abundance_only.rds")
saveRDS(corn16S.y, "./MCM_analysis_output_v3/corncob/metam16S_TOP10_period_nperm10_000_DA_abundance_only.rds")

corn16S.n = differentialTest(formula = ~ period,
                             formula_null = ~ 1,
                             phi.formula = ~ 1,
                             phi.formula_null = ~ 1,
                             test = "LRT",
                             boot = T,
                             B = nboots,
                             full_output = T,
                             data = fum16S.n,
                             try_only = taxa_names(top16S)
)

# saveRDS(corn16S.n, "./MCM_analysis_output_v3/corncob/NOmetam16S_TOP10_period_nperm1000_DA_abundance_only.rds")

saveRDS(corn16S.n, "./MCM_analysis_output_v3/corncob/NOmetam16S_TOP10_period_nperm10_000_DA_abundance_only.rds")

# corn16S.n$significant_taxa
# plot(corn16S.n, level="Genus")

#-----------------------------------------------------------------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ fungi ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
matITS = plotHeat_matrix(mITS.gen.df, groupby = c("period", "metamStatus"))

matITS$dom
matITS$overal


# subset taxa in phyloseq
topITS = subset_taxa(mITS.gen, Genus %in% matITS$dom)

fum.y.ITS = subset_samples(topITS, metamStatus == "Fumigated")
fum.n.ITS = subset_samples(topITS, metamStatus == "Not fumigated")

cornITS.y = differentialTest(formula = ~ period,
                             formula_null = ~ 1,
                             phi.formula = ~ 1,
                             phi.formula_null = ~1,
                             test = "LRT",
                             boot = T,
                             data = fum.y.ITS,
                             B = nboots,
                             full_output = T,
                             try_only = taxa_names(topITS)
)

#saveRDS(cornITS.y, "./MCM_analysis_output_v3/corncob/metamITS_period_nperm1000_DA_abundance_only.rds")
saveRDS(cornITS.y, "./MCM_analysis_output_v3/corncob/metamITS_period_nperm10_000_DA_abundance_only.rds")

cornITS.n = differentialTest(formula = ~ period,
                             formula_null = ~ 1,
                             phi.formula = ~ 1,
                             phi.formula_null = ~1,
                             test = "LRT",
                             boot = T,
                             B = nboots,
                             full_output = T,
                             try_only = taxa_names(topITS),
                             data = fum.n.ITS
)
#saveRDS(cornITS.n, "./MCM_analysis_output_v3/corncob/NOmetamITS_period_nperm1000_DA_abundance_only.rds")
saveRDS(cornITS.n, "./MCM_analysis_output_v3/corncob/NOmetamITS_period_nperm10_000_DA_abundance_only.rds")

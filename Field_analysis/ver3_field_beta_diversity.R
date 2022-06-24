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


# --------------------------------
# glom to genus level
# --------------------------------

fdata16S = readRDS("../Data_prunning/RDS_Field_data_prunning/fdata16S_cleantaxa_bst2.rds")
fdataITS = readRDS("../Data_prunning/RDS_Field_data_prunning/fdataITS_cleantaxa_bst2.rds")

# fdata16S.genus = fdata16S %>%
#   tax_glom("Genus")
# saveRDS(fdata16S.genus, "./RDSdata/fdata16S_genus.rds")
# 
# fdataITS.genus = fdataITS %>%
#   tax_glom("Genus")
# saveRDS(fdataITS.genus, "./RDSdata/fdataITS_genus.rds")
fdata16S.genus = readRDS("./RDSdata/fdata16S_genus.rds")
fdataITS.genus = readRDS("./RDSdata/fdataITS_genus.rds")

norm16S = fdata16S.genus %>%
  transform_sample_counts(function(x)x/sum(x))

normITS = fdataITS.genus %>%
  transform_sample_counts(function(x)x/sum(x))

meta16S = data.frame(sample_data(fdata16S.genus))
metaITS = data.frame(sample_data(fdataITS.genus))

#-----------------------------------------------------------------------------------
otu.16s = data.frame(t(otu_table(norm16S)))
env.16s = data.frame(sample_data(norm16S))

cap.16s = capscale(otu.16s ~ metamStatus+soil_compname+totalCrops+ph, data = env.16s, distance = "bray", add = T)
vif.cca(cap.16s)

# metamStatusFumigated    soil_compnameQuincy  soil_compnameSagehill   soil_compnameTaunton 
# 1.564482               1.512294               1.429474               1.594264 
# soil_compnameTimmerman             totalCrops                     ph 
# 4.765879               1.373700               3.846869 

# check this model
set.seed(123)
anova(cap.16s, permutations = 1000) # 0.001
anova(cap.16s, by="terms", permutations = 1000) # ph is not significant
anova(cap.16s, by="axis", permutations = 1000) # first 3

# get parsimony model
# by AIC
# ?vegan::ordistep
cap.16s.select = ordistep(capscale(otu.16s ~ 1, data = env.16s, distance = "bray", add=T), scope = formula(cap.16s), direction = "forward", permutations = 10000)

# adjust p value
cap.16s.select.padj = cap.16s.select
cap.16s.select.padj$anova$`Pr(>F)` = p.adjust(cap.16s.select$anova$`Pr(>F)`, method = "BH")
cap.16s.select.padj$anova
cap.16s.select.padj

cap.16s.select


# check this model
set.seed(123)
anova(cap.16s, permutations = 1000) # 0.001
anova(cap.16s, by="terms", permutations = 1000) # ph is not significant
anova(cap.16s, by="axis", permutations = 1000) # first 3


cap.16s.select = ordistep(capscale(otu.16s ~ 1, data = env.16s, distance = "bray", add=T), scope = formula(cap.16s), direction = "forward", permutations = 10000)

# adjust p value
cap.16s.select.padj = cap.16s.select
cap.16s.select.padj$anova$`Pr(>F)` = p.adjust(cap.16s.select$anova$`Pr(>F)`, method = "BH")
cap.16s.select.padj$anova
cap.16s.select.padj
cap.16s.select

set.seed(111)
cap.16s.final = capscale(otu.16s ~ soil_compname + metamStatus + totalCrops, 
                         data = env.16s, distance="bray", add=T)


plot(cap.16s.final)


set.seed(104)
anova(cap.16s.final, permutations = 10000, step = 1000) # 0.0000999
anova(cap.16s.final, by="terms", permutations = 10000, step = 1000) # all sig

# Permutation test for capscale under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 10000
# 
# Model: capscale(formula = otu.16s ~ soil_compname + metamStatus + totalCrops, data = env.16s, distance = "bray", add = T)
# Df SumOfSqs      F     Pr(>F)    
# soil_compname  4  1.17767 5.2679 0.00009999 ***
#   metamStatus    1  0.38404 6.8715 0.00009999 ***
#   totalCrops     1  0.21403 3.8295     0.0005 ***
#   Residual      41  2.29144                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


anova(cap.16s.final, by="axis", permutations = 10000, step = 1000) # top 3

# Permutation test for capscale under reduced model
# Forward tests for axes
# Permutation: free
# Number of permutations: 10000
# 
# Model: capscale(formula = otu.16s ~ soil_compname + metamStatus + totalCrops, data = env.16s, distance = "bray", add = T)
# Df SumOfSqs       F     Pr(>F)    
# CAP1      1  0.84162 15.0588 0.00009999 ***
#   CAP2      1  0.53766  9.6202 0.00009999 ***
#   CAP3      1  0.19312  3.4554     0.0029 ** 
#   CAP4      1  0.11277  2.0177     0.1096    
# CAP5      1  0.05780  1.0342     0.7132    
# CAP6      1  0.03277  0.5864     0.9441    
# Residual 41  2.29144                       
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#-----------------------------------------------------------------------------
otu.its = data.frame(t(otu_table(normITS)))
env.its = data.frame(sample_data(normITS))

cap.its = capscale(otu.its ~ metamStatus+soil_compname+totalCrops+ph, data = env.its, distance = "bray", add = T)
cap.its
vif.cca(cap.its) # remove
# metamStatusFumigated    soil_compnameQuincy  soil_compnameSagehill   soil_compnameTaunton 
# 1.564482               1.512294               1.429474               1.594264 
# soil_compnameTimmerman             totalCrops                     ph 
# 4.765879               1.373700               3.846869 

# check this model
set.seed(123)
anova(cap.its, permutations = 1000) # 0.001
anova(cap.its, by="terms", permutations = 1000) # ph is not significant
anova(cap.its, by="axis", permutations = 1000) # first 1

# get parsimony model
# by AIC
cap.its.select = ordistep(capscale(otu.its ~ 1, data = env.its, distance = "bray", add=T), scope = formula(cap.its), direction = "forward", permutations = 10000)

# adjust p value
cap.its.select.padj = cap.its.select
cap.its.select.padj$anova$`Pr(>F)` = p.adjust(cap.its.select$anova$`Pr(>F)`, method = "BH")
cap.its.select.padj$anova
cap.its.select.padj

cap.its.select

set.seed(111)
cap.its.final = capscale(otu.its ~ metamStatus + soil_compname + totalCrops, data = env.its, distance="bray", add=T)

plot(cap.its.final)

set.seed(104)
anova(cap.its.final, permutations = 10000, step = 1000) # 0.0000999
anova(cap.its.final, by="terms", permutations = 10000, step = 1000) # all sig

# Permutation test for capscale under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 10000
# 
# Model: capscale(formula = otu.its ~ metamStatus + soil_compname + totalCrops, data = env.its, distance = "bray", add = T)
# Df SumOfSqs      F     Pr(>F)    
# metamStatus    1   0.7092 3.9101 0.00009999 ***
#   soil_compname  4   1.2371 1.7051     0.0005 ***
#   totalCrops     1   0.2970 1.6375     0.0339 *  
#   Residual      41   7.4370                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

anova(cap.its.final, by="axis", permutations = 10000, step = 1000) #2

# Permutation test for capscale under reduced model
# Forward tests for axes
# Permutation: free
# Number of permutations: 10000
# 
# Model: capscale(formula = otu.its ~ metamStatus + soil_compname + totalCrops, data = env.its, distance = "bray", add = T)
# Df SumOfSqs      F     Pr(>F)    
# CAP1      1   0.9298 5.1262 0.00009999 ***
#   CAP2      1   0.3977 2.1926     0.1265    
# CAP3      1   0.3541 1.9521     0.1375    
# CAP4      1   0.2226 1.2274     0.6863    
# CAP5      1   0.1820 1.0033     0.7545    
# CAP6      1   0.1571 0.8662     0.6354

###############################################################################
# plot
##################################################
######
# plot
###############################################################################
source("../tools/plotCAP2.R")


taxcols.16s = as.data.frame(tax_table(norm16S))
taxcols.16s$ASV = rownames(taxcols.16s)

# plot
taxcols.its = as.data.frame(tax_table(normITS))
taxcols.its$ASV = rownames(taxcols.its)


(graph16s.cap = plotCAP(cap.16s.final, 
                        env = env.16s, 
                        display = 1, 
                        tax_cols = taxcols.16s, 
                        splitDisplay = F, 
                        sa.xadj = 1.2, 
                        sa.yadj = -0.5, 
                        sa.colorSamples = T, 
                        sa.sample_group = "metamStatus", 
                        sa.showCentroidGroup = "all", 
                        sa.biplot = T,
                        sa.ellipse = F, 
                        sa.addSp = F, 
                        sa.labelsp = F, 
                        sa.pointsize = 6,
                        sa,fontsize = 12,
                        sa.tax_level = "Genus",
                        sa.colorLabs = "Fumigation history"))


(graph16s.cap.soil = plotCAP(cap.16s.final, 
                             env = env.16s, 
                             display = 1, 
                             tax_cols = taxcols.16s, 
                             splitDisplay = F, 
                             sa.xadj = 1.2, 
                             sa.yadj = -0.5, 
                             sa.colorSamples = T, 
                             sa.sample_group = "soil_compname", 
                             sa.showCentroidGroup = "all", 
                             sa.biplot = T,
                             sa.ellipse = F, 
                             sa.addSp = F, 
                             sa.labelsp = F, 
                             sa.pointsize = 6,
                             sa,fontsize = 12,
                             sa.tax_level = "Genus",
                             sa.colorLabs = "Soil series"))


#dir.create("./Field_analysis_ver2/beta_diversity/CAP_RDS", recursive = T)
saveRDS(graph16s.cap, "./Field_analysis_output_v3/beta/CAP_rds/ccametam_16S.rds")
saveRDS(graph16s.cap.soil, "./Field_analysis_output_v3/beta/CAP_rds/ccaSoilseries_16S.rds")


(graphITS.cap = plotCAP(cap.its.final, 
                        env = env.its, 
                        display = 1, 
                        tax_cols = taxcols.its, 
                        splitDisplay = F, 
                        sa.xadj = 0.5, 
                        sa.yadj = 0.1, 
                        sa.colorSamples = T, 
                        sa.sample_group = "metamStatus", 
                        sa.showCentroidGroup = "all", 
                        sa.biplot = T,
                        sa.ellipse = F, 
                        sa.addSp = F, 
                        sa.labelsp = F, 
                        sa.pointsize = 6,
                        sa.fontsize = 12,
                        sa.tax_level = "Genus",
                        sa.colorLabs = "Fumigation history"))


(graphITS.cap.soil = plotCAP(cap.its.final, 
                             env = env.its, 
                             display = 1, 
                             tax_cols = taxcols.its, 
                             splitDisplay = F, 
                             sa.xadj = 0.5, 
                             sa.yadj = 0.1, 
                             sa.colorSamples = T, 
                             sa.sample_group = "soil_compname", 
                             sa.showCentroidGroup = "all", 
                             sa.biplot = T,
                             sa.ellipse = F, 
                             sa.addSp = F, 
                             sa.labelsp = F, 
                             sa.pointsize = 6,
                             sa.fontsize = 12,
                             sa.tax_level = "Genus",
                             sa.colorLabs = "Soil series"))

saveRDS(graphITS.cap, "./Field_analysis_output_v3/beta/CAP_rds/ccametam_ITS.rds")
saveRDS(graphITS.cap.soil, "./Field_analysis_output_v3/beta/CAP_rds/ccaSoilseries_ITS.rds")

#-----------------------------------------------------------------------
#
#
#
# asembly figure 1 and figure S7

rm(list=ls())
graph16S.cap.metam = readRDS("./Field_analysis_output_v3/beta/CAP_rds/ccametam_16S.rds")
graphITS.cap.metam = readRDS("./Field_analysis_output_v3/beta/CAP_rds/ccametam_ITS.rds")

graph16S.cap.soil = readRDS("./Field_analysis_output_v3/beta/CAP_rds/ccaSoilseries_16S.rds")
graphITS.cap.soil = readRDS("./Field_analysis_output_v3/beta/CAP_rds/ccaSoilseries_ITS.rds")


alpha_16S = readRDS("./Field_analysis_output_v3/diversity/plot_obs16S_field.rds")
alpha_ITS = readRDS("./Field_analysis_output_v3/diversity/plot_obsITS_field.rds")

(betaField = 
    alpha_16S + alpha_ITS + 
    graph16S.cap.metam + 
    scale_x_continuous(limits = c(-1.8, 2.1)) +
    scale_y_continuous(limits = c(-2.5, 2.5)) +
    scale_color_discrete(labels = c("Not fumigated", "Fumigated")) +
    graphITS.cap.metam + 
    scale_x_continuous(limits = c(-1.8, 2)) +
    scale_y_continuous(limits =  c(-2.5, 2.5))  +
    scale_color_discrete(labels = c("Not fumigated", "Fumigated")) +
    plot_layout(guides = "collect") + 
    plot_annotation(tag_levels = "A"))


graph16S.cap.soil + 
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  graphITS.cap.soil + 
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "A")

#-----------------------------------------------------------

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# permanova
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rm(list = ls())

fdata16.genus = readRDS("./RDSdata/fdata16S_genus.rds")
head(sample_data(fdata16.genus))

norm16S = fdata16.genus %>%
  transform_sample_counts(function(x)x/sum(x))

otu.16s = data.frame(t(otu_table(norm16S)))
env.16s = data.frame(sample_data(norm16S))
colnames(env.16s)
head(env.16s)

set.seed(2)
adon.16s = adonis(otu.16s ~ metamStatus + soil_compname + totalCrops + ph, data = env.16s, method = "bray", permutations = 10000)
adon.16s

# dispersion
?phyloseq::distance
bray.16s = phyloseq::distance(physeq = norm16S, method = "bray")

set.seed(2)
disp16s.metam = betadisper(bray.16s, env.16s$metamStatus, add = T, bias.adjust = T)
permutest(disp16s.metam, pairwise = T, permutations = 10000) # 0.0098

set.seed(2)
disp16s.compname = betadisper(bray.16s, env.16s$soil_compname, add = T, bias.adjust = T)
permutest(disp16s.compname, pairwise = T, permutations = 10000) # 0.05909

plot(disp16s.totalCrops)
#@@@@@@@@@@@@@@@@@@
rm(list = ls())
fdataITS.genus = readRDS("./RDSdata/fdataITS_genus.rds")

head(sample_data(fdataITS.genus))

normITS = fdataITS.genus %>%
  transform_sample_counts(function(x)x/sum(x))

otu.its = data.frame(t(otu_table(normITS)))
env.its = data.frame(sample_data(normITS))
str(env.its)

set.seed(3)
adon.its = adonis(otu.its ~ metamStatus + soil_compname + totalCrops + ph, data = env.its, method = "bray", permutations = 10000)
adon.its


bray.its = phyloseq::distance(physeq = normITS, method = "bray")

set.seed(3)
dispITS.metam = betadisper(bray.its, env.its$metamStatus, add = T, bias.adjust = T)
permutest(dispITS.metam, pairwise = T, permutations = 10000) # 0.32


set.seed(3)
dispITS.soil = betadisper(bray.its, env.its$soil_compname, add = T, bias.adjust = T)
permutest(dispITS.soil, pairwise = T, permutations = 10000) # 0.32

set.seed(2)
dispITS.totalRotYrs = betadisper(bray.its, env.its$totalRotYrs, add = T, bias.adjust = T)
permutest(dispITS.totalRotYrs, pairwise = T, permutations = 10000) # 0.0018

set.seed(2)
dispITS.ph = betadisper(bray.its, env.its$ph, add = T, bias.adjust = T)
dispITS.ph
permutest(dispITS.ph, pairwise = T, permutations = 10000) # 0.002
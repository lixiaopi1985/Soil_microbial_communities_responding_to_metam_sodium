Wdir="Your working directory"
setwd(Wdir)

library(phyloseq)
library(rstatix)
library(tidyverse)
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
library(ggpubr)

# rm(list = ls())
# 
# mdata16_bst = readRDS("../Data_prunning/RDS_MCM_data_prunning/mdata16S_cleantaxa_bst2.rds")
# mdataITS_bst = readRDS("../Data_prunning/RDS_MCM_data_prunning/mdataITS_cleantaxa_bst2.rds")
# 
# 
# mdata16_bst.genus = mdata16_bst %>%
#   tax_glom("Genus")
# 
# mdataITS_bst.genus = mdataITS_bst %>%
#   tax_glom("Genus")
# 
# saveRDS(mdata16_bst.genus, "./RDSdata/mdata16S_bst_genus.rds")
# saveRDS(mdataITS_bst.genus, "./RDSdata/mdataITS_bst_genus.rds")

#-----------------------------------------------------------------------
rm(list = ls())

m16.gen = readRDS("./RDSdata/mdata16S_bst_genus.rds")

meta16S = data.frame(sample_data(m16.gen))
colnames(meta16S)[which(colnames(meta16S) == "soil_ph_value")] = "pH"

sample_data(m16.gen) = meta16S

norm16S = m16.gen %>%
  transform_sample_counts(function(x)x/sum(x))

otu.16S = data.frame(t(otu_table(norm16S)))
env.16S = data.frame(sample_data(norm16S))

cap.16S = capscale(otu.16S ~ soil_compname  + 
                      metamStatus + 
                      totalCrops + 
                      pH + 
                      period, 
                      data = env.16S, distance = "bray", add = T)
vif.cca(cap.16S)
cap.16S

cap.16S.select = ordistep(capscale(otu.16S ~ 1, data = env.16S, distance = "bray", add=T), scope = formula(cap.16S), direction = "forward", permutations = 10000)

saveRDS(cap.16S.select, "./RDSdata/cap_16S_select_perm10_000.rds")

cap.16S.select$call


cap.16S.select.padj = cap.16S.select
cap.16S.select.padj$anova$`Pr(>F)` = p.adjust(cap.16S.select$anova$`Pr(>F)`, method = "BH")
cap.16S.select.padj$anova
cap.16S.select.padj

# Df    AIC      F     Pr(>F)    
# + soil_compname  4 2769.0 9.6807 0.00009999 ***
#   + metamStatus    1 2761.1 9.8679 0.00009999 ***
#   + period         3 2754.7 4.1048 0.00009999 ***
#   + totalCrops     1 2751.3 5.3468 0.00009999 ***
#   + pH             1 2749.2 4.0158 0.00009999 ***

set.seed(111)
# cap.16S.final = capscale(otu.16S ~ soil_compname + metamStatus + totalCrops + ph + soil_depth, data = env.16S, distance = "bray", add = T)
cap.16S.final = capscale(otu.16S ~ soil_compname + metamStatus + period + totalCrops + pH, data = env.16S, distance = "bray", add = T)
cap.16S.final
vif.cca(cap.16S.final)

# plot
taxcols.16s = as.data.frame(tax_table(norm16S))
taxcols.16s$ASV = rownames(taxcols.16s)

source("../tools/plotCAP2.R")
(graph16s.cap = plotCAP(cap.16S.final, 
                        env = env.16S, 
                        display = 1, 
                        tax_cols = taxcols.16s, 
                        splitDisplay = F, 
                        sa.xadj = 2.3, 
                        sa.yadj = -0.6, 
                        sa.colorSamples = T, 
                        sa.sample_group = "metamStatus", 
                        sa.showCentroidGroup = "all", 
                        sa.biplot = T,
                        sa.ellipse = F, 
                        sa.addSp = T, 
                        sa.labelsp = T, 
                        sa.pointsize = 6,
                        sa,fontsize = 12,
                        sa.tax_level = "Genus",
                        sa.colorLabs = "Fumigation history"))

(graph16s.cap.T = plotCAP(cap.16S.final, 
                          env = env.16S, 
                          display = 1, 
                          tax_cols = taxcols.16s, 
                          splitDisplay = F, 
                          sa.xadj = 2.3, 
                          sa.yadj = -0.6, 
                          sa.colorSamples = T, 
                          sa.sample_group = "period", 
                          sa.showCentroidGroup = "all", 
                          sa.biplot = T,
                          sa.ellipse = F, 
                          sa.addSp = T, 
                          sa.labelsp = T, 
                          sa.pointsize = 6,
                          sa,fontsize = 12,
                          sa.tax_level = "Genus",
                          sa.colorLabs = "Sampling time"))

(graph16s.cap.soil = plotCAP(cap.16S.final, 
                             env = env.16S, 
                             display = 1, 
                             tax_cols = taxcols.16s, 
                             splitDisplay = F, 
                             sa.xadj = 2.3, 
                             sa.yadj = -0.6, 
                             sa.colorSamples = T, 
                             sa.sample_group = "soil_compname", 
                             sa.showCentroidGroup = "all", 
                             sa.biplot = T,
                             sa.ellipse = F, 
                             sa.addSp = T, 
                             sa.labelsp = T, 
                             sa.pointsize = 6,
                             sa,fontsize = 12,
                             sa.tax_level = "Genus",
                             sa.colorLabs = "Soil series"))


(graph16s.cap.depth = plotCAP(cap.16S.final, 
                              env = env.16S, 
                              display = 1, 
                              tax_cols = taxcols.16s, 
                              splitDisplay = F, 
                              sa.xadj = 2.3, 
                              sa.yadj = -0.6, 
                              sa.colorSamples = T, 
                              sa.sample_group = "soil_depth", 
                              sa.showCentroidGroup = "all", 
                              sa.biplot = T,
                              sa.ellipse = F, 
                              sa.addSp = T, 
                              sa.labelsp = T, 
                              sa.pointsize = 6,
                              sa,fontsize = 12,
                              sa.tax_level = "Genus",
                              sa.colorLabs = "Soil depth"))

saveRDS(graph16s.cap, "./MCM_analysis_output_v3/beta/mcmbeta16S_triplot.rds")
saveRDS(graph16s.cap.soil, "./MCM_analysis_output_v3/beta/soiserisbeta16S_triplot.rds")
saveRDS(graph16s.cap.T, "./MCM_analysis_output_v3/beta/timebeta16S_triplot.rds")
saveRDS(graph16s.cap.depth, "./MCM_analysis_output_v3/beta/soildepthbeta16S_triplot.rds")


set.seed(104)
anova(cap.16S.final, permutations = 10000, step = 1000) # 0.0000999
anova(cap.16S.final, by="terms", permutations = 10000, step = 1000) # all sig

# Permutation test for capscale under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 10000
# 
# Model: capscale(formula = otu.16S ~ soil_compname + metamStatus + period + totalCrops + pH, data = env.16S, distance = "bray", add = T)
# Df SumOfSqs       F     Pr(>F)    
# soil_compname   4   10.170 10.1392 0.00009999 ***
#   metamStatus     1    2.551 10.1712 0.00009999 ***
#   period          3    3.130  4.1604 0.00009999 ***
#   totalCrops      1    1.348  5.3764 0.00009999 ***
#   pH              1    1.007  4.0158 0.00009999 ***
#   Residual      544  136.420                       
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# anova(cap.16S.final, by="axis", permutations = 10000, step = 1000) #3

set.seed(123)

#perm16S = adonis(otu.16S ~ metamStatus + soil_compname +  totalCrops + pH + period, data = env.16S, method = "bray", permutations = 10000)
#saveRDS(perm16S, "./RDSdata/perm16S_metam_soilcomp_totalcrop_ph_period_perm10_000.rds")
perm16S = readRDS("./RDSdata/perm16S_metam_soilcomp_totalcrop_ph_period_perm10_000.rds")
perm16S
# Permutation: free
# Number of permutations: 10000
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
# metamStatus     1     2.014 2.01401  31.897 0.03974 0.00009999 ***
#   soil_compname   4     9.769 2.44224  38.679 0.19276 0.00009999 ***
#   totalCrops      1     1.165 1.16501  18.451 0.02299 0.00009999 ***
#   pH              1     0.829 0.82875  13.125 0.01635 0.00009999 ***
#   period          3     2.553 0.85110  13.479 0.05038 0.00009999 ***
#   Residuals     544    34.349 0.06314         0.67777               
# Total         554    50.679                 1.00000       


# dispersion
?phyloseq::distance
bray.16s = phyloseq::distance(physeq = norm16S, method = "bray")

set.seed(2)
disp16s.metam = betadisper(bray.16s, env.16S$metamStatus, add = T, bias.adjust = T)
permutest(disp16s.metam, pairwise = T, permutations = 10000) # <0.0001

set.seed(2)
disp16s.compname = betadisper(bray.16s, env.16S$soil_compname, add = T, bias.adjust = T)
permutest(disp16s.compname, pairwise = T, permutations = 10000) # <0.0001

set.seed(2)
disp16s.period = betadisper(bray.16s, env.16S$period, add = T, bias.adjust = T)
permutest(disp16s.period, pairwise = T, permutations = 10000) # <0.0001

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# fungi
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

rm(list = ls())
mITS.gen = readRDS("./RDSdata/mdataITS_bst_genus.rds")

metaITS = data.frame(sample_data(mITS.gen))
colnames(metaITS)[which(colnames(metaITS) == "soil_ph_value")] = "pH"
sample_data(mITS.gen) = metaITS

normITS = mITS.gen %>%
  transform_sample_counts(function(x)x/sum(x))


otu.ITS = data.frame(t(otu_table(normITS)))
env.ITS = data.frame(sample_data(normITS))

cap.ITS = capscale(otu.ITS ~ soil_compname + 
                     metamStatus + 
                     totalCrops + 
                     pH + 
                     period, data = env.ITS, distance = "bray", add = T)
vif.cca(cap.ITS)

cap.ITS.select = ordistep(capscale(otu.ITS ~ 1, data = env.ITS, distance = "bray", add=T), scope = formula(cap.ITS), direction = "forward", permutations = 10000)
saveRDS(cap.ITS.select, "./RDSdata/cap_ITS_select_perm10_000.rds")

cap.ITS.select$call
cap.ITS.select.padj = cap.ITS.select
cap.ITS.select.padj$anova$`Pr(>F)` = p.adjust(cap.ITS.select$anova$`Pr(>F)`, method = "BH")
cap.ITS.select.padj$anova
cap.ITS.select.padj
# + soil_compname  4 3256.8  6.4422 0.00009999 ***
#   + metamStatus    1 3247.6 11.2124 0.00009999 ***
#   + period         3 3243.7  3.2659 0.00009999 ***
#   + totalCrops     1 3242.5  3.1665 0.00009999 ***
#   + pH             1 3241.2  3.2575 0.00009999 ***

set.seed(1)
cap.ITS.final = capscale(otu.ITS ~ soil_compname  + metamStatus + totalCrops + pH + period, data = env.ITS, distance = "bray", add = T)
cap.ITS.final
vif.cca(cap.ITS.final)



# plot
taxcols.its = as.data.frame(tax_table(normITS))
taxcols.its$ASV = rownames(taxcols.its)

source("../tools/plotCAP2.R")
(graphITS.cap = plotCAP(cap.ITS.final, 
                        env = env.ITS, 
                        display = 1, 
                        tax_cols = taxcols.its, 
                        splitDisplay = F, 
                        sa.xadj = 2, 
                        sa.yadj = -0.4, 
                        sa.colorSamples = T, 
                        sa.sample_group = "metamStatus", 
                        sa.showCentroidGroup = "all", 
                        sa.biplot = T,
                        sa.ellipse = F, 
                        sa.addSp = T, 
                        sa.labelsp = T, 
                        sa.pointsize = 6,
                        sa,fontsize = 12,
                        sa.tax_level = "Genus",
                        sa.colorLabs = "Fumigation history"))


(graphITS.cap.T = plotCAP(cap.ITS.final, 
                          env = env.ITS, 
                          display = 1, 
                          tax_cols = taxcols.its, 
                          splitDisplay = F, 
                          sa.xadj = 2, 
                          sa.yadj = -0.4, 
                          sa.colorSamples = T, 
                          sa.sample_group = "period", 
                          sa.showCentroidGroup = "all", 
                          sa.biplot = T,
                          sa.ellipse = F, 
                          sa.addSp = T, 
                          sa.labelsp = T, 
                          sa.pointsize = 6,
                          sa,fontsize = 12,
                          sa.tax_level = "Genus",
                          sa.colorLabs = "Sampling time"))

(graphITS.cap.soil = plotCAP(cap.ITS.final, 
                             env = env.ITS, 
                             display = 1, 
                             tax_cols = taxcols.its, 
                             splitDisplay = F, 
                             sa.xadj = 2, 
                             sa.yadj = -0.4, 
                             sa.colorSamples = T, 
                             sa.sample_group = "soil_compname", 
                             sa.showCentroidGroup = "all", 
                             sa.biplot = T,
                             sa.ellipse = F, 
                             sa.addSp = T, 
                             sa.labelsp = T, 
                             sa.pointsize = 6,
                             sa,fontsize = 12,
                             sa.tax_level = "Genus",
                             sa.colorLabs = "Soil series"))

saveRDS(graphITS.cap, "./MCM_analysis_output_v3/beta/mcmbetaITS_triplot.rds")
saveRDS(graphITS.cap.T, "./MCM_analysis_output_v3/beta/mcm_timeITS_triplot.rds")
saveRDS(graphITS.cap.soil, "./MCM_analysis_output_v3/beta/mcm_soilITS_triplot.rds")


set.seed(104)
anova(cap.ITS.final, permutations = 10000, step = 1000) # 0.001
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 10000
# 
# Model: capscale(formula = otu.ITS ~ soil_compname + metamStatus + totalCrops + ph + period, data = env.ITS, distance = "bray", add = T)
# Df SumOfSqs      F     Pr(>F)    
# Model     10   28.707 5.4517 0.00009999 ***
#   Residual 556  292.774                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

anova(cap.ITS.final, by="terms", permutations = 10000, step = 1000) # all 0.001

# Permutation test for capscale under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 10000
# 
# Model: capscale(formula = otu.ITS ~ soil_compname + metamStatus + totalCrops + ph + period, data = env.ITS, distance = "bray", add = T)
# Df SumOfSqs       F     Pr(>F)    
# soil_compname   4   14.094  6.6916 0.00009999 ***
#   metamStatus     1    6.023 11.4385 0.00009999 ***
#   totalCrops      1    1.675  3.1810 0.00009999 ***
#   ph              1    2.104  3.9963 0.00009999 ***
#   period          3    4.810  3.0451 0.00009999 ***
#   Residual      556  292.774                       
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
anova(cap.ITS.final, by="axis", permutations = 10000, step = 1000) #3
# Permutation test for capscale under reduced model
# Forward tests for axes
# Permutation: free
# Number of permutations: 10000
# 
# Model: capscale(formula = otu.ITS ~ soil_compname + metamStatus + totalCrops + ph + period, data = env.ITS, distance = "bray", add = T)
# Df SumOfSqs       F     Pr(>F)    
# CAP1       1    9.864 18.7333 0.00009999 ***
#   CAP2       1    5.223  9.9187 0.00009999 ***
#   CAP3       1    3.749  7.1199 0.00009999 ***
#   CAP4       1    2.909  5.5253 0.00009999 ***
#   CAP5       1    1.823  3.4626 0.00009999 ***
#   CAP6       1    1.576  2.9936 0.00009999 ***
#   CAP7       1    1.216  2.3091 0.00009999 ***
#   CAP8       1    0.997  1.8938 0.00009999 ***
#   CAP9       1    0.683  1.2966     0.0176 *  
#   CAP10      1    0.666  1.2644     0.0048 ** 
#   Residual 556  292.774                       
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 

set.seed(123)
#permITS = adonis(otu.ITS ~ metamStatus + soil_compname  +  totalCrops + pH + period, data = env.ITS, method = "bray", permutations = 10000)
#saveRDS(permITS, "./RDSdata/permITS_metam_soilcomp_totalcrop_ph_period_perm10_000.rds")
permITS = readRDS("./RDSdata/permITS_metam_soilcomp_totalcrop_ph_period_perm10000.rds")
permITS

# Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
# metamStatus     1     5.515  5.5155  34.116 0.04798 0.00009999 ***
#   soil_compname   4    12.778  3.1944  19.759 0.11116 0.00009999 ***
#   totalCrops      1     1.310  1.3101   8.104 0.01140 0.00009999 ***
#   pH              1     1.739  1.7394  10.759 0.01513 0.00009999 ***
#   period          3     3.716  1.2386   7.661 0.03233 0.00009999 ***
#   Residuals     556    89.888  0.1617         0.78200               
# Total         566   114.946                 1.00000               
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# dispersion
bray.its = phyloseq::distance(physeq = normITS, method = "bray")

set.seed(2)
dispITS.metam = betadisper(bray.its, env.ITS$metamStatus, add = T, bias.adjust = T)
permutest(dispITS.metam, pairwise = T, permutations = 10000) # < 0.0001

set.seed(2)
dispITS.compname = betadisper(bray.its, env.ITS$soil_compname, add = T, bias.adjust = T)
permutest(dispITS.compname, pairwise = T, permutations = 10000) # 

set.seed(2)
dispITS.time = betadisper(bray.its, env.ITS$period, add = T, bias.adjust = T)
permutest(dispITS.time, pairwise = T, permutations = 10000) # 

#################################################################################
# plot all
library(patchwork)
rm(list = ls())

# from metam sodium

alpha16S.obs = readRDS("./MCM_analysis_output_v3/diversity/plot_obs_mcm_16S.rds")
alphaITS.obs = readRDS("./MCM_analysis_output_v3/diversity/plot_obs_mcm_ITS.rds")

m16.beta.graph = readRDS("./MCM_analysis_output_v3/beta/mcmbeta16S_triplot.rds")
mITS.beta.graph = readRDS("./MCM_analysis_output_v3/beta//mcmbetaITS_triplot.rds")

(alpha16S.obs + alphaITS.obs + m16.beta.graph + 
    scale_x_continuous(limits = c(-2.8, 2.5)) + 
    scale_y_continuous(limits = c(-2.2, 3)) +
    # scale_color_discrete(labels=c("Not fumigated", "Fumigated")) +
    # labs(title = "A") + 
    #theme(plot.title = element_text(hjust = -0.05)) + 
    mITS.beta.graph + 
    scale_x_continuous(limits = c(-2.8, 2.5))+ 
    scale_y_continuous(limits = c(-2.2, 3.2)) +
    # scale_color_discrete(labels=c("Not fumigated", "Fumigated")) +
    #labs(title = "B") + 
    #theme(plot.title = element_text(hjust = -0.05)) + 
    plot_layout(guides = "collect", heights = c(1,2)) + plot_annotation(tag_levels = "A"))


bysoil16S = readRDS("./MCM_analysis_output_v3/beta/soiserisbeta16S_triplot.rds")
bysoil16S
bysoilITS = readRDS("./MCM_analysis_output_v3/beta/mcm_soilITS_triplot.rds")
bysoilITS

byTime.16S = readRDS("./MCM_analysis_output_v3/beta/timebeta16S_triplot.rds")
byTime.16S
byTime.ITS = readRDS("./MCM_analysis_output_v3/beta/mcm_timeITS_triplot.rds")
byTime.ITS



(bysoil16S + 
    scale_x_continuous(limits = c(-3, 2.5)) + 
    scale_y_continuous(limits = c(-3, 3.2)) +
    labs(title = "A") + 
    theme(plot.title = element_text(hjust = -0.05)) + 
    bysoilITS + 
    scale_x_continuous(limits = c(-3, 2.5))+ 
    scale_y_continuous(limits = c(-3, 3.2)) +
    labs(title = "B") + 
    theme(plot.title = element_text(hjust = -0.05)) + 
    byTime.16S +
    scale_color_manual(values = c("brown2", "dodgerblue", "olivedrab4", "orchid4")) +
    scale_x_continuous(limits = c(-3, 2.5))+ 
    scale_y_continuous(limits = c(-3, 3.2)) +
    labs(title = "C") + 
    theme(plot.title = element_text(hjust = -0.05)) + 
    byTime.ITS +
    scale_color_manual(values = c("brown2", "dodgerblue", "olivedrab4", "orchid4")) +
    scale_x_continuous(limits = c(-3, 2.5))+ 
    scale_y_continuous(limits = c(-3, 3.2)) +
    labs(title = "D") + 
    theme(plot.title = element_text(hjust = -0.05)) + 
    plot_layout(guides = "collect"))




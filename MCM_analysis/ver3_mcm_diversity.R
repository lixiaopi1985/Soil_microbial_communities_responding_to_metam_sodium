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


rm(list = ls())


mdata16_bst = readRDS("../Data_prunning/RDS_MCM_data_prunning/mdata16S_cleantaxa_bst2.rds")
mdataITS_bst = readRDS("../Data_prunning/RDS_MCM_data_prunning/mdataITS_cleantaxa_bst2.rds")

mdata16_bst

# otu_table()   OTU Table:         [ 44549 taxa and 555 samples ]
# sample_data() Sample Data:       [ 555 samples by 66 sample variables ]
# tax_table()   Taxonomy Table:    [ 44549 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 44549 tips and 43913 internal nodes ]

mdataITS_bst

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6063 taxa and 567 samples ]
# sample_data() Sample Data:       [ 567 samples by 66 sample variables ]
# tax_table()   Taxonomy Table:    [ 6063 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 6063 tips and 5885 internal nodes ]

taxadf16S = data.frame(tax_table(mdata16_bst))
taxadfITS = data.frame(tax_table(mdataITS_bst))


vegan::rarecurve(t(phyloseq::otu_table(mdata16_bst)), step = 100, label = F, ylab = "ASV counts", xlab = "Sequencing depth (Reads)")

vegan::rarecurve(t(phyloseq::otu_table(mdataITS_bst)), step = 100, label = F, ylab = "ASV counts", xlab = "Sequencing depth (Reads)")

meta16S = data.frame(sample_data(mdata16_bst))
meta16S$sampleid = gsub("-", "_", meta16S$sampleid)

metaITS = data.frame(sample_data(mdataITS_bst))
metaITS$sampleid = gsub("-", "_", metaITS$sampleid)

min(sample_sums(mdata16_bst)) #7019
min(sample_sums(mdataITS_bst)) #10476

rare16S = rarefy_even_depth(mdata16_bst, sample.size = min(sample_sums(mdata16_bst)), rngseed = 123, replace = F)
dv16S = estimate_richness(rare16S, measures = c("Observed", "Shannon"))
head(dv16S)
rownames(dv16S) = gsub("\\.", "_", rownames(dv16S))

rareITS = rarefy_even_depth(mdataITS_bst, sample.size = min(sample_sums(mdataITS_bst)), rngseed = 123, replace = F)
dvITS = estimate_richness(rareITS, measures = c("Observed", "Shannon"))
head(dvITS)
rownames(dvITS) = gsub("\\.", "_", rownames(dvITS))
                       
af16S = merge(dv16S, meta16S, by.x = "row.names", by.y = "sampleid")
afITS = merge(dvITS, metaITS, by.x = "row.names", by.y = "sampleid")


af16S.1 = af16S %>%
  gather("alpha_index", "index_measure", c("Observed", "Shannon"))

afITS.1 = afITS %>%
  gather("alpha_index", "index_measure", c("Observed", "Shannon"))

#saveRDS(af16S.1, "./RDSdata/aM16S.rds")
#saveRDS(afITS.1, "./RDSdata/aMITS.rds")

#############################################################################
# statistical analysis
##############################################################################
rm(list=ls())

af16S.1 = readRDS("./RDSdata/aM16S.rds")
afITS.1 = readRDS("./RDSdata/aMITS.rds")

af16S.1$metamStatus


alpha.16s.obs = af16S.1 %>%
  filter(alpha_index == "Observed")

alpha.its.obs = afITS.1 %>%
  filter(alpha_index == "Observed")

# testing on metam history and time and soil depth

mod16S.obs =  aov(index_measure ~ metamStatus*period, data = alpha.16s.obs)
out_mod16S.obs = summary(mod16S.obs)
out_mod16S.obs
plot(resid(mod16S.obs))

modITS.obs = aov(index_measure ~ metamStatus*period, data = alpha.its.obs)
out_modITS.obs = summary(modITS.obs)
out_modITS.obs
plot(resid(modITS.obs))

# boxplots

timePatte_16S = c("#001933", "#004C99", "#0080FF", "#99CCFF")

(box16S.obs = ggboxplot(alpha.16s.obs, x = "metamStatus", y="index_measure", fill="period", alpha=0.7, size=1, ylab = "Observed", xlab = "Fumigation history"))

(stat16S.obs = alpha.16s.obs %>%
    dplyr::group_by(metamStatus) %>%
    rstatix::t_test(index_measure ~ period, ref.group = "0wk", p.adjust.method = "hochberg", detailed = T) %>%
    add_xy_position(x = "metamStatus", dodge = 0.8))

stat16S.obs


(b16S.obs = box16S.obs +
    stat_pvalue_manual(
      stat16S.obs, 
      label = "p.adj.signif", tip.length = 0.01, bracket.nudge.y = 0.2,
      step.increase = 0
    ) +
    # scale_y_continuous(limits = c(minObs, maxObs+500)) +
    scale_fill_manual(values = timePatte_16S, labels =c("pre", "1 week", "3 weeks", "6 weeks")) +
    labs(fill = "Time") +
    theme(legend.position = "right",
          text = element_text(size=12),
          axis.text.x = element_text(size=12, color="black"),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_text(size=12, margin = margin(t=5)),
          axis.title.y = element_text(size=12, margin = margin(r=5))))

#---------------------------------------------------------------------------
# ITS
#---------------------------------------------------------------------------

(boxITS.obs = ggboxplot(alpha.its.obs, x = "metamStatus", y="index_measure", fill="period", alpha=0.7, size=1, ylab = "", xlab = "Fumigation history"))


(statITS.obs = alpha.its.obs %>%
    dplyr::group_by(metamStatus) %>%
    rstatix::t_test(index_measure ~ period, ref.group = "0wk", p.adjust.method = "hochberg", detailed = T) %>%
    add_xy_position(x = "metamStatus", dodge = 0.8))

statITS.obs %>%
  dplyr::select(metamStatus, estimate, p.adj)

(bITS.obs = boxITS.obs +
    stat_pvalue_manual(
      statITS.obs, 
      label = "p.adj.signif", 
      tip.length = 0.01, 
      bracket.nudge.y = 0.2,
      step.increase = 0.05
    ) +
    # scale_y_continuous(limits = c(minObs, maxObs)) +
    scale_fill_manual(values = timePatte_16S, labels =c("pre", "1 week", "3 weeks", "6 weeks")) +
    labs(fill = "Time") +
    theme(legend.position = "right",
          text = element_text(size=12),
          axis.text.x = element_text(size=12, color="black"),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_text(size=12, margin = margin(t=5)),
          axis.title.y = element_text(size=12, margin = margin(r=5))))

#saveRDS(b16S.obs, "./MCM_analysis_output_v3/diversity/plot_obs_mcm_16S.rds")
#saveRDS(bITS.obs, "./MCM_analysis_output_v3/diversity/plot_obs_mcm_ITS.rds")


#----------------------------------------------------------
#----------------------------------------------------------
# Shannon 
#---------------------------------------------------------
#----------------------------------------------------------
rm(list = ls())
af16S.1 = readRDS("./RDSdata/aM16S.rds")
afITS.1 = readRDS("./RDSdata/aMITS.rds")

af16S.1$metamStatus

alpha.16s.shan = af16S.1 %>%
  filter(alpha_index == "Shannon")
head(alpha.16s.shan)

alpha.its.shan = afITS.1 %>%
  filter(alpha_index == "Shannon")

# testing on metam history and time and soil depth

mod16S.shan =  aov(index_measure ~ metamStatus*period, data = alpha.16s.shan)
out_mod16S.shan = summary(mod16S.shan)
out_mod16S.shan

plot(resid(mod16S.shan))


modITS.shan = aov(index_measure ~ metamStatus*period, data = alpha.its.shan)
out_modITS.shan = summary(modITS.shan)
out_modITS.shan

timePatte_16S = c("#001933", "#004C99", "#0080FF", "#99CCFF")

#--------------
(box16S.shan = ggboxplot(alpha.16s.shan, x = "metamStatus", y="index_measure", fill="period", alpha=0.7, size=1, ylab = "Shannon", xlab = "Fumigation history"))

(stat16S.shan = alpha.16s.shan %>%
    dplyr::group_by(metamStatus) %>%
    rstatix::t_test(index_measure ~ period, ref.group = "0wk", p.adjust.method = "hochberg", detailed = T) %>%
    add_xy_position(x = "metamStatus", dodge = 0.8))

stat16S.shan$p.adj


(b16S.shan = box16S.shan +
    stat_pvalue_manual(
      stat16S.shan, 
      label = "p.adj.signif", tip.length = 0.01, bracket.nudge.y = 0.2,
      step.increase = 0
    ) +
    # scale_y_continuous(limits = c(minObs, maxObs+500)) +
    scale_fill_manual(values = timePatte_16S, labels =c("pre", "1 week", "3 weeks", "6 weeks")) +
    labs(fill = "Time") +
    theme(legend.position = "right",
          text = element_text(size=12),
          axis.text.x = element_text(size=12, color="black"),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_text(size=12, margin = margin(t=5)),
          axis.title.y = element_text(size=12, margin = margin(r=5))))

#-----------------------------------------------------------------------------
(boxITS.shan = ggboxplot(alpha.its.shan, x = "metamStatus", y="index_measure", fill="period", alpha=0.7, size=1, ylab = "", xlab = "Fumigation history"))


(statITS.shan = alpha.its.shan %>%
    dplyr::group_by(metamStatus) %>%
    rstatix::t_test(index_measure ~ period, ref.group = "0wk", p.adjust.method = "hochberg", detailed = T) %>%
    add_xy_position(x = "metamStatus", dodge = 0.8))

statITS.shan %>%
  dplyr::select(metamStatus, estimate, p.adj)

(bITS.shan = boxITS.shan +
    stat_pvalue_manual(
      statITS.shan, 
      label = "p.adj.signif", 
      tip.length = 0.01, 
      bracket.nudge.y = 0.2,
      step.increase = 0.05
    ) +
    # scale_y_continuous(limits = c(minObs, maxObs)) +
    scale_fill_manual(values = timePatte_16S, labels =c("pre", "1 week", "3 weeks", "6 weeks")) +
    labs(fill = "Time") +
    theme(legend.position = "right",
          text = element_text(size=12),
          axis.text.x = element_text(size=12, color="black"),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_text(size=12, margin = margin(t=5)),
          axis.title.y = element_text(size=12, margin = margin(r=5))))

saveRDS(b16S.shan, "./MCM_analysis_output_v3/diversity/plot_shan_mcm_16S.rds")
saveRDS(bITS.shan, "./MCM_analysis_output_v3/diversity/plot_shan_mcm_ITS.rds")

#--------------------------------------------------------------------------------

b16S.shan = readRDS("./MCM_analysis_output_v3/diversity/plot_shan_mcm_16S.rds")
bITS.shan = readRDS("./MCM_analysis_output_v3/diversity/plot_shan_mcm_ITS.rds")


b16S.shan + bITS.shan + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")

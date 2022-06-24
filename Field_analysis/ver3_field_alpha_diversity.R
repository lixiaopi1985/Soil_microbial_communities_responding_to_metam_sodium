WDir="E:/OSU_project_2017_2019/mcm_paper/reviewer_comments/MCM_phytobiome/ver3/Field_analysis"
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
library(rstatix)


rm(list = ls())


fdata16S = readRDS("../Data_prunning/RDS_Field_data_prunning/fdata16S_cleantaxa_bst2.rds")
fdataITS = readRDS("../Data_prunning/RDS_Field_data_prunning/fdataITS_cleantaxa_bst2.rds")

meta16S = data.frame(sample_data(fdata16S))
metaITS = data.frame(sample_data(fdataITS))

fdata16S
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 21899 taxa and 48 samples ]
# sample_data() Sample Data:       [ 48 samples by 61 sample variables ]
# tax_table()   Taxonomy Table:    [ 21899 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 21899 tips and 21641 internal nodes ]
fdataITS
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 3028 taxa and 48 samples ]
# sample_data() Sample Data:       [ 48 samples by 61 sample variables ]
# tax_table()   Taxonomy Table:    [ 3028 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 3028 tips and 2985 internal nodes ]

taxdf.16S = data.frame(tax_table(fdata16S))
length(unique(taxdf.16S$Genus)) #1241


taxdf.ITS = data.frame(tax_table(fdataITS))
length(unique(taxdf.ITS$Genus)) #481


min(sample_sums(fdata16S)) #39275
min(sample_sums(fdataITS)) #29042

# figure S1
vegan::rarecurve(t(phyloseq::otu_table(fdata16S)), step = 100, label = F, ylab = "ASV counts", xlab = "Sequencing depth (Reads)")
vegan::rarecurve(t(phyloseq::otu_table(fdataITS)), step = 100, label = F, ylab = "ASV counts", xlab = "Sequencing depth (Reads)")


rare16S = rarefy_even_depth(fdata16S, sample.size = min(sample_sums(fdata16S)), rngseed = 123, replace = F)
dv16S = estimate_richness(rare16S, measures = c("Observed", "Shannon"))

rareITS = rarefy_even_depth(fdataITS, sample.size = min(sample_sums(fdataITS)), rngseed = 123, replace = F)
dvITS = estimate_richness(rareITS, measures = c("Observed", "Shannon"))

af16S = merge(dv16S, meta16S, by.x = "row.names", by.y = "SampleID")
afITS = merge(dvITS, metaITS, by.x = "row.names", by.y = "SampleID")

af16S.1 = af16S %>%
  gather("alpha_index", "index_measure", c("Observed", "Shannon"))

afITS.1 = afITS %>%
  gather("alpha_index", "index_measure", c("Observed", "Shannon"))


af16S.1$metamStatus

# reorder
af16S.1 = af16S.1[c(which(af16S.1$metamStatus == "Fumigated"), which(af16S.1$metamStatus == "Not fumigated")),]
afITS.1 = afITS.1[c(which(afITS.1$metamStatus == "Fumigated"), which(afITS.1$metamStatus == "Not fumigated")),]

af16S.1$field_alias = factor(af16S.1$field_alias, levels = unique(af16S.1$field_alias))
afITS.1$field_alias = factor(afITS.1$field_alias, levels = unique(afITS.1$field_alias))

af16S.1$metamRates
afITS.1$metamRates


af16S.1$metamAvg = ifelse(af16S.1$metamRates != 0, af16S.1$metamRates / af16S.1$metamApps, 0)
af16S.1$metamAvg
afITS.1$metamAvg = ifelse(afITS.1$metamRates != 0, afITS.1$metamRates / afITS.1$metamApp, 0)
afITS.1$metamAvg

#saveRDS(af16S.1, "./RDSdata/af16S.rds")
#saveRDS(afITS.1, "./RDSdata/afITS.rds")

#-----------------------------------------------------------------------------
# statistical analysis
rm(list=ls())

af16S.1 = readRDS("./RDSdata/af16S.rds")
afITS.1 = readRDS("./RDSdata/afITS.rds")

testON = "Observed"
testON = "Shannon"

data.obs.16s = af16S.1 %>%
  filter(alpha_index == testON) 
# 13 groups

head(data.obs.16s)

(stat16S = data.obs.16s %>%
  group_by(field_alias) %>%
  summarise(avg = mean(index_measure), n=n(), var = var(index_measure)))

kwpower(stat16S$n, stat16S$avg, "normal")



effsize16S = kruskal_effsize(data = data.obs.16s, index_measure ~ field_alias, ci=T)
effsize16S

#  Observed
# .y.               n effsize conf.low conf.high method  magnitude
# * <chr>         <int>   <dbl>    <dbl>     <dbl> <chr>   <ord>    
#   1 index_measure    48   0.625      0.5      0.84 eta2[H] large   

# Shannon
# .y.               n effsize conf.low conf.high method  magnitude
# * <chr>         <int>   <dbl>    <dbl>     <dbl> <chr>   <ord>    
#   1 index_measure    48   0.649     0.51      0.88 eta2[H] large    

data.obs.its = afITS.1 %>%
  filter(alpha_index == testON)

effsizeITS = kruskal_effsize(data = data.obs.its, index_measure ~ field_alias, ci=T)
effsizeITS

# observed
# .y.               n effsize conf.low conf.high method  magnitude
# * <chr>         <int>   <dbl>    <dbl>     <dbl> <chr>   <ord>    
#   1 index_measure    48   0.155     0.07      0.69 eta2[H] large  

# Shannon
# .y.               n effsize conf.low conf.high method  magnitude
# * <chr>         <int>   <dbl>    <dbl>     <dbl> <chr>   <ord>    
#   1 index_measure    48  0.0374     0.02       0.6 eta2[H] small 


head(data.obs.16s)

kw16S.obs = kruskal.test(index_measure ~ field_alias, data = data.obs.16s)
kw16S.obs
kw16S.obs$p.value # obs 0.0007037 shannon: 0.0005165
dun16S.obs = FSA::dunnTest(index_measure~field_alias,  data = data.obs.16s, method = "bh")
dun16S_res = as.data.frame(dun16S.obs$res)
str(dun16S_res)

# reorder the factor levels in the data set so that the group with the largest mean or median is first
dun16S_res_order = dun16S_res[order(dun16S_res$Z),]
dun16S_res_order
dun16S_res_order[dun16S_res_order$P.adj<0.05,]
letters16S = cldList(P.adj ~ Comparison, data = dun16S_res_order, threshold = 0.05, remove.zero = F, print.comp = T)
letters16S
letters16S.merge = merge(letters16S, aggregate(index_measure ~ field_alias, data=data.obs.16s, max), by.x = "Group", by.y = "field_alias")
head(letters16S.merge)

kwITS.obs = kruskal.test(index_measure ~ field_alias, data = data.obs.its) # 0.1344
kwITS.obs
kwITS.obs$p.value # OTU 0.134, Shannon 0.347
dunITS.obs = FSA::dunnTest(index_measure~field_alias,  data = data.obs.its, method = "bh")
dunITS.obs
lettersITS = cldList(comparison = dunITS.obs$res$Comparison, p.value = dunITS.obs$res$P.adj, threshold = 0.05, remove.zero = F, decreasing =T)

# plot alpha diversity
plabel16S = paste0("italic(P)==", round(kw16S.obs$p.value, 4))
plabelITS = paste0("italic(P)==", round(kwITS.obs$p.value, 4))

if(testON == "Observed"){
  gap = 55
  ylabs = "Observed richness"
} else if(testON == "Shannon") {
  gap = 0.2
  ylabs = testON
}

(obs16S = data.obs.16s %>%
    ggplot(aes(x = field_alias, y = index_measure, fill = metamStatus)) +
    geom_boxplot( alpha = 0.7, size=1) +
    annotate("text", x = letters16S.merge$Group, y = letters16S.merge$index_measure+gap, label=letters16S.merge$Letter) +
    annotate("text", x = 2, y=max(letters16S.merge$index_measure)+gap, label = plabel16S, parse=T) +
    scale_y_continuous(limits = c(min(data.obs.its$index_measure), max(letters16S.merge$index_measure)+gap)) +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2") +
    labs(fill = "Fumigation history", x = "Fields", y = ylabs) +
    theme(
      text =  element_text(size=12, color="black"),
      legend.text = element_text(size=12),
      axis.title.x = element_text(margin = margin(t=5)),
      axis.title.y = element_text(margin = margin(r=5)),
      axis.text.x = element_text(color="black", size = 12),
      axis.text.y = element_text(color="black", size = 12)))

(obsITS = data.obs.its %>%
    ggplot(aes(x = field_alias, y = index_measure, fill = metamStatus)) +
    geom_boxplot(alpha = 0.7, size=1) +
    annotate("text", x = 2, y= max(letters16S.merge$index_measure)+gap, label = plabelITS, parse=T) +
    scale_y_continuous(limits = c(min(data.obs.its$index_measure), max(letters16S.merge$index_measure)+gap)) +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2") +
    labs(fill = "Fumigation history", x = "Fields", y = "") +
    theme(
      text =  element_text(size=12, color="black"),
      legend.text = element_text(size=12),
      axis.title.x = element_text(margin = margin(t=5)),
      axis.title.y = element_text(margin = margin(r=5)),
      axis.text.x = element_text(color="black", size = 12),
      axis.text.y = element_text(color="black", size = 12)))


# saveRDS(obs16S, "./Field_analysis_output_v3/diversity/plot_obs16S_field.rds")
# saveRDS(obsITS, "./Field_analysis_output_v3/diversity/plot_obsITS_field.rds")

saveRDS(obs16S, "./Field_analysis_output_v3/diversity/plot_shan16S_field.rds")
saveRDS(obsITS, "./Field_analysis_output_v3/diversity/plot_shanITS_field.rds")


data.obs.16s %>%
  group_by(field_alias) %>%
  summarise(mean.index = mean(index_measure, na.rm=T),
            sd.index = sd(index_measure, na.rm=T),
            n.index = n()) %>%
  mutate(se.index = sd.index / sqrt(n.index),
         lower.ci = mean.index - qnorm(0.975)*se.index,
         upper.ci = mean.index + qnorm(0.975)*se.index) %>%
  select(field_alias, mean.index, lower.ci, upper.ci)




#----------------------------------------------------------------
# other factors
#--------------------------------------------------------------
rm(list = ls())

af16S.1 = readRDS("./RDSdata/af16S.rds")
afITS.1 = readRDS("./RDSdata/afITS.rds")

testON = "Shannon"
data.obs.16s = af16S.1 %>%
  filter(alpha_index == testON)

data.obs.its = afITS.1 %>%
  filter(alpha_index == testON)

################
# metam sodium
################

kw16S.obs.metam = kruskal.test(index_measure ~ metamStatus, data = data.obs.16s)
kw16S.obs.metam$p.value # 0.5977 # 0.0859

kwITS.obs.metam = kruskal.test(index_measure ~ metamStatus, data = data.obs.its) # 0.1344
kwITS.obs.metam$p.value # 0.185 #0.341


# plot alpha diversity
plabel16S = paste0("italic(P)==", round(kw16S.obs.metam$p.value, 4))
plabelITS = paste0("italic(P)==", round(kwITS.obs.metam$p.value, 4))

if(testON == "Observed"){
  gap = 50
  ylabs = "Observed richness"
} else if(testON == "Shannon") {
  gap = 0.1
  ylabs = testON
}

source("../tools/plotReg.R")
minObs = min(c(data.obs.its$index_measure, data.obs.16s$index_measure))
maxObs = max(c(data.obs.its$index_measure, data.obs.16s$index_measure))

(obs16S.metam =  plotBox(data.obs.16s, x = "metamStatus", y="index_measure", xLabs = "Fumigation history", yLabs = ylabs, minLimit = minObs, maxLimit = maxObs, gap = 0.1, formulaLabs = plabel16S))
(obsITS.metam =  plotBox(data.obs.its, x = "metamStatus", y="index_measure", xLabs = "Fumigation history",yLabs = ylabs,  minLimit = minObs, maxLimit = maxObs, gap = 0.1, formulaLabs = plabelITS))


# by soil series

kw16S.obs.soil = kruskal.test(index_measure ~ soil_compname, data = data.obs.16s)
kw16S.obs.soil$p.value # obs 0.2589, shannon 0.1725

plabel16S.soil = paste0("italic(P)==", round(kw16S.obs.soil$p.value, 4))

(obs16S.soil =  plotBox(data.obs.16s, x = "soil_compname", y="index_measure", xLabs = "Soil series", yLabs = ylabs, minLimit = minObs, maxLimit = maxObs, gap = 0.05, formulaLabs = plabel16S.soil, label.pos.x = 1))


kwITS.obs.soil = kruskal.test(index_measure ~ soil_compname, data = data.obs.its) # 0.01186
kwITS.obs.soil$p.value  # obs: 0.01187, shannon: 0.4568
plabelITS.soil = paste0("italic(P)==", round(kwITS.obs.soil$p.value, 4))

if(testON == "Shannon") {
  (obsITS.soil =  plotBox(data.obs.its, x = "soil_compname", y="index_measure", xLabs = "Soil series", yLabs = ylabs, minLimit = minObs, maxLimit = maxObs, gap = 0.1, formulaLabs = plabelITS.soil, extraAnno = F, Anno2 = lettersITS.soil.merge, label.pos.x = 1))
} else {
  dunITS.soil = FSA::dunnTest(index_measure~soil_compname,  data = data.obs.its, method = "bh")
  dunITS_res.soil = as.data.frame(dunITS.soil$res)
  str(dunITS_res.soil)
  # reorder the factor levels in the data set so that the group with the largest mean or median is first
  dunITS_res_order.soil = dunITS_res.soil[order(dunITS_res.soil$Z),]
  dunITS_res_order.soil
  lettersITS.soil = cldList(P.adj ~ Comparison, data = dunITS_res_order.soil, threshold = 0.05, remove.zero = F, print.comp = T)
  lettersITS.soil
  
  lettersITS.soil.merge = merge(lettersITS.soil, aggregate(index_measure ~ soil_compname, data=data.obs.its, max), by.x = "Group", by.y = "soil_compname")
  head(lettersITS.soil.merge)
  
  (obsITS.soil =  plotBox(data.obs.its, x = "soil_compname", y="index_measure", xLabs = "Soil series", yLabs = ylabs, minLimit = minObs, maxLimit = maxObs, gap = 10, formulaLabs = plabelITS.soil, extraAnno = T, Anno2 = lettersITS.soil.merge, label.pos.x = 1))
}

obsITS.soil

if(testON == "Observed"){
  lingap = 500
  extragap = 150
}else{
  lingap = 1
  extragap = 0.5
}

# by pH
(div.ph.16 = plotReg(data.obs.16s, 
                     x="ph", y="index_measure", minLimit = minObs, maxLimit = maxObs, gap = lingap, xLabs = "Soil pH", yLabs = ylabs, labelpos.y = minObs + lingap + extragap ))

(div.ph.its = plotReg(data.obs.its, x="ph", y="index_measure", minLimit = minObs, maxLimit = maxObs, gap = lingap, xLabs = "Soil pH", yLabs = ylabs, labelpos.y = minObs + lingap + 4))


# crop diversity
(div.cropdiv.16 = plotReg(data.obs.16s, 
                          x="totalCrops", y="index_measure", minLimit = minObs, maxLimit = maxObs, gap = lingap, xLabs = "Crop diversity", yLabs = ylabs, labelpos.y = minObs + lingap + extragap))

(div.cropdiv.its = plotReg(data.obs.its, 
                           x="totalCrops", y="index_measure", minLimit = minObs, maxLimit = maxObs, gap = lingap, xLabs = "Crop diversity", yLabs = ylabs, labelpos.y = minObs + 5))


library(patchwork)
library(ggpubr)

(obs16S.metam/obs16S.soil | (div.ph.16 / div.cropdiv.16 )) + plot_annotation(tag_levels = "A")

(obsITS.metam/obsITS.soil | (div.ph.its / div.cropdiv.its ))+ plot_annotation(tag_levels = "A")


shan16S = readRDS("./Field_analysis_output_v3/diversity/plot_shan16S_field.rds")
shanITS = readRDS("./Field_analysis_output_v3/diversity/plot_shanITS_field.rds")

ggarrange(shan16S, shanITS, labels = c("A", "B"), align = "hv", common.legend = T, legend = "bottom")

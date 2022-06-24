WDir="Your working directory"
setwd(WDir)

library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(tidyverse)


rm(list = ls())

# import data

# # # import biom files from qiime2
mydata_16s = import_biom("../source_data/mcm/16s/feature-table-tax.biom ")
tree.16s = read_tree("../source_data/mcm/16s/tree.nwk")
phy_tree(mydata_16s) = tree.16s
meta_16s = data.frame(sample_data(mydata_16s))

# # taxonomy header converted to bASV+number
tax_table(mydata_16s) = gsub("D_[0-9]__", "", tax_table(mydata_16s)) # remove prefix
colnames(tax_table(mydata_16s)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa_names(mydata_16s) = paste0("bASV", seq(1, nrow(otu_table(mydata_16s))))


pathogen_cols = c('Meloidogyne_hapla',
                  'Meloidogyne_chitwoodi',
                  'Pratylenchus_penetrans',
                  'Pratylenchus_neglectus',
                  'Pratylenchus_thornei',
                  'Paratrichodorus_sp',
                  'Tylenchorhynchus_sp',
                  'Paratylenchus_sp',
                  'Helicotylenchus_sp',
                  'Mesocriconema_sp',
                  'Ditylenchus_sp',
                  'Xiphinema_sp',
                  'Hemicycliophora_sp',
                  'vert_count',
                  'bd_count')

numeric_cols = c(pathogen_cols, "Period_numeric", "Soil_pH_Value")

colnames(meta_16s)[which(colnames(meta_16s) == "MesocriconemaÂ.sp")] ="Mesocriconema_sp"

# # reorder
#
cols.16s = colnames(meta_16s)
newcols.16s = c("SampleID", cols.16s)
meta_16s$SampleID = rownames(meta_16s)
meta_16s = meta_16s[newcols.16s]


# change pathogen count to numeric
for(i in numeric_cols){
  meta_16s[i] = lapply(meta_16s[i], as.numeric)
}

factor_cols = c("Period", "Soil_depth", "cul_type", "soil_compname","soil_taxclass","soil_taxgrtgroup","soil_taxorder","soil_taxsubgrp","soil_taxsuborder")

for(i in factor_cols){
  meta_16s[i] = lapply(meta_16s[i], as.factor)
}


# change field and microcosm naming
meta_16s$Field = gsub("Boardman", "B", meta_16s$Field)
meta_16s$Field = gsub("Organic", "O", meta_16s$Field)

meta_16s$microcosms = gsub("Boardman", "B", meta_16s$microcosms)
meta_16s$microcosms = gsub("Organic", "O", meta_16s$microcosms)


# all lower case column names
colnames(meta_16s) = tolower(colnames(meta_16s))

# calculate rotation years

yrs = paste0("y", seq(1998, 2018))

totalCrops = apply(meta_16s[, yrs], 1, function(x){
  
  x = as.vector(unlist(x))
  x[x==""] = NA
  # some rotation has - in them
  allcrops = unlist(str_split(x, "-"))
  uniqueCrops = unique(na.omit(allcrops))
  if(uniqueCrops[1] == "unknown"){
    return(-999)
  } else {
    return(length(uniqueCrops))
  }
})



totalRotYrs = apply(meta_16s[, yrs], 1, function(x){
  
  # return not ''
  x[x==""] = NA
  allyears = na.omit(unlist(x))
  if(unique(allyears)[1] == "unknown"){
    return(-999)
  }else{
    return(length(allyears))
  }
  
})


totalRotYrs_Pota = apply(meta_16s[, yrs], 1, function(x){
  
  x[x==""] = NA
  # some rotation has - in them
  allcrops = na.omit(unlist(str_split(x, "-")))
  if(unique(allcrops)[1] == "unknown"){
    return(-999)
  } else {
    return(sum(tolower(allcrops)=="potato"))
  }
  
})


percPota = totalRotYrs_Pota / totalRotYrs

meta_16s = cbind(meta_16s, totalCrops, totalRotYrs, totalRotYrs_Pota, percPota)

# new meta data
sample_data(mydata_16s) = meta_16s

saveRDS(mydata_16s, "./RDS_MCM_data_prunning/mydata16s_origin.rds")


###############################################################################
###############################################################################
###############################################################################

mydata_its = import_biom("../source_data/mcm/its/feature-table-tax.biom")
tree.its = read_tree("../source_data/mcm/its/tree.nwk")
phy_tree(mydata_its) =  tree.its
meta_its = data.frame(sample_data(mydata_its))

tax_table(mydata_its) = gsub("k__|p__|c__|o__|f__|g__|s__", "", tax_table(mydata_its))
colnames(tax_table(mydata_its)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

taxa_names(mydata_its) = paste0("fASV", seq(1, nrow(otu_table(mydata_its))))


# change to numbers or factors
pathogen_cols = c('Meloidogyne_hapla',
                  'Meloidogyne_chitwoodi',
                  'Pratylenchus_penetrans',
                  'Pratylenchus_neglectus',
                  'Pratylenchus_thornei',
                  'Paratrichodorus_sp',
                  'Tylenchorhynchus_sp',
                  'Paratylenchus_sp',
                  'Helicotylenchus_sp',
                  'Mesocriconema_sp',
                  'Ditylenchus_sp',
                  'Xiphinema_sp',
                  'Hemicycliophora_sp',
                  'vert_count',
                  'bd_count')
#
numeric_cols = c(pathogen_cols, "Period_numeric", "Soil_pH_Value")


colnames(meta_its)[which(colnames(meta_its) == "MesocriconemaÂ.sp")] = "Mesocriconema_sp"

# reorder put SampleID first
cols.its = colnames(meta_its)
newcols.its = c("SampleID", cols.its)
meta_its$SampleID = rownames(meta_its)
meta_its = meta_its[newcols.its]


# change those pathogen counts to numeric
for(i in numeric_cols){
  meta_its[i] = lapply(meta_its[i], as.numeric)
}

factor_cols = c("Period", "Soil_depth", "cul_type", "soil_compname","soil_taxclass","soil_taxgrtgroup","soil_taxorder","soil_taxsubgrp","soil_taxsuborder")

# change those to factors
for(i in factor_cols){
  meta_its[i] = lapply(meta_its[i], as.factor)
}


meta_its$Field = gsub("Boardman", "B", meta_its$Field)
meta_its$Field = gsub("Organic", "O", meta_its$Field)

meta_its$microcosms = gsub("Boardman", "B", meta_its$microcosms)
meta_its$microcosms = gsub("Organic", "O", meta_its$microcosms)

# all lower case column names, except for the first
colnames(meta_its) = tolower(colnames(meta_its))

totalCrops = apply(meta_its[, yrs], 1, function(x){
  
  x = as.vector(unlist(x))
  x[x==""] = NA
  # some rotation has - in them
  allcrops = unlist(str_split(x, "-"))
  uniqueCrops = unique(na.omit(allcrops))
  if(uniqueCrops[1] == "unknown"){
    return(-999)
  } else {
    return(length(uniqueCrops))
  }
})


totalRotYrs = apply(meta_its[, yrs], 1, function(x){
  
  # return not ''
  x[x==""] = NA
  allyears = na.omit(unlist(x))
  if(unique(allyears)[1] == "unknown"){
    return(-999)
  }else{
    return(length(allyears))
  }
  
})


totalRotYrs_Pota = apply(meta_its[, yrs], 1, function(x){
  
  x[x==""] = NA
  # some rotation has - in them
  allcrops = na.omit(unlist(str_split(x, "-")))
  if(unique(allcrops)[1] == "unknown"){
    return(-999)
  } else {
    return(sum(tolower(allcrops)=="potato"))
  }
  
})

percPota = totalRotYrs_Pota / totalRotYrs
meta_its = cbind(meta_its, totalCrops, totalRotYrs, totalRotYrs_Pota, percPota)

sample_data(mydata_its) = meta_its

# save to orignal files folder
saveRDS(mydata_its, "./RDS_MCM_data_prunning/mydataITS_origin.rds")

##############################################################################
# more pruning
##############################################################################

mydata_16s #124696 taxa 574 samples
mydata_its # 16417 taxa 576 samples

mydata_16s = readRDS("./RDS_MCM_data_prunning/mydata16s_origin.rds")
mydata_its = readRDS("./RDS_MCM_data_prunning/mydataITS_origin.rds")


sample_data(mydata_16s)


# total reads
sum(sample_sums(mydata_16s)) # 15,166,470
sum(sample_sums(mydata_its)) # 25,541,470


table(tax_table(mydata_16s)[, "Phylum"], exclude= F) # 1670 NA, 
table(tax_table(mydata_its)[, "Phylum"], exclude = F) # 5012 NA, 1380 unidentified


# remove NA, and rare species for downstream (some algorithm maybe sensitive to rare species)
# set at 10

mydata_16s.prune = subset_taxa(mydata_16s, !is.na(Phylum))
mydata_16s.prune = prune_taxa(taxa_sums(mydata_16s.prune)>10, mydata_16s.prune)

mydata_16s.prune
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 44593 taxa and 574 samples ]
# sample_data() Sample Data:       [ 574 samples by 62 sample variables ]
# tax_table()   Taxonomy Table:    [ 44593 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 44593 tips and 43957 internal nodes ]

mydata_its.prune = subset_taxa(mydata_its, !is.na(Phylum))
mydata_its.prune = prune_taxa(taxa_sums(mydata_its.prune)>10, mydata_its.prune)

mydata_its.prune
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6063 taxa and 576 samples ]
# sample_data() Sample Data:       [ 576 samples by 62 sample variables ]
# tax_table()   Taxonomy Table:    [ 6063 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 6063 tips and 5885 internal nodes ]

#---------------
# for 16s we will remove samples less than 7000 reads
# for ITS we will remove smaples less than 10000 reads

mydata_16s.prune = prune_samples(sample_sums(mydata_16s.prune) >= 7000, mydata_16s.prune)
min(sample_sums(mydata_16s.prune)) # 7019

mydata_its.prune = prune_samples(sample_sums(mydata_its.prune) >= 10000, mydata_its.prune)
min(sample_sums(mydata_its.prune)) # 10476

# use microbiomeutility to find the best taxonomy rank for NA
mbst16s = format_to_besthit(mydata_16s.prune)
mbstITS = format_to_besthit(mydata_its.prune)


mbst16s_taxtb = tax_table(mbst16s)[, 1:7]
rownames(mbst16s_taxtb) = gsub("OTU-|:.*$", "", rownames(mbst16s_taxtb))

mbstits_taxtb = tax_table(mbstITS)[, 1:7]
rownames(mbstits_taxtb) = gsub("OTU-|:.*$", "", rownames(mbstits_taxtb))

tax_table(mydata_16s.prune) = mbst16s_taxtb
tax_table(mydata_its.prune) = mbstits_taxtb

saveRDS(mydata_16s.prune, file.path("./RDS_MCM_data_prunning/", "mydata_16s_pruned.rds"))
saveRDS(mydata_its.prune, file.path("./RDS_MCM_data_prunning/", "mydata_its_pruned.rds"))

#--------------------------------------------------------------------------------
# more prunning
#--------------------------------------------------------------------------------

rm(list = ls())

mdata_16s.prune = readRDS("./RDS_MCM_data_prunning/mydata_16s_pruned.rds")
mdata_its.prune = readRDS("./RDS_MCM_data_prunning/mydata_its_pruned.rds")

# fumigation history
meta16S = data.frame(sample_data(mdata_16s.prune))
metaITS = data.frame(sample_data(mdata_its.prune))

meta16S$field
metaRate = c(122, 161, 0, 152, 190, 75, 120, 130, 0, 0, 0, 0, 0)
metaApp = c(3, 4, 0, 4, 5, 2, 3,3, 0, 0, 0, 0, 0)

names(metaRate) = unique(meta16S$field)
metaRate
names(metaApp) = unique(meta16S$field)
metaApp
all_metaRate = metaRate[meta16S$field]
all_metaRate
all_metaApp = metaApp[meta16S$field]
all_metaApp

meta16S$metamRates = all_metaRate
meta16S$metamApps = all_metaApp
meta16S$metamStatus = ifelse(meta16S$metamRates > 0, "fumigated", "Not-fumigated")
meta16S$metamStatus = as.factor(meta16S$metamStatus)

sample_data(mdata_16s.prune) = meta16S

all_metaRate = metaRate[metaITS$field]
all_metaRate
all_metaApp = metaApp[metaITS$field]
all_metaApp

metaITS$metamRates = all_metaRate
metaITS$metamApp = all_metaApp
metaITS$metamStatus = ifelse(metaITS$metamRates > 0, "fumigated", "Not-fumigated")
metaITS$metamStatus = as.factor(metaITS$metamStatus)



sample_data(mdata_its.prune) = metaITS

saveRDS(mdata_16s.prune, "./RDS_MCM_data_prunning/mydata_16s_pruned.rds")
saveRDS(mdata_its.prune, "./RDS_MCM_data_prunning/mydata_its_pruned.rds")

taxmtrx16S = as.matrix(tax_table(mdata_16s.prune))
taxmtrx16S[, 1:7]
taxmtrx16S[taxmtrx16S == "uncultured"] = NA
taxmtrx16S[taxmtrx16S == "uncultured soil bacterium"] = NA
taxmtrx16S[taxmtrx16S == "uncultured bacterium"] = NA
taxmtrx16S[taxmtrx16S == "metagenome"] = NA
taxmtrx16S[taxmtrx16S == "wastewater metagenome"] = NA
taxmtrx16S[taxmtrx16S == "groundwater metagenome"] = NA
taxmtrx16S[taxmtrx16S == "uncultured organism"] = NA
taxmtrx16S[taxmtrx16S == "uncultured sediment bacterium"] = NA
taxmtrx16S[taxmtrx16S == "bacterium enrichment culture clone auto84_4W"] = NA
taxmtrx16S[taxmtrx16S == "permafrost metagenome"] = NA
taxmtrx16S[taxmtrx16S == "uncultured forest bacterium"] = NA
taxmtrx16S[taxmtrx16S == "microbial mat metagenome"] = NA
taxmtrx16S[taxmtrx16S == "unidentified"] = NA
taxmtrx16S[taxmtrx16S == "g__uncultured"] = NA
taxmtrx16S[taxmtrx16S == "f__uncultured"] = NA
taxmtrx16S[taxmtrx16S == "p__uncultured"] = NA
taxmtrx16S[taxmtrx16S == "c__uncultured"] = NA
taxmtrx16S[taxmtrx16S == "f__uncultured"] = NA
unique(taxmtrx16S[,"Species"])
unique(taxmtrx16S[,"Domain"])

tax_table(mdata_16s.prune) = taxmtrx16S
mdata_16s.prune = subset_taxa(mdata_16s.prune, Domain != "Archaea")

unique(data.frame(tax_table(mdata_16s.prune))$Domain) # Bacteria

saveRDS(mdata_16s.prune, "./RDS_MCM_data_prunning/mdata16S_cleantaxa_updated.rds")
#
mdata16_bst = format_to_besthit(mdata_16s.prune)
taxmtrx16S.2 = as.matrix(tax_table(mdata16_bst))
taxmtrx16S.2[taxmtrx16S.2 == "d__uncultured"] = NA
taxmtrx16S.2[taxmtrx16S.2 == "p__uncultured"] = NA
taxmtrx16S.2[taxmtrx16S.2 == "c__uncultured"] = NA
taxmtrx16S.2[taxmtrx16S.2 == "o__uncultured"] = NA
taxmtrx16S.2[taxmtrx16S.2 == "f__uncultured"] = NA
taxmtrx16S.2[taxmtrx16S.2 == "g__uncultured"] = NA
taxmtrx16S.2[taxmtrx16S.2 == "s__uncultured"] = NA
#
taxmtrx16S.2[1:2, 1:5]
rownames(taxmtrx16S.2) = gsub(":.*", "", rownames(taxmtrx16S.2))
taxmtrx16S.2[1:2, 1:5]
tax_table(mdata_16s.prune) = taxmtrx16S.2[, -ncol(taxmtrx16S.2)]

mdata16_bst2 = format_to_besthit(mdata_16s.prune)

tax_table(mdata16_bst2)[1:100, c("Class", "Family")]

meta16S = data.frame(sample_data(mdata16_bst2))
#
meta16S$cul_type
meta16S$cul_type = gsub("Green", "Conventional", as.character(meta16S$cul_type))
#
fieldnames = sort(unique(meta16S$field))
fieldnames
fieldseq = paste0("f", seq(fieldnames))
fieldseq
names(fieldseq) = fieldnames
newcol = fieldseq[meta16S$field]
newcol
#
meta16S$field_alias = newcol
meta16S = meta16S[c(which(meta16S$cul_type == "Conventional"), which(meta16S$cul_type == "Organic")),]
meta16S$cul_type = as.factor(meta16S$cul_type)
sample_data(mdata16_bst2) = meta16S
saveRDS(mdata16_bst2, "./RDS_MCM_data_prunning/mdata16S_cleantaxa_bst.rds")

#-------++++++++++++++++++++++++++++++++++++++++++++++++++

taxmtrxITS = as.matrix(tax_table(mdata_its.prune))
unique(taxmtrxITS[, "Species"])
taxmtrxITS[taxmtrxITS == "uncultured"] = NA
taxmtrxITS[taxmtrxITS == "unidentified"] = NA
tax_table(mdata_its.prune) = taxmtrxITS

unique(data.frame(tax_table(mdata_its.prune))$Domain) # Fungi

saveRDS(mdata_its.prune, "./RDS_MCM_data_prunning/mdataITS_cleantaxa.rds")

mdataITS_bst = format_to_besthit(mdata_its.prune)

metaITS = data.frame(sample_data(mdataITS_bst))
metaITS$cul_type = gsub("Green", "Conventional", as.character(metaITS$cul_type))

fieldnames = sort(unique(metaITS$field))
fieldnames
fieldseq = paste0("f", seq(fieldnames))
fieldseq
names(fieldseq) = fieldnames
newcol = fieldseq[metaITS$field]
newcol

metaITS$field_alias = newcol

metaITS = metaITS[c(which(metaITS$cul_type == "Conventional"), which(metaITS$cul_type == "Organic")),]
metaITS$cul_type = as.factor(metaITS$cul_type)

sample_data(mdataITS_bst) = metaITS
saveRDS(mdataITS_bst, "./RDS_MCM_data_prunning/mdataITS_cleantaxa_bst.rds")


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

rm(list = ls())

mdata16S = readRDS("./RDS_MCM_data_prunning/mdata16S_cleantaxa_bst.rds")
mdataITS = readRDS("./RDS_MCM_data_prunning/mdataITS_cleantaxa_bst.rds")


mdata16S
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 44549 taxa and 555 samples ]
# sample_data() Sample Data:       [ 555 samples by 66 sample variables ]
# tax_table()   Taxonomy Table:    [ 44549 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 44549 tips and 43913 internal nodes ]

mdataITS
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6063 taxa and 567 samples ]
# sample_data() Sample Data:       [ 567 samples by 66 sample variables ]
# tax_table()   Taxonomy Table:    [ 6063 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 6063 tips and 5885 internal nodes ]


meta16S = data.frame(sample_data(mdata16S))
str(meta16S)
meta16S$metamStatus
#?relevel
meta16S$metamStatus = relevel(meta16S$metamStatus, "Not-fumigated")
meta16S$metamStatus
levels(meta16S$metamStatus) = c("Not fumigated", "Fumigated")
sample_data(mdata16S) = meta16S

rank_names(mdata16S)
# mdata16S@otu_table
taxa_names(mdata16S) = gsub(":.*", "", rownames(mdata16S@otu_table))
# mdata16S@otu_table
mdata16S@tax_table = mdata16S@tax_table[, -ncol(mdata16S@tax_table)]
colnames(mdata16S@tax_table) = c("Kingdom", "Phylum",   "Class",    "Order",    "Family",   "Genus",    "Species")

saveRDS(mdata16S, "./RDS_MCM_data_prunning/mdata16S_cleantaxa_bst2.rds")


metaITS = data.frame(sample_data(mdataITS))
str(metaITS)
metaITS$metamStatus

metaITS$metamStatus = relevel(metaITS$metamStatus, "Not-fumigated")
metaITS$metamStatus
levels(metaITS$metamStatus) = c("Not fumigated", "Fumigated")
sample_data(mdataITS) = metaITS

rank_names(mdataITS)
# mdata16S@otu_table
taxa_names(mdataITS) = gsub(":.*", "", rownames(mdataITS@otu_table))
# mdata16S@otu_table
mdataITS@tax_table = mdataITS@tax_table[, -ncol(mdataITS@tax_table)]
colnames(mdataITS@tax_table) = c("Kingdom", "Phylum",   "Class",    "Order",    "Family",   "Genus",    "Species")


saveRDS(mdataITS, "./RDS_MCM_data_prunning/mdataITS_cleantaxa_bst2.rds")

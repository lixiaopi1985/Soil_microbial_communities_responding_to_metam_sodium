WDir="Your working directory"
setwd(WDir)


library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(tidyverse)

# import data from QIIME2 output

#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
# 16S
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

fdata_16s = import_biom("../source_data/onsite/16s/feature-table-tax.biom")
tree.16s = read_tree("../source_data/onsite/16s/tree.nwk")
phy_tree(fdata_16s) = tree.16s
meta_16s = data.frame(sample_data(fdata_16s))

# change the taxonomy rank
tax_table(fdata_16s) = gsub("D_[0-9]__", "", tax_table(fdata_16s)) # remove prefix to make it look better
colnames(tax_table(fdata_16s)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # fixed the rank names

numeric_cols = c('Meloidogyne_hapla',	'Meloidogyne_chitwoodi',	'Pratylenchus_penetrans',	'Pratylenchus_neglectus',	'Pratylenchus_thornei',	'Paratrichodorus_sp',	'Tylenchorhynchus_sp','Paratylenchus_sp',	'Helicotylenchus_sp','Mesocriconema_sp',	'Ditylenchus_sp',	'Xiphinema_sp',	'Hemicycliophora_sp',	'vert_count',	'bd_count', 'pH')


colnames(meta_16s)[which(colnames(meta_16s) == "MesocriconemaÂ.sp")]="Mesocriconema_sp"

# change pathogen count to numeric
for(i in numeric_cols){
  meta_16s[i] = lapply(meta_16s[i], as.numeric)
}

factor_cols = c("cul_type", "soil_compname","soil_taxclass","soil_taxgrtgroup","soil_taxorder","soil_taxsubgrp","soil_taxsuborder", paste0("y", seq(1998, 2018)))

for(i in factor_cols){
  meta_16s[i] = lapply(meta_16s[i], as.factor)
}

# colum names to lower case
colnames(meta_16s) = tolower(colnames(meta_16s))
cols_16s = colnames(meta_16s)
meta_16s$SampleID = row.names(meta_16s)
meta_16s = meta_16s[c("SampleID", cols_16s)]
sample_data(fdata_16s) = meta_16s
# add bASV to indicate bacterial community
taxa_names(fdata_16s) = paste0("bASV", seq_along(taxa_names(fdata_16s)))

fdata_16s

saveRDS(fdata_16s, "./RDS_Field_data_prunning/fdata16S_origin.rds")

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 38974 taxa and 48 samples ]
# sample_data() Sample Data:       [ 48 samples by 53 sample variables ]
# tax_table()   Taxonomy Table:    [ 38974 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 38974 tips and 38653 internal nodes ]


#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
# Fungi
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

fdata_its = import_biom("../source_data/onsite/its/feature-table-tax.biom")
tree.its = read_tree("../source_data/onsite/its/tree.nwk")
phy_tree(fdata_its) = tree.its

meta_its = data.frame(sample_data(fdata_its))

tax_table(fdata_its) = gsub("k__|p__|c__|o__|f__|g__|s__", "", tax_table(fdata_its))
colnames(tax_table(fdata_its)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

numeric_cols = c('Meloidogyne_hapla',	'Meloidogyne_chitwoodi',	'Pratylenchus_penetrans',	'Pratylenchus_neglectus',	'Pratylenchus_thornei',	'Paratrichodorus_sp',	'Tylenchorhynchus_sp','Paratylenchus_sp',	'Helicotylenchus_sp','Mesocriconema_sp',	'Ditylenchus_sp',	'Xiphinema_sp',	'Hemicycliophora_sp',	'vert_count',	'bd_count', 'pH')


colnames(meta_its)[which(colnames(meta_its) == "MesocriconemaÂ.sp")]="Mesocriconema_sp"

# change pathogen count to numeric
for(i in numeric_cols){
  meta_its[i] = lapply(meta_its[i], as.numeric)
}

factor_cols = c("cul_type", "soil_compname","soil_taxclass","soil_taxgrtgroup","soil_taxorder","soil_taxsubgrp","soil_taxsuborder", paste0("y", seq(1998, 2018)))

for(i in factor_cols){
  meta_its[i] = lapply(meta_its[i], as.factor)
}


colnames(meta_its) = tolower(colnames(meta_its))
cols_its = colnames(meta_its)
meta_its$SampleID = row.names(meta_its)
meta_its = meta_its[c("SampleID", cols_its)]
sample_data(fdata_its) = meta_its
taxa_names(fdata_its) = paste0("fASV", seq_along(taxa_names(fdata_its)))

fdata_its

saveRDS(fdata_its, "./RDS_Field_data_prunning/fdataITS_origin.rds")

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7127 taxa and 48 samples ]
# sample_data() Sample Data:       [ 48 samples by 53 sample variables ]
# tax_table()   Taxonomy Table:    [ 7127 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 7127 tips and 7006 internal nodes ]

#----------------------------------------------------------------------------

##################################################
# Get summary of the reads and reads distribution
##################################################

# total counts

sum(sample_sums(fdata_16s)) # 3724679
sum(sample_sums(fdata_its)) # 7043461

# min depth
min(sample_sums(fdata_16s)) # min 40026
min(sample_sums(fdata_its)) # min 31323


# rare taxa
min(taxa_sums(fdata_16s)) # min 1
min(taxa_sums(fdata_its)) # min 1

# prune data
# remove NA phylum

fdata_16s.prune = subset_taxa(fdata_16s, !is.na(Phylum))
fdata_its.prune = subset_taxa(fdata_its, !is.na(Phylum))

# remove singleton and doubleton and OTU abundance less than 10 across all samples

fdata_16s.prune = prune_taxa(taxa_sums(fdata_16s.prune) > 10, fdata_16s.prune)
fdata_its.prune = prune_taxa(taxa_sums(fdata_its.prune) > 10, fdata_its.prune)


# add crop rotation data
meta_16s.prune = as(sample_data(fdata_16s.prune), "data.frame")
meta_its.prune = as(sample_data(fdata_its.prune), "data.frame")

crops.16s = meta_16s.prune %>%
  select(SampleID, field, "y1998":"y2018") %>%
  mutate_all(list(~na_if(., ""))) %>%
  gather(key="rotYear", value="crops", "y1998":"y2018") %>%
  group_by(field) %>%
  mutate(totalCrops = n_distinct(crops, na.rm = T)) %>%
  mutate(totalRotYrs = n_distinct(rotYear[!is.na(crops)]), totalYrs_Pota = sum(crops=="POTATO", na.rm = T)/n_distinct(SampleID)) %>%
  mutate(percPota = totalYrs_Pota / totalRotYrs ) %>%
  spread(rotYear, crops)


new_meta.16s.prune = merge(meta_16s.prune, crops.16s[, c("SampleID", "totalCrops", "totalRotYrs", "totalYrs_Pota", "percPota")], by = "SampleID")

new_meta.its.prune = merge(meta_its.prune, crops.16s[, c("SampleID", "totalCrops", "totalRotYrs", "totalYrs_Pota", "percPota")], by = "SampleID")

rownames(new_meta.16s.prune) = new_meta.16s.prune$SampleID
sample_data(fdata_16s.prune) = new_meta.16s.prune

rownames(new_meta.its.prune) = new_meta.its.prune$SampleID
sample_data(fdata_its.prune) = new_meta.its.prune

# generate pathogen data
path_meta = new_meta.16s.prune %>%
  select(-c("field_coord", "coll_coord", "cords_note", "quant_reading", "sample_or_control", "soil_collected_date", "mean_day_air_temp"))

pathop = file.path("./Pathogen_data/", "meta_pathogens.tsv")
write.table(path_meta, file=pathop, sep = "\t", row.names = F)


# use microbiomeutility to find the best taxonomy rank for NA
fbst16s = format_to_besthit(fdata_16s.prune)
fbstITS = format_to_besthit(fdata_its.prune)

fbst16s_taxtb = tax_table(fbst16s)[, 1:7]
fbst16s_taxtb
rownames(fbst16s_taxtb) = gsub("OTU-|:.*$", "", rownames(fbst16s_taxtb))

fbstits_taxtb = tax_table(fbstITS)[, 1:7]
rownames(fbstits_taxtb) = gsub("OTU-|:.*$", "", rownames(fbstits_taxtb))


tax_table(fdata_16s.prune) = fbst16s_taxtb
tax_table(fdata_its.prune) = fbstits_taxtb


saveRDS(fdata_16s.prune, file.path("./RDS_Field_data_prunning/", "fdata_16s_pruned.rds"))
saveRDS(fdata_its.prune, file.path("./RDS_Field_data_prunning/", "fdata_its_pruned.rds"))


fdata_16s.prune
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 21915 taxa and 48 samples ]
# sample_data() Sample Data:       [ 48 samples by 57 sample variables ]
# tax_table()   Taxonomy Table:    [ 21915 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 21915 tips and 21657 internal nodes ]

fdata_its.prune
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 3028 taxa and 48 samples ]
# sample_data() Sample Data:       [ 48 samples by 57 sample variables ]
# tax_table()   Taxonomy Table:    [ 3028 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 3028 tips and 2985 internal nodes ]

#===============================================================================
# more prunning
#===============================================================================

rm(list = ls())

fdata_16s.prune = readRDS("./RDS_Field_data_prunning/fdata_16s_pruned.rds")
fdata_its.prune = readRDS("./RDS_Field_data_prunning/fdata_its_pruned.rds")


#fumigation history
meta16S = data.frame(sample_data(fdata_16s.prune))
metaITS = data.frame(sample_data(fdata_its.prune))

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

all_metaRate = metaRate[metaITS$field]
all_metaRate
all_metaApp = metaApp[metaITS$field]
all_metaApp

metaITS$metamRates = all_metaRate
metaITS$metamApp = all_metaApp

sample_data(fdata_16s.prune) = meta16S
sample_data(fdata_its.prune) = metaITS

saveRDS(fdata_16s.prune, "./RDS_Field_data_prunning/fdata_16s_pruned.rds")
saveRDS(fdata_its.prune, "./RDS_Field_data_prunning/fdata_its_pruned.rds")

taxmtrx16S = as.matrix(tax_table(fdata_16s.prune))
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
unique(taxmtrx16S[,"Species"])

tax_table(fdata_16s.prune) = taxmtrx16S
fdata_16s.prune
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 21915 taxa and 48 samples ]
# sample_data() Sample Data:       [ 48 samples by 59 sample variables ]
# tax_table()   Taxonomy Table:    [ 21915 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 21915 tips and 21657 internal nodes ]

unique(data.frame(tax_table(fdata_16s.prune))$Domain)

fdata_16s.prune = subset_taxa(fdata_16s.prune, Domain != "Archaea")
fdata_16s.prune # 16 taxa were Archaea
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 21899 taxa and 48 samples ]
# sample_data() Sample Data:       [ 48 samples by 59 sample variables ]
# tax_table()   Taxonomy Table:    [ 21899 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 21899 tips and 21641 internal nodes ]

saveRDS(fdata_16s.prune, "./RDS_Field_data_prunning/fdata16S_cleantaxa.rds")

fdata16_bst = format_to_besthit(fdata_16s.prune)

taxmtrx16S.2 = as.matrix(tax_table(fdata16_bst))
taxmtrx16S.2[taxmtrx16S.2 == "d__uncultured"] = NA
taxmtrx16S.2[taxmtrx16S.2 == "p__uncultured"] = NA
taxmtrx16S.2[taxmtrx16S.2 == "c__uncultured"] = NA
taxmtrx16S.2[taxmtrx16S.2 == "o__uncultured"] = NA
taxmtrx16S.2[taxmtrx16S.2 == "f__uncultured"] = NA
taxmtrx16S.2[taxmtrx16S.2 == "g__uncultured"] = NA
taxmtrx16S.2[taxmtrx16S.2 == "s__uncultured"] = NA

taxmtrx16S.2[1:2, 1:5]
rownames(taxmtrx16S.2) = gsub(":.*", "", rownames(taxmtrx16S.2))
taxmtrx16S.2[1:2, 1:5]
tax_table(fdata_16s.prune) = taxmtrx16S.2[, -ncol(taxmtrx16S.2)]

fdata16_bst2 = format_to_besthit(fdata_16s.prune)

meta16S = data.frame(sample_data(fdata16_bst2))
meta16S$cul_type
meta16S$cul_type = gsub("Green", "Conventional", as.character(meta16S$cul_type))

fieldnames = unique(meta16S$field)
fieldnames
fieldseq = paste0("f", seq(fieldnames))
fieldseq
names(fieldseq) = fieldnames
newcol = fieldseq[meta16S$field]
newcol

meta16S$field_alias = newcol
meta16S = meta16S[c(which(meta16S$cul_type == "Conventional"), which(meta16S$cul_type == "Organic")),]
meta16S$cul_type = as.factor(meta16S$cul_type)
meta16S$metamStatus = ifelse(meta16S$metamRates > 0, "Fumigated", "No-fumigated")
meta16S$metamStatus = as.factor(meta16S$metamStatus)

sample_data(fdata16_bst2) = meta16S
fdata16_bst2

saveRDS(fdata16_bst2, "./RDS_Field_data_prunning/fdata16S_cleantaxa_bst.rds")

###################################################################
#
###################################################################


taxmtrxITS = as.matrix(tax_table(fdata_its.prune))
unique(taxmtrxITS[, "Species"])
taxmtrxITS[taxmtrxITS == "uncultured"] = NA
taxmtrxITS[taxmtrxITS == "unidentified"] = NA

tax_table(fdata_its.prune) = taxmtrxITS
fdata_its.prune
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 3028 taxa and 48 samples ]
# sample_data() Sample Data:       [ 48 samples by 59 sample variables ]
# tax_table()   Taxonomy Table:    [ 3028 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 3028 tips and 2985 internal nodes ]


saveRDS(fdata_its.prune, "./RDS_Field_data_prunning/fdataITS_cleantaxa.rds")

fdataITS_bst = format_to_besthit(fdata_its.prune)

taxmtrxITS.2 = as.matrix(tax_table(fdataITS_bst))
unique(taxmtrxITS.2[,"Domain"])

taxmtrxITS.2[1:2,1:2]
taxmtrxITS.2[taxmtrxITS.2 == "d__uncultured"] = NA
taxmtrxITS.2[taxmtrxITS.2 == "p__uncultured"] = NA
taxmtrxITS.2[taxmtrxITS.2 == "c__uncultured"] = NA
taxmtrxITS.2[taxmtrxITS.2 == "o__uncultured"] = NA
taxmtrxITS.2[taxmtrxITS.2 == "f__uncultured"] = NA
taxmtrxITS.2[taxmtrxITS.2 == "g__uncultured"] = NA
taxmtrxITS.2[taxmtrxITS.2 == "s__uncultured"] = NA

taxmtrxITS.2[1:2, 1:5]
rownames(taxmtrxITS.2) = gsub(":.*", "", rownames(taxmtrxITS.2))
taxmtrxITS.2[1:2, 1:5]
tax_table(fdata_its.prune) = taxmtrxITS.2[, -ncol(taxmtrxITS.2)]

fdataITS_bst2 = format_to_besthit(fdata_its.prune)

metaITS = data.frame(sample_data(fdataITS_bst2))
metaITS$cul_type
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
metaITS$metamStatus = ifelse(metaITS$metamRates > 0, "Fumigated", "No-fumigated")
metaITS$metamStatus = as.factor(metaITS$metamStatus)
sample_data(fdataITS_bst2) = metaITS

saveRDS(fdataITS_bst2, "./RDS_Field_data_prunning/fdataITS_cleantaxa_bst.rds")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

rm(list = ls())
fdata16S = readRDS("./RDS_Field_data_prunning/fdata16S_cleantaxa_bst.rds")
fdataITS = readRDS("./RDS_Field_data_prunning/fdataITS_cleantaxa_bst.rds")

meta16S = data.frame(sample_data(fdata16S))
metaITS = data.frame(sample_data(fdataITS))

meta16S$metamStatus

meta16S$metamStatus = relevel(meta16S$metamStatus, "No-fumigated")
meta16S$metamStatus
levels(meta16S$metamStatus) = c("Not fumigated", "Fumigated")
sample_data(fdata16S) = meta16S
meta16S$metamStatus

rank_names(fdata16S)
#fdata16S@otu_table
taxa_names(fdata16S) = gsub(":.*", "", rownames(fdata16S@otu_table))
rownames(fdata16S@otu_table)
fdata16S@tax_table = fdata16S@tax_table[, -ncol(fdata16S@tax_table)]
colnames(fdata16S@tax_table) = c("Kingdom", "Phylum",   "Class",    "Order",    "Family",   "Genus",    "Species")

saveRDS(fdata16S, "./RDS_Field_data_prunning/fdata16S_cleantaxa_bst2.rds")

metaITS$metamStatus
metaITS$metamStatus = relevel(metaITS$metamStatus, "No-fumigated")
metaITS$metamStatus
levels(metaITS$metamStatus) = c("Not fumigated", "Fumigated")
sample_data(fdataITS) = metaITS

rank_names(fdataITS)

taxa_names(fdataITS) = gsub(":.*", "", rownames(fdataITS@otu_table))

fdataITS@tax_table = fdataITS@tax_table[, -ncol(fdataITS@tax_table)]
colnames(fdataITS@tax_table) = c("Kingdom", "Phylum",   "Class",    "Order",    "Family",   "Genus",    "Species")

saveRDS(fdataITS, "./RDS_Field_data_prunning/fdataITS_cleantaxa_bst2.rds")

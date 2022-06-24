WDir="Your working directory"
setwd(WDir)


library(NetCoMi)
library(phyloseq)
library(RColorBrewer)
rm(list = ls())


fdata16S.genus = readRDS("./RDSdata/fdata16S_genus.rds")
fdataITS.genus = readRDS("./RDSdata/fdataITS_genus.rds")


levels(phyloseq::get_variable(fdata16S.genus, "metamStatus"))



f16S.metam = subset_samples(fdata16S.genus, metamStatus=="Fumigated")
taxtab.metam = f16S.metam@tax_table@.Data

# assign OTU names
rownames(f16S.metam@otu_table@.Data) = taxtab.metam[, "Genus"]

f16S.no = subset_samples(fdata16S.genus, metamStatus=="Not fumigated")
taxtab.no = f16S.no@tax_table@.Data
rownames(f16S.no@otu_table@.Data) = taxtab.no[, "Genus"]

net16S.spie = netConstruct(data = f16S.metam,
                       data2 = f16S.no,
                       measure = "spieceasi",
                       measurePar = list(lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50)),
                       filtTax = "highestFreq",
                       filtTaxPar = list(highestFreq=80),
                       zeroMethod = "multRepl",
                       normMethod = "clr",
                       sparsMethod = "threshold",
                       cores = parallel::detectCores()/2,
                       thresh = 0.6,
                       verbose = 3,
                       seed = 123456)

# saveRDS(net16S.spie, "./Field_analysis_output_v3/network/net16S_spiecEasi.rds")

netAna16S = netAnalyze(net16S.spie,
                       centrLCC = TRUE,
                       hubPar = c("degree", "between", "closeness", "eigenvector"),
                       clustMethod = "cluster_fast_greedy",
                       lnormFit = T,
                       normDeg = T,
                       normBetw = T,
                       normClose = T,
                       normEigen = T)

# saveRDS(netAna16S, "./Field_analysis_output_v3/network/netAnalysis16S.rds")

netAna16S = readRDS("./Field_analysis_output_v3/network/netAnalysis16S.rds")


plot(netAna16S,
     layout = "spring",
     sameLayout = T,
     layoutGroup = "union",
     repulsion = 0.95,
     rmSingles = "inboth",
     posCol = "darkturquoise", 
     negCol = "orange",
     edgeTranspLow = 0,
     edgeTranspHigh = 40,
     labels = g16S,
     shortenLabels = "none",
     charToRm = "g__",
     labelScale = F,
     nodeSize = "degree",
     nodeSizeSpread = 4,
     nodeColor = "cluster",
     nodeTransp = 10,
     cexNodes = 5,
     cexLabels = 0.9,
     cexHubLabels = 1,
     cexTitle = 1.2,
     groupNames = c("Fumigated", "Not fumigated"),
     showTitle = T,
     hubBorderCol = "black",
     hubBorderWidth = 5, 
     mar=c(2,3,2,3))
# mar=c(4,6,4,6))

legend("bottom", cex = 1, title = "estimated correlation:",
       legend = c("+","-"), inset=0.1,lty = 1, lwd = 2, col = c("darkturquoise","orange"),
       bty = "n", horiz = TRUE)


#----------------------------------------------------------------------------------------
# fungi

levels(phyloseq::get_variable(fdataITS.genus, "metamStatus"))

fITS.metam = subset_samples(fdataITS.genus, metamStatus=="Fumigated")
taxtab.metam = fITS.metam@tax_table@.Data

# assign OTU names
rownames(fITS.metam@otu_table@.Data) = taxtab.metam[, "Genus"]


fITS.no = subset_samples(fdataITS.genus, metamStatus=="Not fumigated")
taxtab.no = fITS.no@tax_table@.Data
rownames(fITS.no@otu_table@.Data) = taxtab.no[, "Genus"]

netITS.spie = netConstruct(data = fITS.metam,
                       data2 = fITS.no,
                       measure = "spieceasi",
                       measurePar = list(lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50)),
                       filtTax = "highestFreq",
                       filtTaxPar = list(highestFreq=80),
                       zeroMethod = "multRepl",
                       normMethod = "clr",
                       sparsMethod = "threshold",
                       cores = parallel::detectCores()/2,
                       thresh = 0.6,
                       verbose = 3,
                       seed = 123456)

#saveRDS(netITS.spie, "./Field_analysis_output_v3/network/netITS_spiecEasi.rds")

netAnaITS = netAnalyze(netITS.spie,
                       centrLCC = TRUE,
                       hubPar = c("degree", "between", "closeness", "eigenvector"),
                       clustMethod = "cluster_fast_greedy",
                       lnormFit = T,
                       normDeg = T,
                       normBetw = T,
                       normClose = T,
                       normEigen = T)

#saveRDS(netAnaITS, "./Field_analysis_output_v3/network/netAnalysisITS.rds")

plot(netAnaITS,
     layout = "spring",
     sameLayout = T,
     layoutGroup = "union",
     repulsion = 0.95,
     rmSingles = "inboth",
     posCol = "darkturquoise", 
     negCol = "orange",
     edgeTranspLow = 0,
     edgeTranspHigh = 40,
     shortenLabels = "none",
     charToRm = "g__",
     labelScale = F,
     nodeSize = "degree",
     nodeSizeSpread = 4,
     nodeColor = "cluster",
     nodeTransp = 10,
     cexNodes = 5,
     cexLabels = 0.9,
     cexHubLabels = 1,
     cexTitle = 1.2,
     groupNames = c("Fumigated", "Not fumigated"),
     showTitle = T,
     hubBorderCol = "black",
     hubBorderWidth = 5, 
     mar=c(2,3,2,3))
# mar=c(4,6,4,6))

legend("bottom", cex = 1, title = "estimated correlation:",
       legend = c("+","-"), inset=0.1,lty = 1, lwd = 2, col = c("darkturquoise","orange"),
       bty = "n", horiz = TRUE)

#-------------------------------------------------------------------------------------

rm(list = ls())
netAna16S = readRDS("./Field_analysis_output_v3/network/netAnalysis16S.rds")
netAnaITS = readRDS("./Field_analysis_output_v3/network/netAnalysisITS.rds")

length(unique(netAna16S$clustering$clust1)) #11
table(netAna16S$clustering$clust1) # cluster 1, 2 have 9 and 7 genera
# 0  1  2  3  4  5  6  7  8  9 10 
# 23  9  7  4  5  4  3  3  3  2  2 
netAna16S$clustering$clust1[netAna16S$clustering$clust1 == 1]
netAna16S$clustering$clust1[netAna16S$clustering$clust1 == 2]

length(unique(netAna16S$clustering$clust2)) #12
table(netAna16S$clustering$clust2) # cluster 1, 2, 3 have 9,11, 7 genera
# 0  1  2  3  4  5  6  7  8  9 10 11 
# 11  9 11  7  5  2  5  5  2  2  4  2 
netAna16S$clustering$clust2[netAna16S$clustering$clust2 == 2]

netAna16S$hubs

netAna16S$centralities
netAna16S$input$assoMat1

#$$$$$$$$$$$$$$$$$$$$$$$
length(unique(netAnaITS$clustering$clust1)) #10
length(unique(netAnaITS$clustering$clust2)) #8

table(netAnaITS$clustering$clust1)
netAnaITS$clustering$clust1[netAnaITS$clustering$clust1 == 1]

table(netAnaITS$clustering$clust2)
netAnaITS$clustering$clust2[netAnaITS$clustering$clust2 == 1]

#------------------------------------------------------------------------
netCompare16S = netCompare(netAna16S,
                        permTest = T,
                        storeAssoPerm = T,
                        fileStoreAssoPerm = "./Field_analysis_output_v3/network/f16S_assoPerm_nperm1000_updated",
                        nPerm = 1000,
                        cores = 2,
                        seed = 123)
saveRDS(netCompare16S, "./Field_analysis_output_v3/network/netCompare16S.rds")

netCompareITS = netCompare(netAnaITS,
                           permTest = T,
                           storeAssoPerm = T,
                           fileStoreAssoPerm = "./Field_analysis_output_v3/network/ITS_assoPerm_nperm1000_updated",
                           nPerm = 1000,
                           cores = 2,
                           seed = 123)
saveRDS(netCompareITS, "./Field_analysis_output_v3/network/netCompareITS.rds")

#------------------------------------------------------------------------------------------------
# check stats
rm(list=ls())


netCompare16S = readRDS("./Field_analysis_output_v3/network/netCompare16S.rds")
netCompareITS = readRDS("./Field_analysis_output_v3/network/netCompareITS.rds")


summary(netCompare16S)
summary(netCompareITS)

netAna16S = readRDS("./Field_analysis_output_v3/network/netAnalysis16S.rds")
netAnaITS = readRDS("./Field_analysis_output_v3/network/netAnalysisITS.rds")

summary(netAna16S)

netAna16S$clustering$clust1[netAna16S$clustering$clust1 == 1]
netAna16S$clustering$clust2[netAna16S$clustering$clust2 == 2]

table(netAnaITS$clustering$clust1)
table(netAnaITS$clustering$clust2)

netAnaITS$clustering$clust1[netAnaITS$clustering$clust1 == 1]
netAnaITS$clustering$clust2[netAnaITS$clustering$clust2 == 1]

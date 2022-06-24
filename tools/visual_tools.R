plotHeat_matrix = function(psmelt, taxlevel="Genus", groupby = c("period", "cul_type"), ntop=20, nshow=10, printlevel=0){
  # make heatmap of this (top 10?)
  library(tidyr)
  
  top20 = psmelt %>%
    group_by(!!!syms(groupby)) %>%
    mutate(total_abund = sum(Abundance)) %>%
    mutate(new_taxa =  !!sym(taxlevel)) %>%
    group_by(!!!syms(c(groupby, "new_taxa"))) %>%
    mutate(rel_abund = 100*Abundance/total_abund) %>%# change cul_type to groupby
    summarise(accum_abund = sum(rel_abund)) %>%
    slice_max(order_by = accum_abund, n=ntop)
  #----------------------------------------------------------------------


  # overall top abundant genera
  domdf = psmelt %>%
    mutate(total_abund = sum(Abundance), rel_abund = 100*Abundance/total_abund) %>%
    mutate(new_taxa = !!sym(taxlevel)) %>%
    group_by(new_taxa) %>%
    mutate(rel_abund = 100*Abundance/total_abund) %>%
    summarise(accum_abund = sum(rel_abund)) %>%
    arrange(desc(accum_abund))
  
  doms = domdf$new_taxa[1:nshow] #ordered
  
  # some top genus might be not in other group, so find them as well
  top_accum = psmelt %>%
    group_by(!!!syms(groupby)) %>%
    mutate(total_abund = sum(Abundance), rel_abund = 100*Abundance/total_abund) %>%
    mutate(new_taxa = !!sym(taxlevel)) %>%
    group_by(!!!syms(c(groupby, "new_taxa"))) %>%
    mutate(rel_abund = 100*Abundance/total_abund) %>%
    summarise(accum_abund = sum(rel_abund)) %>%
    filter(new_taxa %in% unique(top20$new_taxa)) %>%
    arrange(desc(accum_abund))

  top.meta.2 = top_accum %>%
    unite(col=Samples, !!!syms(groupby)) %>%
    spread(Samples, accum_abund) %>%
    as.data.frame()
  
  rownames(top.meta.2) = top.meta.2$new_taxa
  
  # matrix
  top.matrix = as.matrix(top.meta.2[, -1])
  
  if(printlevel==0){
    return(list(m = top.matrix, topbygroup=top_accum, dom = doms, overal=domdf))
  } else if(printlevel == 1){
    
    prinL.df = psmelt %>%
      dplyr::select(Phylum, !!sym(taxlevel)) %>%
      mutate(taxa_cbind = paste(Phylum, !!sym(taxlevel), sep = "-")) %>%
      filter(!!sym(taxlevel) %in% doms) %>%
      dplyr::distinct()
    
    # order it
    new_doms = prinL.df$taxa_cbind[match(doms, prinL.df[[taxlevel]])]
    
    # keep the order right
    
    return(list(m = top.matrix, topbygroup = top_accum, dom = new_doms, overal=domdf))
  }
  
  
}

#---------------------------------------------------------------------------------------------------

# plotHeat_plot = function(tops, COLOR=NA, new_break = NA,fontsize=12, out=T, outfile="heatmap.png", clusterRows=F, clusterCols=F, Legend=T, anno_legend=T, custom_anno=T, annoInput=NA, annoColor=NA, cellw = 40, cellh = 40, ...){
#   
#   library(stringr)
#   library(RColorBrewer)
#   library(pheatmap)
#   
#   top.matrix = tops$m
#   doms = tops$dom
#   
#   org_colname = colnames(top.matrix)
#   
#   color_set = brewer.pal(9, "BuGn")
#   colorTime = colorRampPalette(color_set)
#   
#   if(custom_anno){
#     if(!is.na(annoInput)){
#       
#       top.matrix.ord2 = top.matrix[doms,]
#       anno_colors = annoColor
#       
#       if(out){
#         g = pheatmap::pheatmap(top.matrix.ord2, 
#                                annotation_col = annoInput,
#                                annotation_colors = anno_colors, 
#                                cellwidth = cellw, 
#                                cellheight = cellh, 
#                                color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(length(new_break)), 
#                                cluster_cols = clusterCols, 
#                                cluster_rows = clusterRows, 
#                                breaks = new_break, 
#                                fontsize = fontsize, 
#                                filename = outfile, 
#                                legend = Legend,
#                                annotation_legend = anno_legend, 
#                                # column_split = annoCol$Fumigation_history, 
#                                ...)
#         
#         
#       } else {
#         g = pheatmap::pheatmap(top.matrix.ord2, 
#                                annotation_col = annoCol,
#                                cellwidth = cellw, 
#                                cellheight = cellh,
#                                # color = colorRampPalette(c("navy", "white", "firebrick3"))(length(new_break)),
#                                color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(length(new_break)), 
#                                annotation_colors = anno_colors, 
#                                cluster_cols = clusterCols, 
#                                cluster_rows = clusterRows, 
#                                breaks = new_break, 
#                                fontsize = fontsize, 
#                                # filename = outfile, 
#                                legend = Legend,
#                                annotation_legend = anno_legend, 
#                                # column_split = annoCol$Fumigation_history, 
#                                ...)
#       }
#       
#     } else {
#       stop("Please provide annotation")
#     }
#     
#   } else {
#     
#     # no annotation
#     top.matrix.ord2 = top.matrix[dom, ]
#     
#     if(out){
#       g = pheatmap::pheatmap(top.matrix.ord2, 
#                              # annotation_col = annoCol,
#                              cellwidth = cellw, 
#                              cellheight = cellh, 
#                              color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(length(new_break)), 
#                              # annotation_colors = anno_colors, 
#                              cluster_cols = clusterCols, 
#                              cluster_rows = clusterRows, 
#                              breaks = new_break, 
#                              fontsize = fontsize, 
#                              filename = outfile, 
#                              legend = Legend,
#                              annotation_legend = anno_legend,
#                              # column_split = annoCol$Fumigation_history, 
#                              ...)
#     } else {
#       g = pheatmap::pheatmap(top.matrix.ord2, 
#                              # annotation_col = annoCol,
#                              cellwidth = cellw, 
#                              cellheight = cellh, 
#                              color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(length(new_break)), 
#                              # annotation_colors = anno_colors, 
#                              cluster_cols = clusterCols, 
#                              cluster_rows = clusterRows, 
#                              breaks = new_break, 
#                              fontsize = fontsize, 
#                              # filename = outfile, 
#                              legend = Legend,
#                              annotation_legend = anno_legend,
#                              # column_split = annoCol$Fumigation_history, 
#                              ...)
#     }
#    
#   }
#   
#   return(g)
# }


# add significant * to the row labels
addSig = function(origlabel, cond1, cond2){
  
  # italicized the organisms
  new_labs = c()
  for(i in origlabel){
    if((i %in% cond1) & (i %in% cond2)){
      sigLab = paste0(i, " */*")
      new_labs = c(new_labs, sigLab)
    } else if((i %in% cond1) & !(i %in% cond2)){
      sigLab = paste0(i, " */")
      new_labs = c(new_labs, sigLab)
    } else if(!(i %in% cond1) & (i %in% cond2)){
      sigLab = paste0(i, " /*")
      new_labs = c(new_labs, sigLab)
    } else {
      sigLab = i
      new_labs = c(new_labs, sigLab)
    }
  }
  
  return(new_labs)
}

scientific_name_formatter <- function(raw_name, plotdevice="ggplot") { 
  library(stringr)
  # containing p__, c__, o__, f__, g__, s__
  # name contain numbers and other 
  
  if(plotdevice=="ggplot"){
    formatted_names = sapply(raw_name, function(x){
      
      print(x)
      nsplit = str_split(x, "__")
      
      if(length(nsplit[[1]]) == 2){
        # if the second name contains numbers, -, length > 1, all capitalized do not italicize
        if(grepl("[0-9]+|-|uncultured", nsplit[[1]][2])){
          spname = paste0(nsplit[[1]][1], "__", nsplit[[1]][2])
        } else {
          spname = paste0(nsplit[[1]][1], "__*", nsplit[[1]][2], "*")
        }
        
      } else if(length(nsplit[[1]] == 1)) {
        
        if(grepl("[0-9]+|-|uncultured", nsplit[[1]][1])){
          spname = nsplit[[1]][1]
        } else {
          spname = paste0("*", nsplit[[1]][1], "*")
        }
      }
      return(spname)
    })
    
    
  } else {
    # other format
    
    formatted_names = sapply(raw_name, function(x){
      
      nsplit = str_split(x, "__")
      print(x)
      
      if(length(nsplit[[1]]) == 2){

        if(grepl("[0-9]+|-|uncultured", nsplit[[1]][2])){
          
          # if the second name contains numbers, -, length > 1, all capitalized do not italicize
          spname = toString(paste0(nsplit[[1]][1], "__", nsplit[[1]][2]))
          
          
          
        } else {
          
          # italicize
          # to extract "*/*" from f__XXXXX */*
          
          if(grepl("\\*/$|\\*/\\*$|/\\*$", nsplit[[1]][2])){
            
            # separate out "/*" from the string first
            nsplit_space = str_split(nsplit[[1]][2], " ")
            
            if(length(nsplit_space[[1]]) == 2){
              # contains "/*" etc
              spname = bquote(.(nsplit[[1]][1]) * "__" * italic(.(nsplit_space[[1]][1])) ~ .(nsplit_space[[1]][2]))
              
              
            } else if(length(nsplit_space[[1]]) == 1) {
              # does not contain "/*
              spname = bquote(.(nsplit[[1]][1]) * "__" * italic(.(nsplit_space[[1]][1])))
              
            }
          } else {
            # if f__XXXX does not contain "*/", "*/*", "/*"
            spname = bquote(.(nsplit[[1]][1]) * "__" * italic(.(nsplit[[1]][2])))
            
          }


        }
        
      } else if(length(nsplit[[1]]) == 1) {
        print("Length 1")
        if(grepl("[0-9]+|-|uncultured", nsplit[[1]][1])){
          spname = nsplit[[1]][1]
        } else {
          
          # situation like XXXX */
          if(grepl("\\*/$|\\*/\\*$|/\\*$", nsplit[[1]][1])){
            
            # separate out "/*" from the string first
            nsplit_space = str_split(nsplit[[1]][1], " ")
            
            if(length(nsplit_space[[1]]) == 2){
              # contains "/*" etc
              spname = bquote(italic(.(nsplit_space[[1]][1])) ~ .(nsplit_space[[1]][2]))
              
            } else if(length(nsplit_space[[1]]) == 1) {
              # does not contain "/*
              spname = bquote(italic(.(nsplit[[1]][1])))
            }
          } else {
            # if f__XXXX does not contain "*/", "*/*", "/*"
            spname = bquote(italic(.(nsplit[[1]][1])))
          }
        }
      }
      return(spname)
    })
  }

  
  return(unname(formatted_names))
}


plotBars = function(glomdf, taxlevel="Genus", groupby = c("Season", "time_label", "funlab"), ntop=20, nshow=10, customColors = NULL, choosecolorset = NULL, otherColor = "grey45", plotx = "funlab", facet_group = "~ Season + time_label", xlab="Fungicides", ylab="Relative abundance (%)", labfill="", plottitle="", xlabmargin_t=10, ylabmargin_r=10, fontsize=12, printlevels=0, fontangle=0, fontHjust=1, fontVjust=0.5, DFonly=F, plotMargin_t = 2, plotMargin_r = 2, plotMargin_b = 2, plotMargin_l = 2){
  
  library(ggplot2)
  library(ggh4x)
  
  top = plotHeat_matrix(glomdf, taxlevel = taxlevel, groupby = groupby, ntop = ntop, nshow = nshow, printlevel = printlevels)

  doms = top$dom[!is.na(top$dom)]
  
  cat("Top taxa selected:\n")
  print(paste(doms, sep = "\n"))

  if(printlevels==0){
    
    Togroup = c(groupby, taxlevel)
    
    # sp. not italic, p__, o__, c__, f__, g__, c__Subgroup 6, c__KD4-96 not italic
    topn = glomdf %>%
      group_by(!!!syms(Togroup)) %>%
      summarise(accum_abund = sum(Abundance)) %>%
      mutate(new_taxa = ifelse(!!sym(taxlevel) %in% doms, !!sym(taxlevel), "Others"))%>%
      mutate(new_taxa2 = ifelse(new_taxa == "Others", "Others", scientific_name_formatter(new_taxa))) %>% # format the name
      mutate(new_taxa_itali = ifelse(new_taxa2 == "Others", 
                                     "Others",
                                     ifelse(grepl(" sp.*", new_taxa2),
                                            paste0(paste0(gsub(" sp.*$", "", new_taxa2), "*"), " spp."),
                                            new_taxa2)))
      
      # mutate(new_taxa_itali = ifelse(new_taxa2 == "Others", "Others",
      #                                ifelse(grepl(" sp.", new_taxa2),  
      #                                       paste(paste0("*", gsub(" sp.$", "", new_taxa2), "*"), " spp."), 
      #                                       new_taxa2)))
      
    
    # if new_taxa == other, no change, than look for sp in the species names
    

    topn$new_taxa = factor(topn$new_taxa, levels = c(doms, "Others"))
    
    # changed it to italic
    doms_italic = scientific_name_formatter(doms)
    

    dom_level = c(ifelse(grepl(" sp.*", doms_italic),
                       paste(paste0(gsub(" sp.*$", "", doms_italic), "*"), " spp."),
                       doms_italic), "Others")

    
    # dom_level = c(ifelse(grepl(" sp.", doms), 
    #                    paste(paste0("*", gsub(" sp.$", "", doms), "*"), " spp."), 
    #                    paste0("*", doms, "*")), "Others")

    topn$new_taxa_itali = factor(topn$new_taxa_itali, levels = dom_level)

    
  } else if(printlevels == 1){
    # add Phylum to it
    Togroup = c(groupby, "taxa_cbind", taxlevel)

    topn = glomdf %>%
      mutate(taxa_cbind = paste(Phylum, !!sym(taxlevel), sep = "-")) %>%
      group_by(!!!syms(Togroup)) %>%
      summarise(accum_abund = sum(Abundance)) %>%
      mutate(new_taxa = ifelse(taxa_cbind %in% doms, taxa_cbind, "Others")) %>%
      mutate(new_taxa2 = ifelse(new_taxa == "Others", "Others", scientific_name_formatter(new_taxa))) %>%
      mutate(new_taxa_itali = ifelse(new_taxa2 == "Others", "Others",
                                     ifelse(grepl(" sp.*", new_taxa2),
                                            paste0(paste0(gsub(" sp.*$", "", new_taxa2), "*"), " spp."),
                                     new_taxa2)))
      
      # mutate(new_taxa_itali = ifelse(new_taxa == "Others", "Others",
      #                                ifelse(grepl(" sp.", new_taxa),  
      #                                       paste(paste0("*", gsub(" sp.$", "", new_taxa), "*"), " spp."), 
      #                                       paste0("*", new_taxa, "*")))) 
      #mutate(new_taxa_itali = ifelse(new_taxa != "Others", paste0("*", new_taxa, "*"), new_taxa))
    #!!!!!!!!!!!!!!!!!!! new_taxa_itali has not tested

    # order new_taxa
    topn$new_taxa = factor(topn$new_taxa, levels = c(doms, "Others"))
    #topn$new_taxa_itali = factor(topn$new_taxa_itali, levels = c(paste0("*", doms, "*"), "Others"))
    doms_italic = scientific_name_formatter(doms)
    dom_level = c(ifelse(grepl(" sp.*", doms_italic),
                         paste(paste0(gsub(" sp.*$", "", doms_italic), "*"), " spp."),
                         doms_italic), "Others")
    
    # dom_level = c(ifelse(grepl(" sp.", doms), 
    #                      paste(paste0("*", gsub(" sp.$", "", doms), "*"), " spp."), 
    #                      paste0("*", doms, "*")), "Others")
    
    topn$new_taxa_itali = factor(topn$new_taxa_itali, levels = dom_level)
  }
  
  
  if(is.null(customColors)){
    if(is.null(choosecolorset)){
      colors = ggsci::pal_npg()(10)
    } else {
      colors = choosecolorset
    }
    
    colorPal = colorRampPalette(colors)(length(unique(topn$new_taxa_itali)))
    names(colorPal) = levels(topn$new_taxa_itali)
    colorPal["Others"] = otherColor
    
  } else {
    colorPal = customColors
  }


  
  if(is.null(labfill)){
    labfill = NULL 
  } else if(labfill==""){
    labfill = taxlevel
  }
  
  
  if(DFonly){
    return(topn)
  }
  
  if(facet_group == ""){

    g = topn %>%
      ggplot(aes(x = !!sym(plotx), y=accum_abund, fill=new_taxa_itali))+
      geom_col(position = position_fill()) +
      scale_fill_manual(values = colorPal) +
      # scale_y_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0), labels = scales::percent_format(accuracy =1, suffix = "")) +
      theme_classic() +
      labs(x=xlab, y=ylab, fill=labfill, title=plottitle) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.text = ggtext::element_markdown(),
        text = element_text(size=fontsize, color="black"),
        #plot.margin = margin(t=plotMargin_t, r=plotMargin_r, b=plotMargin_b, l=plotMargin_l),
        legend.margin = margin(t=plotMargin_t, r=plotMargin_r, b=plotMargin_b, l=plotMargin_l),
        axis.text.x = element_text(size=fontsize, color="black", angle = fontangle, hjust = fontHjust, vjust = fontVjust),
        axis.text.y = element_text(size=fontsize, color="black"),
        axis.title.x = element_text(size=fontsize, margin = margin(t=xlabmargin_t), color="black"),
        axis.title.y = element_text(size=fontsize, margin = margin(r=ylabmargin_r), color="black"))
  } else {
    g = topn %>%
      ggplot(aes(x = !!sym(plotx), y=accum_abund, fill=new_taxa_itali))+
      geom_col(position = position_fill()) +
      facet_nested(formula(facet_group), scales = "free") +
      scale_fill_manual(values = colorPal) +
      # scale_y_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0), labels = scales::percent_format(accuracy =1, suffix = "")) +
      theme_classic() +
      labs(x=xlab, y=ylab, fill=labfill, title = plottitle) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.text = ggtext::element_markdown(),
        text = element_text(size=fontsize, color="black"),
        #plot.margin = margin(t=plotMargin_t, r=plotMargin_r, b=plotMargin_b, l=plotMargin_l),
        legend.margin = margin(t=plotMargin_t, r=plotMargin_r, b=plotMargin_b, l=plotMargin_l),
        axis.text.x = element_text(size=fontsize, color="black",angle = fontangle, hjust = fontHjust, vjust = fontVjust),
        axis.text.y = element_text(size=fontsize, color="black"),
        axis.title.x = element_text(size=fontsize, margin = margin(t=xlabmargin_t), color="black"),
        axis.title.y = element_text(size=fontsize, margin = margin(r=ylabmargin_r), color="black"))
  }
  return(g)
}


plotORDS = function(phylo, tranform_method = "hellinger", ordmethod="PCoA", distmethod="bray", fs=12, plottype="sample", colorby="funlab", colorLab = "Fungicides", fontcolor="black", alpha=1, ncolor = 4, colorPal = ggsci::pal_npg(), Ellipse=F, ell_type="t", ell_linetype=2, returnDist = F, gtitle=""){

  
  # ?transform
  phylo.norm = transform(phylo, transform = tranform_method)
  
  phylo.ord = ordinate(phylo.norm, method = ordmethod, distmethod)

  
  g =  plot_ordination(phylo.norm, phylo.ord, type=plottype, color=colorby) +
    geom_vline(aes(xintercept=0, ), color = "grey80") +
    geom_hline(aes(yintercept=0), color = "grey80") +
    geom_point(size=5, alpha=alpha) +
    scale_color_manual(values = colorPal(ncolor)) +
    theme_bw() +
    labs(color=colorLab, title=gtitle) +
    theme(text = element_text(size=fs, color= fontcolor),
          axis.title.x = element_text(size=fs, color = fontcolor, margin=margin(t=10)),
          axis.title.y = element_text(size=fs, color=fontcolor, margin=margin(r=10)),
          axis.text.x = element_text(size=fs, color=fontcolor),
          axis.text.y = element_text(size=fs, color=fontcolor))
  
  if(Ellipse){
    g = g+stat_ellipse(type=ell_type, linetype=ell_linetype)
  }
  
  if(returnDist){
    return(phyloseq::distance(phylo.norm, method = distmethod))
  } else {
    return(g)
  }

}


DEanalysis = function(phylo, model, fittype="local", Contrast=NULL, alpha=0.05, abs_log2cutoff=1, Colns = c("log2FoldChange", "padj", "Phylum", "Class", "Genus", "group")){
  # create a significant df
  
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  
  diagdds = phyloseq_to_deseq2(phylo, model)
  
  geoMeans = apply(counts(diagdds), 1, gm_mean)
  diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
  diagdds = DESeq(diagdds, fitType=fittype)
  print(resultsNames(diagdds))
  res = results(diagdds, contrast = Contrast)
  summary(res)
  res = res[order(res$padj, na.last = NA),]
  
  sigtab = res[(res$padj < alpha),]
  
  if(nrow(sigtab)>0){
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phylo)[rownames(sigtab), ], "matrix"))
    sigtab = sigtab[(abs(sigtab$log2FoldChange)>=abs_log2cutoff), ]
    sigtab$group = paste(Contrast[2], Contrast[3], sep="_vs_")
    sigtab = sigtab[, Colns]
    return(sigtab)
  } else {
    return(data.frame())
  }

}


plotDE = function(df, topn = 10, xlab="log2 Fold Change", ylab = "Taxa (Phylum-Genus)", labfill = "Abundance"){
  
  # df = rbind(...)
  # 
  # df_spread = df %>%
  #   dplyr::mutate(taxa=paste(Phylum, Genus, sep="-")) %>%
  #   dplyr::select(taxa, log2FoldChange, group) %>%
  #   tidyr::spread(key = group, value = log2FoldChange)
  # 
  # df_shared = df_spread[complete.cases(df_spread),]
  # 
  # df_gather = df_shared %>%
  #   tidyr::gather(key=group, value = log2FoldChange, colnames(df_shared)[2:length(colnames(df_shared))])

  g = df %>%
    dplyr::mutate(taxa=paste(Phylum, Genus, sep="-")) %>%
    dplyr::select(taxa, log2FoldChange, group) %>%
    mutate(RAchanges = ifelse(log2FoldChange < 0, "Decreased", "Increased")) %>%
    mutate(RAchanges = factor(RAchanges, levels = c("Increased", "Decreased"))) %>%
    dplyr::select(log2FoldChange, taxa, group, RAchanges) %>%
    group_by(RAchanges) %>%
    slice_max(order_by = abs(log2FoldChange), n=topn) %>%
    ggplot(aes(x=log2FoldChange, y = reorder(taxa, log2FoldChange), fill=RAchanges)) +
    geom_bar(stat = "identity") +
    labs(x=xlab, y=ylab, fill=labfill) +
    theme_classic() +
    scale_fill_manual(values = c("#B03A2E", "#1B4F72")) +
    theme(text=element_text(size=12),
              axis.title.x = element_text(size=12, margin=margin(t=10), color="black"),
              axis.title.y = element_text(size=12, margin=margin(r=10), color="black"),
              axis.text.x = element_text(size=12, color="black"),
              axis.text.y = element_text(size=12, color="black"))
  
  return(g)
}


plotANCOM = function(microM_ancom_out, searchtop=NA, ntop=10, out=T, prefix="ANCOM_", savepath="./", xlab="W", ylab="Taxa", grp1="Group1", grp2="Group2", labfills="Enrich groups", setFactorLevels=F, levels=""){
  
  library(ggpubr)
  library(ggplot2)
  #output barplot and pieplot
  mrtb = marker_table(microM_ancom_out)
  mrtb.df = data.frame(mrtb)
  
  if(setFactorLevels){
    mrtb.df$enrich_group = factor(mrtb.df$enrich_group, levels=levels)
  }
  
  titlegroup = mrtb.df$enrich_group[order(mrtb.df$enrich_group, decreasing = T)]
  
  if(is.na(unique(titlegroup)[1]) ||  is.na(unique(titlegroup)[2])){
    warning("One of group does not exist")
    if(is.na(unique(titlegroup)[1])){
      print("Group1 is NA")
      grp1 = grp1
      comparison = paste0(unique(titlegroup)[1], "_vs_", grp1)
    } else if(is.na(unique(titlegroup)[2])) {
      print("Group2 is NA")
      grp1 = grp2
      comparison = paste0(grp2, "_vs_", unique(titlegroup)[1])
    } else {
      comparison = paste0(grp2, "_vs_", grp1)
    }
    
    
  } else {
    comparison = paste0(unique(titlegroup)[2], "_vs_", unique(titlegroup)[1])
  }
  
  
  if(!is.na(searchtop)){
    if(length(searchtop) < 0){
      warning("Please input search top vector")
    }
    
    mrtb.df$predom = ifelse(mrtb.df$feature %in% searchtop, paste0("*", mrtb.df$feature), mrtb.df$feature)
    mrtb.df$Ifpredom = ifelse(mrtb.df$feature %in% searchtop, "Y", "N")
    # mrtb.df$fontface = ifelse(mrtb.df$Ifpredom == "Y", "bold", "plain")
  }

  if(out){
    filename = paste0(prefix, "_", comparison, ".csv") 
    write.csv(mrtb.df, paste(savepath, filename, sep="/"), row.names = F, quote = F)
  }

  if("predom" %in% colnames(mrtb.df)){
    bar = mrtb.df %>%
      group_by(enrich_group) %>%
      slice_max(order_by = abs(ef_W), n = ntop) %>%
      ggplot(aes(x=ef_W, y=reorder(predom, ef_W), fill=enrich_group)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values = c("#2874A6", "#EC7063")) +
      theme_classic2() +
      labs(title = comparison, x=xlab, y=ylab, fill=labfills) +
      theme(text=element_text(size=12),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(size=12, margin=margin(t=10), color="black"),
            axis.title.y = element_text(size=12, margin=margin(r=10), color="black"),
            axis.text.x = element_text(size=12, color="black"),
            axis.text.y = element_text(size=12, color="black"))
  } else {
    bar = mrtb.df %>%
      group_by(enrich_group) %>%
      slice_max(order_by = abs(ef_W), n = ntop) %>%
      ggplot(aes(x=ef_W, y=reorder(feature, ef_W), fill=enrich_group)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values = c("#2874A6", "#EC7063")) +
      theme_classic2() +
      labs(title = comparison, x=xlab, y=ylab, fill=labfills) +
      theme(text=element_text(size=12),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(size=12, margin=margin(t=10), color="black"),
            axis.title.y = element_text(size=12, margin=margin(r=10), color="black"),
            axis.text.x = element_text(size=12, color="black"),
            axis.text.y = element_text(size=12, color="black"))
  }

  return(list(barp = bar, sig_features = mrtb.df$feature, df=mrtb.df))
}



plotPIE = function(phylo_df, mrtb.df, pie_cutoff=0.01, provideColor = NA, defaultColor="Paired", grp1="Group1", grp2="Group2"){
  
  
  library(RColorBrewer)
  
  getMeta = phylo_df %>%
    dplyr::select(Phylum, Genus) %>%
    dplyr::distinct() %>%
    filter(Genus %in% mrtb.df$feature) %>%
    dplyr::select(Phylum) %>%
    dplyr::group_by(Phylum) %>%
    dplyr::count(sort=T) %>%
    ungroup() %>%
    mutate(perc = n/sum(n)) %>%
    mutate(Phylum2 = ifelse(perc<pie_cutoff, paste0("<", 100*pie_cutoff, "%"), Phylum)) %>%
    group_by(Phylum2) %>%
    summarise(n2=sum(n)) %>%
    mutate(total = sum(n2), perc = n2/total, labels=scales::percent(perc, accuracy = 1)) 
  
  
  if(is.na(provideColor)){

    phylumColors = unique(getMeta$Phylum2)
    phylumColors2 = phylumColors[phylumColors != paste0("<", 100*pie_cutoff, "%")]
    phylumColors2 = c(phylumColors2, paste0("<", 100*pie_cutoff, "%"))
    cols = colorRampPalette(brewer.pal(9, defaultColor))(length(phylumColors2))
    names(cols) = phylumColors2

  } else {
    
    # select used colors
    cols = provideColor[names(provideColor) %in% unique(getMeta$Phylum2)]
  }

  # deal with plot title
  titlegroup = mrtb.df$enrich_group[order(mrtb.df$enrich_group, decreasing = T)]
  
  if(is.na(unique(titlegroup)[1]) ||  is.na(unique(titlegroup)[2])){
    warning("One of group does not exist")
    if(is.na(unique(titlegroup)[1])){
      print("Group1 is NA")
      grp1 = grp1
      comparison = paste0(unique(titlegroup)[1], "_vs_", grp1)
    } else if(is.na(unique(titlegroup)[2])) {
      print("Group2 is NA")
      grp1 = grp2
      comparison = paste0(grp2, "_vs_", unique(titlegroup)[1])
    } else {
      comparison = paste0(grp2, "_vs_", grp1)
    }
    
    
  } else {
    comparison = paste0(unique(titlegroup)[2], "_vs_", unique(titlegroup)[1])
  }
  
  
  pie = ggpie(getMeta, x="n2", label = "labels", fill="Phylum2", palette = cols, color = "white", legend="right") + labs(fill="Phylum", title = comparison) + theme(
    plot.title = element_text(hjust = 0.5)
  )
    

  return(list(df=getMeta, pie=pie))
}



plotDE2 = function(..., topn = 10, xlab="log2 Fold Change", ylab = "Taxa (Phylum-Genus)", labfill = "Abundance"){
  
  df = rbind(...)

  df_spread = df %>%
    dplyr::mutate(taxa=paste(Phylum, Genus, sep="-")) %>%
    dplyr::select(taxa, log2FoldChange, group) %>%
    tidyr::spread(key = group, value = log2FoldChange)

  df_shared = df_spread[complete.cases(df_spread),]

  df_gather = df_shared %>%
    tidyr::gather(key=group, value = log2FoldChange, colnames(df_shared)[2:length(colnames(df_shared))])
  
  g = df_gather %>%
    dplyr::select(taxa, log2FoldChange, group) %>%
    mutate(RAchanges = ifelse(log2FoldChange < 0, "Decreased", "Increased")) %>%
    mutate(RAchanges = factor(RAchanges, levels = c("Increased", "Decreased"))) %>%
    dplyr::select(log2FoldChange, taxa, group, RAchanges) %>%
    group_by(RAchanges) %>%
    slice_max(order_by = abs(log2FoldChange), n=topn) %>%
    ggplot(aes(x=log2FoldChange, y = reorder(taxa, log2FoldChange), fill=RAchanges)) +
    geom_bar(stat = "identity") +
    facet_grid(~ group) +
    labs(x=xlab, y=ylab, fill=labfill) +
    theme_classic() +
    scale_fill_manual(values = c("#B03A2E", "#1B4F72")) +
    theme(text=element_text(size=12),
          panel.grid.major.y = element_line(),
          axis.title.x = element_text(size=12, margin=margin(t=10), color="black"),
          axis.title.y = element_text(size=12, margin=margin(r=10), color="black"),
          axis.text.x = element_text(size=12, color="black"),
          axis.text.y = element_text(size=12, color="black"))
  
  return(g)
}

plotCoef_mcm = function(result_diffTest, phylo, fdr = 0.05, taxLevel="Genus", nlevel = 3,orderList = NA){
  
  # specific for mcm
  all_fdr_values = result_diffTest$p_fdr[!is.na(result_diffTest$p_fdr)]
  sigTaxa = corncob::otu_to_taxonomy(result_diffTest$significant_taxa, data=phylo, level=taxLevel)
  
  try_species = names(all_fdr_values)
  taxaNames = corncob::otu_to_taxonomy(try_species, data = phylo, level = taxLevel) # OTU name - genus 
  
  all_models = result_diffTest[["all_models"]][!is.na(result_diffTest[["all_models"]])]
  names(all_models) = taxaNames


  if(!is.na(orderList)){
    
    if((length(orderList) == length(all_models)) | (length(orderList) == length(taxaNames))){
      all_models = all_models[as.character(orderList)]
      taxaNames = taxaNames[order(match(taxaNames, orderList))]
      all_fdr_values = all_fdr_values[names(taxaNames)]
    } else {
      stop("orderList length is different to the model or taxaNames")
    }
    
  }
  
  
  var_levels = c()
  
  coef = c()
  maxse = c()
  minse = c()
  
  vari = c()
  vari_maxse = c()
  vari_minse = c()
  taxaNames_rep = c()
  fdr_rep = c()
  sigTaxa_w_symbol = c()
  
  
  for(i in 1:length(all_models)){
    
    print(i)
    print(paste("Order list:", orderList[i]))
    print(paste("Model name:", names(all_models[i])))
    
    if(taxaNames[i] %in% sigTaxa){
      sigTaxa_w_symbol = c(sigTaxa_w_symbol, rep("*", nlevel))
    } else {
      sigTaxa_w_symbol = c(sigTaxa_w_symbol, rep("ns", nlevel))
    }
    
    taxaNames_rep = c(taxaNames_rep, rep(taxaNames[i], nlevel))
    fdr_rep = c(fdr_rep, rep(all_fdr_values[i], nlevel))
    
    # how many coefficients
    n.mu = length(grep("mu", rownames(all_models[[i]]$coefficients)))
    n.phi = length(grep("phi", rownames(all_models[[i]]$coefficients)))

    
    if(n.mu > 1 & n.phi > 1){
      
      print("Mu and phi both tested")
      
      var_levels = c(var_levels, gsub("mu.", "", names(all_models[[i]]$coefficients[2:n.mu,1])))
      
      # coefficients
      coef = c(coef, all_models[[i]]$coefficients[2:n.mu,1]) 
      maxse = c(maxse, all_models[[i]]$coefficients[2:n.mu,1]+qnorm(0.975)*all_models[[i]]$coefficients[2:n.mu,2])
      minse = c(minse, all_models[[i]]$coefficients[2:n.mu,1]-qnorm(0.975)*all_models[[i]]$coefficients[2:n.mu,2])
      # coefficients pvalu
      
      # variability
      vari = c(vari, all_models[[i]]$coefficients[6:8,1])
      vari_maxse = c(vari_maxse, all_models[[i]]$coefficients[6:8,1]+qnorm(0.975)*all_models[[i]]$coefficients[6:8,2])
      vari_minse = c(vari_minse, all_models[[i]]$coefficients[6:8,1]-qnorm(0.975)*all_models[[i]]$coefficients[6:8,2])
      
      df = data.frame(coef = coef, facet_var = var_levels, variability = vari, xmin = minse, xmax = maxse, var_xmin = vari_minse, var_xmax = vari_maxse, taxa = taxaNames_rep, siglevels=sigTaxa_w_symbol, fdr = fdr_rep)

      
    } else if (n.mu > 1 & n.phi == 1){
      
      
      print("Only Mu tested")
      var_levels = c(var_levels, gsub("mu.", "", names(all_models[[i]]$coefficients[2:n.mu,1])))
      
      # coefficients
      coef = c(coef, all_models[[i]]$coefficients[2:n.mu,1]) 
      maxse = c(maxse, all_models[[i]]$coefficients[2:n.mu,1]+qnorm(0.975)*all_models[[i]]$coefficients[2:n.mu,2])
      minse = c(minse, all_models[[i]]$coefficients[2:n.mu,1]-qnorm(0.975)*all_models[[i]]$coefficients[2:n.mu,2])
      # coefficients pvalu
      
      # variability
      # vari = c(vari, all_models[[i]]$coefficients[6:8,1])
      # vari_maxse = c(vari_maxse, all_models[[i]]$coefficients[6:8,1]+qnorm(0.975)*all_models[[i]]$coefficients[6:8,2])
      # vari_minse = c(vari_minse, all_models[[i]]$coefficients[6:8,1]-qnorm(0.975)*all_models[[i]]$coefficients[6:8,2])
      
      df = data.frame(coef = coef, facet_var = var_levels, xmin = minse, xmax = maxse, taxa = taxaNames_rep, siglevels=sigTaxa_w_symbol, fdr = fdr_rep)
      

    } else if(n.mu == 1) {

      stop("No test")
    }
    }
    
  return(list(DF=df, mod=all_models))
}



plot_BCAsp = function(sp_df, speciesname, formulas, maxIT=25, SN="ES", padjust="tukey", rds_save = F, rds_path = "", posi_dodge=0.9, alpha=0.9, err_dodge_width=0.2, addColor="none", ncolor=4, fontcolor="black", fontsize=12, labfill="Fungicides", xlab="Sampling time (day)", ylab="Square root of abundance", setRange = NULL, sameYmax=F, nudgeY=0){
  

  cat("\n-------------------------------------\n")
  cat(SN, "\n")
  cat(speciesname)
  cat("\n-------------------------------------\n")
  # this function is hard coded
  select_sp = sp_df %>% filter(Species == speciesname & Season == SN)
  
  # 0 distribution
  print(scales::percent(sum(select_sp$Abundance == 0)/length(select_sp$Abundance)), digits = 4)
  
  # negative binomial
  model = glm.nb(formula(formulas), data=select_sp, control = glm.control(maxit = maxIT))
  
  # comparision within each time
  em = emmeans(model, pairwise ~ funlab | time_label, p.adjust= padjust)
  # cld.em = cld(em, Letters = letters, alpha=0.05, sort=T)
  cld.em = cld(em, Letters = letters)
  cld_df = as.data.frame(cld.em)

  # cld_df$group = gsub("3", "c", gsub("2", "b", gsub("1", "a", stringr::str_trim(cld_df$.group))))
  # 
  cld_df$joint = paste(cld_df$time_label, cld_df$funlab, sep = "_")
  select_sp$joint = paste(select_sp$time_label, select_sp$funlab, sep = "_")
  
  (select_sp.agg = select_sp %>%
      group_by(joint, time_label, funlab) %>%
      mutate(sqrtAbund = sqrt(Abundance+1)) %>%
      summarise(avg = mean(sqrtAbund), se = plotrix::std.error(sqrtAbund)))
  
  plotsp.df.ltrs = merge(select_sp.agg, cld_df, by="joint", no.dups = T)
  
  if(is.null(setRange)){
    rangeMax = sqrt(max(sp_df %>% filter(Species == speciesname) %>% pull(Abundance)))
  } else {
    rangeMax =setRange
  }
  
  
  if(addColor=="none"){
    colorpal =ggsci::pal_jco()(ncolor)
  } else {
    colorpal = addColor
  }
  
  if(sameYmax){
    plotsp_plot = plotsp.df.ltrs %>%
      ggplot(aes(x = time_label.x, y=avg, fill=funlab.x)) +
      geom_bar(stat="identity", position = position_dodge(posi_dodge), alpha=alpha, color="black") +
      geom_text(aes(label=.group, y=rangeMax), position = position_dodge(posi_dodge)) +
      geom_errorbar(aes(ymin=avg-se, ymax = avg+se), width=err_dodge_width, position = position_dodge(posi_dodge)) +
      # scale_y_continuous(limits = c(0, rangeMax)) +
      scale_fill_manual(values = colorpal) +
      theme_classic() +
      labs(fill = labfill, x = xlab, y= ylab, title=speciesname) +
      theme(text=element_text(size=12),
            axis.title.x = element_text(size=fontsize, margin=margin(t=10), color=fontcolor),
            axis.title.y = element_text(size=fontsize, margin=margin(r=10), color=fontcolor),
            axis.text.x =  element_text(size=fontsize, color=fontcolor),
            axis.text.y =  element_text(size=fontsize, color=fontcolor))
  } else {
    plotsp_plot = plotsp.df.ltrs %>%
      ggplot(aes(x = time_label.x, y=avg, fill=funlab.x)) +
      geom_bar(stat="identity", position = position_dodge(posi_dodge), alpha=alpha, color="black") +
      geom_text(aes(label=.group, y=avg+se+nudgeY), position = position_dodge(posi_dodge)) +
      geom_errorbar(aes(ymin=avg-se, ymax = avg+se), width=err_dodge_width, position = position_dodge(posi_dodge)) +
      # scale_y_continuous(limits = c(0, rangeMax)) +
      scale_fill_manual(values = colorpal) +
      theme_classic() +
      labs(fill = labfill, x = xlab, y= ylab, title=speciesname) +
      theme(text=element_text(size=12),
            axis.title.x = element_text(size=fontsize, margin=margin(t=10), color=fontcolor),
            axis.title.y = element_text(size=fontsize, margin=margin(r=10), color=fontcolor),
            axis.text.x =  element_text(size=fontsize, color=fontcolor),
            axis.text.y =  element_text(size=fontsize, color=fontcolor))
  }

  
  if(rds_save){
    if(rds_path == ""){
      rds_path = paste0("./BCA/BCA_plots/", speciesname, "_", SN, ".rds")
    }
    saveRDS(plotsp_plot, rds_path)
  }
  
  return(list(plot=plotsp_plot, summary_model=summary(model), statistics_anova=car::Anova(model)))

}


plotTukey = function(obj, diversity_data, plotx="funlab", ploty="value", facet.by="~ Season", mod = "pairwise ~ funlab | Season", padj = "tukey", labfill="", xlab="Fungicides", ylab="", title="", fontcolor="black", fontsize=12, labelsize=6){
  
  mod0 = as.formula(mod)
  tukeyout = emmeans::emmeans(obj, mod0, p.adjust.methods = padj)
  
  tukeydf = as.data.frame(multcomp::cld(tukeyout, Letters = letters))
  
  ym = min(diversity_data$value)
  
  g = ggplot(diversity_data, aes(x= !!sym(plotx), y = !!sym(ploty))) +
    geom_boxplot() +
    facet_grid(facet.by) +
    geom_text(data= tukeydf, aes(y=ym, label=.group), size=labelsize) +
    theme_bw() +
    labs(fill = labfill, x = xlab, y= ylab, title=title) +
    theme(text=element_text(size=12),
          axis.title.x = element_text(size=fontsize, margin=margin(t=10), color=fontcolor),
          axis.title.y = element_text(size=fontsize, margin=margin(r=10), color=fontcolor),
          axis.text.x =  element_text(size=fontsize, color=fontcolor),
          axis.text.y =  element_text(size=fontsize, color=fontcolor))

  
  return(g)
  
}


plot_lefse = function(run_lefse_out, phy=NULL, topn=10, colorSet = NULL,relevel=F, new_levels=NULL, padj_cutoff=0.05, xlab="LDA score(log10)", ylab="Genus", labfill = "Enriched group", fontcolor="black", fontsize=12, bartrans= 0.9, label_print=0){
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(ggtext)
  

  
  df = data.frame(marker_table(run_lefse_out))
  
  df$feature = scientific_name_formatter(df$feature)
  
  df$enrich_group = as.factor(df$enrich_group)

  if(relevel){
    if(is.null(new_levels)){
      stop("Levels need to be specific")
    } else {
      levels(df$enrich_group) = new_levels
    }
  }
  
  if(is.null(colorSet)){
    colorInput = colorRampPalette(brewer.pal(8, "Dark2"))(length(levels(df$enrich_group)))
  } else {
    colorInput = colorSet
  }
  
  df_order = df%>%
    arrange(.data$enrich_group, .data$ef_lda)
  
  feat = df_order$feature
  df_order$feature = factor(feat, levels = feat)
  
  df_top = df_order %>%
    group_by(enrich_group) %>%
    filter(padj < padj_cutoff) %>%
    slice_max(order_by = ef_lda, n = topn)
  
  if(label_print == 1){
    if(is.null(phy)){
      stop("Need phyloseq object")
    }
    tax_m = as.data.frame(tax_table(phy))
    taxa  = tax_m %>% 
      filter(Genus %in% as.character(df_top$feature)) %>%
      mutate(Phylum = scientific_name_formatter(Phylum), Genus = scientific_name_formatter(Genus)) %>%
      mutate(new_feature = paste(Phylum, Genus, sep = "-")) %>%
      select(new_feature, Genus) %>%
      distinct()
    
    rownames(taxa) = taxa$Genus
    
    df_top_merge = merge(df_top, taxa, by.x = "feature", by.y = "Genus")
    
    df_top_merge = df_top_merge %>%
      arrange(.data$enrich_group, .data$ef_lda)
    
    feat = df_top_merge$new_feature
    df_top_merge$new_feature = factor(feat, levels = feat)
    
    print(df_top_merge)
    
    g = df_top_merge %>%
      group_by(enrich_group) %>%
      ggplot(aes(x=ef_lda, y=new_feature, fill=enrich_group)) +
      geom_col(alpha=bartrans) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_discrete(expand = c(0.01, 0.01))+ 
      scale_fill_manual(values = colorInput) + 
      labs(fill=labfill, x=xlab, y=ylab) +
      theme(
        text=element_text(size=12),
        axis.title.x = element_text(size=fontsize, margin=margin(t=10), color=fontcolor),
        axis.title.y = element_text(size=fontsize, margin=margin(r=10), color=fontcolor),
        axis.text.x =  element_text(size=fontsize, color=fontcolor),
        axis.text.y =  element_markdown(size=fontsize, color=fontcolor),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(colour = "black", linetype = 8)) +
      guides(guide_legend())
    

     
  } else {

    
    g = df_top %>%
      ggplot(aes(x=ef_lda, y=feature, fill=enrich_group)) +
      geom_col(alpha=bartrans) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_discrete(expand = c(0.01, 0.01))+ 
      scale_fill_manual(values = colorInput) + 
      labs(fill=labfill, x=xlab, y=ylab) +
      theme(
        text=element_text(size=12),
        axis.title.x = element_text(size=fontsize, margin=margin(t=10), color=fontcolor),
        axis.title.y = element_text(size=fontsize, margin=margin(r=10), color=fontcolor),
        axis.text.x =  element_text(size=fontsize, color=fontcolor),
        axis.text.y =  element_markdown(size=fontsize, color=fontcolor),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(colour = "black", linetype = 8)) +
      guides(guide_legend())
  }


  


    
  return(list(graph=g, df=df))
  
}

extractCornCob = function(moddf, outputfile=T, outname=""){
  
  library(dplyr)
  
  startdf = data.frame(matrix(ncol = 7, nrow = 0))
  colnames(startdf) = c("taxa","mu.time_label1", "mu.time_label7", "mu.time_label14", "mu.funlabDL", "mu.funlabBM", "mu.funlabC2")
  
  spnames = names(moddf$mod)
  
  print(spnames)
  
  for(i in spnames){
    print(i)
    
    spdf = as.data.frame(moddf$mod[[i]]$coefficients)[, c("Estimate", "Pr(>|t|)")]
    
    colnames(spdf) = c("Est", "Pval")
    spdf2 = spdf %>%
      mutate(OUTtest = ifelse(Est > 0, "+", ifelse(Est < 0, "-", "0"))) %>%
      dplyr::select(OUTtest, Pval)
    
    spdf2$Pval = round(spdf2$Pval,4)
    tsp_df = t(spdf2)[, c("mu.time_label1", "mu.time_label7", "mu.time_label14", "mu.funlabDL", "mu.funlabBM", "mu.funlabC2")]
    
    tsp_df2 = as.data.frame(tsp_df)
    tsp_df2$taxa = i
    
    print(head(tsp_df2))
    tsp_df2 = tsp_df2[,c("taxa","mu.time_label1", "mu.time_label7", "mu.time_label14", "mu.funlabDL", "mu.funlabBM", "mu.funlabC2")]
    
    startdf = rbind(startdf, tsp_df2)
  }
  
  if(outputfile){
    write.csv(startdf, outname, quote = F, row.names = F)
  }
  
  return(startdf)
}

plotLinda = function(lindaplotobj, ylab = "Taxa", xlab = "Log2 Fold Changes", titles=NULL, pointsize=3, pointalpha=0.8, slice_max_n = 50, dotcolor="blue4", dotlabel="Corrected"){
  library(ggplot2)
  library(ggtext)
  
  figs = list()
  
  for(i in 1:length(lindaplotobj$plot.lfc)){
    
    if( (!is.null(titles)) & (length(titles) > 1)){
        Title = titles[i]
      } else {
        Title = paste0("fig",i)
      }
      
    print(Title)
    
    theData = lindaplotobj$plot.lfc[[i]]$data
    
    fig_n = theData %>%
      mutate(OTU = gsub("(^OTU[0-9]+?:).*", "\\1", Taxa)) %>%
      mutate(only_taxa = gsub("(^OTU[0-9]+?:)(.*)", "\\2", Taxa)) %>%
      mutate(italic_taxa = ifelse(grepl(" sp.", only_taxa),  
                              paste(paste0("*", gsub(" sp.$", "", only_taxa), "*"), " sp."), paste0("*", only_taxa, "*"))) %>%
      mutate(new_taxa = paste0(OTU, italic_taxa)) %>%
      dplyr::filter(bias == "Debiased") %>%
      mutate(sign=ifelse(Log2FoldChange > 0, "pos", "neg")) %>%
      group_by(sign) %>%
      slice_max(order_by = abs(Log2FoldChange), n=slice_max_n) %>%
      ggplot(aes(x = Log2FoldChange, y=reorder(new_taxa, Log2FoldChange))) +
      geom_errorbar(aes(xmin = Log2FoldChange - 1.96 * lfcSE,
                        xmax = Log2FoldChange + 1.96 * lfcSE), width = .2) +
      geom_point(aes(color = bias), size = pointsize, alpha=pointalpha) +
      geom_vline(xintercept = 0, color = 'gray', linetype = 'dashed') +
      geom_vline(xintercept = 2, color = 'gray', linetype = 'dashed') +
      geom_vline(xintercept = -2, color = 'gray', linetype = 'dashed') +
      scale_color_manual(values = dotcolor, labels=dotlabel) +
      labs(title = Title, y=ylab, x=xlab, color="Bias correction") +
      theme_classic() +
      theme(
            legend.position = "none",
            text = element_text(size=12),
            panel.grid.major.y = element_line(),
            axis.title.x = element_text(margin = margin(t=10), color="black", size=12),
            axis.title.y = element_text(margin = margin(r=10), color="black", size=12),
            axis.text.x = element_text(color="black", size=12),
            axis.text.y = element_markdown(color="black", size=12))
     
     figs[[Title]] = fig_n
  }
  
  return(figs)
}

summLinda = function(df, pv=0.05, log2=1, slicemax=5){
  
  df$taxa = gsub("OTU[0-9]+?:", "", rownames(df))
  
  df2 = df %>%
    filter(padj < pv & abs(log2FoldChange) >= log2) %>%
    mutate(signs = ifelse(log2FoldChange > 0, "+", "-")) %>%
    group_by(signs,taxa) %>%
    summarise(howmany = n()) %>%
    arrange(desc(howmany), .by_group = T) %>%
    slice_max(order_by = howmany, n=slicemax)
  
  return(df2)
}

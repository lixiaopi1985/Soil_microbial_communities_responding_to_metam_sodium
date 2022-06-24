

plotHeat_matrix = function(psmelt, groupby = c("period", "cul_type"), ntop=20, nshow=10){
  # make heatmap of this (top 10?)
  library(tidyr)
  
  top20 = psmelt %>%
    group_by(!!!syms(groupby)) %>%
    mutate(total_abund = sum(Abundance)) %>%
    mutate(new_genus =  Genus) %>%
    group_by(!!!syms(c(groupby, "new_genus"))) %>%
    mutate(rel_abund = 100*Abundance/total_abund) %>%# change cul_type to groupby
    summarise(accum_abund = sum(rel_abund)) %>%
    slice_max(order_by = accum_abund, n=ntop)
  #----------------------------------------------------------------------
  # top20 = psmelt %>%
  #   group_by(!!!syms(groupby)) %>%
  #   mutate(total_abund = sum(Abundance)) %>% # total in that grouping
  #   mutate(new_genus =  Genus) %>%
  #   group_by(!!!syms(c(groupby, "new_genus"))) %>% # change cul_type to groupby
  #   mutate(total_genus_in_group = sum(Abundance), rel_abund = 100*total_genus_in_group/total_abund) %>%
  #   summarise(accum_abund = mean(rel_abund)) %>%
  #   slice_max(order_by = accum_abund, n=ntop)
  # 
  # overall top abundant genera
  domdf = psmelt %>%
    mutate(total_abund = sum(Abundance), rel_abund = 100*Abundance/total_abund) %>%
    mutate(new_genus = Genus) %>%
    group_by(new_genus) %>%
    mutate(rel_abund = 100*Abundance/total_abund) %>%
    summarise(accum_abund = sum(rel_abund)) %>%
    arrange(desc(accum_abund))

  ##--------------------------------------------------------------
  # domdf = psmelt %>%
  #   mutate(total_abund = sum(Abundance)) %>%
  #   mutate(new_genus = Genus) %>%
  #   group_by(new_genus) %>%
  #   mutate(total_genus_in_group = sum(Abundance), rel_abund = 100*total_genus_in_group/total_abund) %>%
  #   summarise(accum_abund = mean(rel_abund)) %>%
  #   arrange(desc(accum_abund))
  
  doms = domdf$new_genus[1:nshow]
  
  # some top genus might be not in other group, so find them as well
  top_accum = psmelt %>%
      group_by(!!!syms(groupby)) %>%
      mutate(total_abund = sum(Abundance), rel_abund = 100*Abundance/total_abund) %>%
      mutate(new_genus = Genus) %>%
      group_by(!!!syms(c(groupby, "new_genus"))) %>%
      mutate(rel_abund = 100*Abundance/total_abund) %>%
      summarise(accum_abund = sum(rel_abund)) %>%
      filter(new_genus %in% unique(top20$new_genus)) %>%
      arrange(desc(accum_abund))
  
  #----------------------------------------------------------------------
  # top_accum = psmelt %>%
  #   group_by(!!!syms(groupby)) %>%
  #   mutate(total_abund = sum(Abundance)) %>%
  #   mutate(new_genus = Genus) %>%
  #   group_by(!!!syms(c(groupby, "new_genus"))) %>%
  #   mutate(total_genus_in_group = sum(Abundance), rel_abund = 100*total_genus_in_group/total_abund) %>%
  #   summarise(accum_abund = mean(rel_abund)) %>%
  #   filter(new_genus %in% unique(top20$new_genus)) %>%
  #   arrange(desc(accum_abund))

  
  # hard coded
  # top.meta = top_accum %>%
  #   unite(col = samples, period, metamStatus, remove = F)
  # 
  # top.meta.2 = top_accum %>%
  #   unite(col = Samples,period, metamStatus) %>%
  #   spread(Samples, accum_abund) %>%
  #   as.data.frame()
  top.meta.2 = top_accum %>%
    unite(col=Samples, !!!syms(groupby)) %>%
    spread(Samples, accum_abund) %>%
    as.data.frame()
  
  rownames(top.meta.2) = top.meta.2$new_genus
  
  # matrix
  top.matrix = as.matrix(top.meta.2[, -1])
  
  return(list(m = top.matrix, dom = doms, overal=domdf))
}

#---------------------------------------------------------------------------------------------------

plotHeat_plot = function(tops, COLOR=NA, new_break = NA,fontsize=12, out=T, outfile="heatmap.png", Legend=T, anno_legend=T, cellw = 40, cellh = 40, ...){
  
  library(stringr)
  library(RColorBrewer)
  library(pheatmap)
  
  top.matrix = tops$m
  doms = tops$dom
  
  org_colname = colnames(top.matrix)
  
  annoCol = data.frame(Fumigation_history= gsub("[0-9]wk_", "", org_colname),
                       Time= str_extract(org_colname, "[0-9]wk"))
  
  print(annoCol)
  pasteName = function(x){
    return(paste0("Microcosm T", x))
  }
  
  

  # if(is.na(COLOR)){
  #   color_set_cat = brewer.pal(7, "Accent")
  #   
  # } else {
  #   color_set_cat = COLOR
  # }
  
  color_set_time = brewer.pal(9, "BuGn")
  
  # set category manually
  # colorPat = colorRampPalette(color_set_cat)
  colorTime = colorRampPalette(color_set_time)
  
  
  rownames(annoCol) = org_colname
  annoCol
  annoCol$Fumigation_history= factor(annoCol$Fumigation_history, levels = c("Not fumigated", "Fumigated"))
  levels(annoCol$Fumigation_history) = c("Not fumigated", "Fumigated")
  annoCol$Time = factor(annoCol$Time, levels = c("0wk", "1wk", "3wk", "6wk"))
  
  top.matrix.ord = top.matrix[,order(annoCol$Fumigation_history)]
  top.matrix.ord2 = top.matrix.ord[doms,]
  
  # crop_color = colorPat(length(unique(annoCol$Fumigation_history)))
  crop_color = c("red4", "blue4")
  time_color = colorTime(length(unique(annoCol$Time)))
  names(crop_color) = unique(annoCol$Fumigation_history)
  names(time_color) = unique(annoCol$Time)
  
  
  anno_colors = list(Time=time_color, Fumigation_history=crop_color)
  

    # without italicize the label
  if(out){
      g = pheatmap::pheatmap(top.matrix.ord2, 
                             annotation_col = annoCol,
                             cellwidth = cellw, 
                             cellheight = cellh, 
                             color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(length(new_break)), 
                             annotation_colors = anno_colors, 
                             cluster_cols = F, 
                             cluster_rows = F, 
                             breaks = new_break, 
                             fontsize = fontsize, 
                             filename = outfile, 
                             legend = Legend,
                             annotation_legend = anno_legend, 
                             # column_split = annoCol$Fumigation_history, 
                             ...)
      
      
  } else {
      g = pheatmap::pheatmap(top.matrix.ord2, 
                             annotation_col = annoCol,
                             cellwidth = cellw, 
                             cellheight = cellh,
                             # color = colorRampPalette(c("navy", "white", "firebrick3"))(length(new_break)),
                             color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(length(new_break)), 
                             annotation_colors = anno_colors, 
                             cluster_cols = F, 
                             cluster_rows = F, 
                             breaks = new_break, 
                             fontsize = fontsize, 
                             # filename = outfile, 
                             legend = Legend,
                             annotation_legend = anno_legend, 
                             # column_split = annoCol$Fumigation_history, 
                             ...)
    }


  
  return(g)
}


# add significant * to the row labels
addSig = function(origlabel, cond1, cond2){
  
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










aggAbund_barplot_ver2 = function(phylo_melt, taxlevel, colors = NA, colorSet="Paired", ncolor=12, facet_ = "cul_type", ggplotx="field", groupby="field", xLab=NULL, yLab="Relative abundance (%)", perc = 1, order_tax_by=NA, legpos = "bottom", fontsize=12, guide_col=1, leg_fill="", percLabeltop = F, OrderByAbund = T, greyColor = "grey90", ...){
  
  
  require(ggplot2)
  require(pals)
  require(tidyverse)
  require(rlang)
  

  
  perclabel = paste0("<", perc, "%")

  if(length(facet_) > 1){
    facet_arg = paste0(facet_[1], "~", paste(facet_[2:length(facet_)], sep = "+", collapse = "+"))
  } else if(length(facet_) == 1){
    facet_arg = paste0("~", facet_[1])
  }
  
  print(facet_arg)
  
  # old parts of code for reference
  # df = phylo_melt %>%
  #   dplyr::group_by(!!!syms(groupby)) %>%
  #   mutate(total_abund = sum(Abundance), rel_abund = 100*Abundance/total_abund) %>%
  #   mutate(new_genus = ifelse(rel_abund < perc, perclabel, !!sym(taxlevel)))
  
  df = phylo_melt %>%
    dplyr::group_by(!!!syms(groupby)) %>%
    mutate(total_group = sum(Abundance)) %>% # calculate total abundance in the groups but the genus abundance is still at OTU level we need to get that too
    dplyr::group_by(!!!syms(c(groupby, taxlevel))) %>%
    mutate(total_group_tax = sum(Abundance), rel_abund = 100*total_group_tax / total_group) %>% # here is the total abundance of a taxon in that grouping
    mutate(new_genus = ifelse(rel_abund < perc, perclabel, !!sym(taxlevel))) 
  
  # sort by most abundant (sum up) in that group
  topAbund =  df %>%
    ungroup() %>%
    group_by(new_genus) %>%
    summarize(acc_abund = mean(total_group_tax)) %>% # summ abundance by new genus
    arrange(desc(acc_abund)) %>%
    pull(new_genus)
  
  if(!is.na(order_tax_by)){
    orderTaxa = topAbund[topAbund == order_tax_by]
  } else {
    orderTaxa = topAbund[1]
  }
  # topAbund2 = c(perclabel, topAbund[-which(topAbund == perclabel)])
  # print(topAbund2)
  # by total counts
  # topAbund =  phylo_melt %>%
  #   group_by(!!sym(taxlevel)) %>%
  #   summarize(total_abund = sum(Abundance)) %>%
  #   arrange(desc(total_abund)) %>%
  #   pull(!!sym(taxlevel))
  # 
  # print(topAbund)
  
  
  # topTaxa = topAbund[[taxlevel]][1]
  
  

  # sum order by a taxa, default the most abundant
  dforder = df %>%
      group_by(!!!syms(groupby), new_genus) %>%
      summarize(acc_abund = mean(rel_abund)) %>%
      filter(new_genus == orderTaxa) %>%
      arrange(desc(acc_abund)) %>%
      pull(!!sym(ggplotx))
 
  
  n_tax = length(topAbund)
  
  if(is.na(colors)){
    colPatte = colorRampPalette(brewer.pal(ncolor, colorSet))
  } else {
    colPatte = colorRampPalette(colors)
    
  }
  
  colorValues = colPatte(n_tax)
  names(colorValues) = unique(topAbund)
  
  colorValues[perclabel] = greyColor
  
  if(percLabeltop){
    removeNames = names(colorValues)[-which(names(colorValues) == perclabel)]
    colorValues = colorValues[c(perclabel, removeNames)]
  }


  if(OrderByAbund){
    outplot = df %>%
      # group_by(!!!syms(c(groupby, "new_genus"))) %>%
      # mutate(accu_counts = sum(total_group_tax)) %>%
      mutate(order_g = factor(new_genus, levels = topAbund)) %>%
      mutate(order_f = factor(!!sym(ggplotx), levels = dforder)) %>%
      ggplot(aes(x = order_f, y = total_group_tax, fill = order_g)) + # rel_abund == Abund
      geom_bar(stat = "identity", position = "fill") +
      facet_grid(facet_arg, scales = "free_x", space="free_x", ...) +
      scale_fill_manual(values = colorValues) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0,0)) +
      labs(x=xLab, y=yLab, fill=leg_fill) +
      theme_bw()+
      theme(axis.text.x = element_text(size=fontsize, color="black"),
            axis.text.y = element_text(size=fontsize, color="black"),
            legend.text = element_text(size=fontsize, color="black"),
            axis.title.x = element_text(margin = margin(t=5), color="black"),
            axis.title.y = element_text(margin = margin(r=5), color="black"),
            legend.position = legpos, 
            text = element_text(size=fontsize, color="black")) +
      guides(fill=guide_legend(ncol=guide_col))
  } else {
    outplot = df %>%
      mutate(order_g = factor(new_genus, levels = topAbund)) %>%
      # mutate(order_f = factor(!!sym(ggplotx), levels = unique(dforder))) %>%
      ggplot(aes(x = !!sym(ggplotx), y = total_group_tax, fill = order_g)) +
      geom_bar(stat = "identity", position = "fill") +
      facet_grid(facet_arg, scales = "free_x", space="free_x", ...) +
      scale_fill_manual(values = colorValues) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0,0)) +
      labs(x=xLab, y=yLab, fill=leg_fill) +
      theme_bw()+
      theme(axis.text.x = element_text(size=fontsize, color="black"),
            axis.text.y = element_text(size=fontsize, color="black"),
            legend.text = element_text(size=fontsize, color="black"),
            axis.title.x = element_text(margin = margin(t=5), color="black"),
            axis.title.y = element_text(margin = margin(r=5), color="black"),
            legend.position = legpos, 
            text = element_text(size=fontsize, color="black")) +
      guides(fill=guide_legend(ncol=guide_col))
  }

  
  return(graph=outplot)
  
}
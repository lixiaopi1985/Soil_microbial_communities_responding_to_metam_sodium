plotCoef = function(result_diffTest, phylo, fdr = 0.05, taxLevel="Genus",orderList = NA){
  
  all_fdr_values = result_diffTest$p_fdr[!is.na(result_diffTest$p_fdr)]
  try_species = names(all_fdr_values)
  taxaNames = otu_to_taxonomy(try_species, data = phylo, level = taxLevel) # this order is different from dominant 
  all_models = result_diffTest[["all_models"]][!is.na(result_diffTest[["all_models"]])]
  
  names(all_models) = taxaNames
  
  if(!is.na(orderList)){
    
    if( (length(orderList) == length(all_models)) | (length(orderList) == length(taxaNames))){
      all_models = all_models[orderList]
      all_fdr_values = all_fdr_values[order(match(taxaNames, orderList))]
      taxaNames = taxaNames[order(match(taxaNames, orderList))]
      
    } else {
      stop("orderList length is different to the model or taxaNames")
    }
    
  }
  coef = c()
  maxse = c()
  minse = c()
  
  vari = c()
  vari_maxse = c()
  vari_minse = c()
  
  
  for(i in 1:length(all_models)){
    
    print(i)
    print(paste("Order list:", orderList[i]))
    print(paste("Model name:", names(all_models[i])))
    
    coef = c(coef, all_models[[i]]$coefficients[2,1]) 
    maxse = c(maxse, all_models[[i]]$coefficients[2,1]+qnorm(0.975)*all_models[[i]]$coefficients[2,2])
    minse = c(minse, all_models[[i]]$coefficients[2,1]-qnorm(0.975)*all_models[[i]]$coefficients[2,2])
    
    vari = c(vari, all_models[[i]]$coefficients[4,1])
    vari_maxse = c(vari_maxse, all_models[[i]]$coefficients[4,1]+qnorm(0.975)*all_models[[i]]$coefficients[4,2])
    vari_minse = c(vari_minse, all_models[[i]]$coefficients[4,1]-qnorm(0.975)*all_models[[i]]$coefficients[4,2])
    
  }
  
  
  df = data.frame(coef = coef, variability = vari, xmin = minse, xmax = maxse, var_xmin = vari_minse, var_xmax = vari_maxse, taxa = taxaNames, fdr = all_fdr_values)
  
  return(list(DF=df, mod=all_models))
}

plotCoef_mcm = function(result_diffTest, phylo, fdr = 0.05, taxLevel="Genus", nlevel = 3,orderList = NA){
  
  # specific for mcm
  all_fdr_values = result_diffTest$p_fdr[!is.na(result_diffTest$p_fdr)]
  sigTaxa = otu_to_taxonomy(result_diffTest$significant_taxa, data=phylo, level=taxLevel)
  
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
    var_levels = c(var_levels, gsub("mu.", "", names(all_models[[i]]$coefficients[2:4,1])))
    
    # coefficients
    coef = c(coef, all_models[[i]]$coefficients[2:4,1]) 
    maxse = c(maxse, all_models[[i]]$coefficients[2:4,1]+qnorm(0.975)*all_models[[i]]$coefficients[2:4,2])
    minse = c(minse, all_models[[i]]$coefficients[2:4,1]-qnorm(0.975)*all_models[[i]]$coefficients[2:4,2])
    # coefficients pvalu
    
    # variability
    vari = c(vari, all_models[[i]]$coefficients[6:8,1])
    vari_maxse = c(vari_maxse, all_models[[i]]$coefficients[6:8,1]+qnorm(0.975)*all_models[[i]]$coefficients[6:8,2])
    vari_minse = c(vari_minse, all_models[[i]]$coefficients[6:8,1]-qnorm(0.975)*all_models[[i]]$coefficients[6:8,2])
    
  }

  
  df = data.frame(coef = coef, facet_var = var_levels, variability = vari, xmin = minse, xmax = maxse, var_xmin = vari_minse, var_xmax = vari_maxse, taxa = taxaNames_rep, siglevels=sigTaxa_w_symbol, fdr = fdr_rep)
  
  return(list(DF=df, mod=all_models))
}

# add significant * to the row labels
addSig_1cond = function(origlabel, cond1){
  
  new_labs = c()
  for(i in origlabel){
    if( i %in% cond1 ){
      sigLab = paste0("*", i)
      new_labs = c(new_labs, sigLab)
    } else {
      sigLab = i
      new_labs = c(new_labs, sigLab)
    }
  }
  
  return(new_labs)
}

# plotDE = function(mergeddf, x="coef", y="genus_order", xLab = "Coefficients", yLab = "", plottitle = "", facetLabeller = NA, vline_color = "gray50", vline_lty = "dashed"){
#   
#   theplot = mergeddf %>%
#      ggplot(aes(y = reorder(!!sym(y, desc(y)), x = x)) +
#      geom_vline(xintercept = 0, color = vline_color, lty=vline_lty, alpha=0.75, lwd = 1) +
#      geom_point() +
#      geom_errorbarh(aes(xmin = xmin, xmax=xmax), height=0.2) +
#      theme_minimal() +
#      labs(title = plottitle, x = xLab, y = yLab) +
#      facet_grid(~ facet_var, labeller = facetLabeller) +
#      # scale_y_discrete(expand = c(0,0)) +
#      scale_x_continuous(breaks = scales::pretty_breaks(n=6)) +
#      ggplot2::theme(axis.title.x = element_text(size=12,color="black", margin = margin(t=10)),
#                     axis.title.y = element_text(size=12,color="black", margin = margin(r=10)),
#                     axis.text.x = element_text(color="black", size=12),
#                     axis.text.y = element_text(size=12,color="black", margin = margin(r=10)),
#                     axis.line.x = element_line(),
#                     axis.ticks.y = element_blank(),
#                     strip.text.x = element_text(size=12, face = "bold", margin = margin(b=10)),
#                     panel.grid.major.y = element_blank())+
#      plot_layout(tag_level = "new"))
#   
# }
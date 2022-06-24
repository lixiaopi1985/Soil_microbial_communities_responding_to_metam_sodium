plotReg = function(df, x, y, xLabs = "", yLabs = "", linecolor = "red", regFont=5, fontsize=12, minLimit = 0, maxLimit = 100, gap = 100, labelpos.y=0){
  library(ggpubr)
  
  if(labelpos.y==0){
    labelpos = maxLimit
  } else {
    labelpos = labelpos.y
  }
  
  outgraph = df %>%
  ggplot(aes(x = !!sym(x), y = !!sym(y))) +
    geom_point() +
    geom_smooth(method = "lm", formula = y~x, color=linecolor) +
    scale_y_continuous(limits = c(minLimit,maxLimit)) +
    stat_regline_equation(aes(label = ..eq.label..), label.y = labelpos - gap*0.5, size=regFont) +
    stat_regline_equation(aes(label = ..rr.label..), label.y = labelpos - gap, size=5) +
    labs(x = xLabs, y=yLabs) +
    theme_classic() +
    theme(text = element_text(size=fontsize),
          axis.title.x = element_text(size=fontsize, color="black", margin = margin(t=5)),
          axis.title.y = element_text(size=fontsize, color="black", margin = margin(r=5)),
          axis.text.x = element_text(size=fontsize, color="black"),
          axis.text.y=element_text(size=fontsize, color="black"))
  
  return(outgraph)

}

plotBox = function(df, x, y, xLabs = "", yLabs = "", linecolor = "red", regFont=5, fontsize=12, minLimit = 0, maxLimit = 100, gap = 100, moregap = 100, formulaLabs = "", label.pos.x = 0.8, extraAnno = F, Anno2 = NA){
  
  if(extraAnno){
    outgraph = df %>%
      ggplot(aes(x = !!sym(x), y = !!sym(y))) +
      geom_boxplot( alpha = 0.7, fill = "gray80") +
      annotate("text", x = Anno2[["Group"]], y = Anno2[[y]]+moregap, label=Anno2$Letter) +
      annotate("text", x = label.pos.x, y=maxLimit-gap*0.5, label = formulaLabs, parse=T) +
      # scale_x_discrete(labels=c("Conventional", "Organic")) +
      scale_y_continuous(limits = c(minLimit, maxLimit+gap*0.5)) +
      theme_bw() +
      labs(fill = "", x = xLabs, y = yLabs) +
      theme(
        text =  element_text(size=fontsize, color="black"),
        legend.text = element_text(size=fontsize),
        axis.title.x = element_text(margin = margin(t=5)),
        axis.title.y = element_text(margin = margin(r=5)),
        axis.text.x = element_text(color="black", size = fontsize),
        axis.text.y = element_text(color="black", size = fontsize))
  } else {
    
    outgraph = df %>%
      ggplot(aes(x = !!sym(x), y = !!sym(y))) +
      geom_boxplot( alpha = 0.7, fill = "gray80") +
      # annotate("text", x = letters16S.merge2$Group, y = letters16S.merge2$index_measure+gap, label=letters16S.merge2$Letter) +
      annotate("text", x = label.pos.x, y=maxLimit-gap*0.5, label = formulaLabs, parse=T) +
      # scale_x_discrete(labels=c("Conventional", "Organic")) +
      scale_y_continuous(limits = c(minLimit, maxLimit+gap*0.5)) +
      theme_bw() +
      labs(fill = "", x = xLabs, y = yLabs) +
      theme(
        text =  element_text(size=fontsize, color="black"),
        legend.text = element_text(size=fontsize),
        axis.title.x = element_text(margin = margin(t=5)),
        axis.title.y = element_text(margin = margin(r=5)),
        axis.text.x = element_text(color="black", size = fontsize),
        axis.text.y = element_text(color="black", size = fontsize))
  }

  
  
  return(outgraph)
  
}

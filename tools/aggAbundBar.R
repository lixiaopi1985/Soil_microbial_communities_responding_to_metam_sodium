aggAbund_barplot = function(phylo, taxlevel, colors, phylodf=F, x = "Sample", y = "Abundance", facet_ = "cul_type", summarizeby=TRUE, summarizeVariable = c("Genus", "field"), summarizemethod=mean(), ggplotx="field", xLab=NULL, yLab="Relative abundance (%)", perc = T, abund = 1/100, NArm = T, legpos = "bottom", fontsize=12, guide_col=1, leg_fill="", ...){
        
        
        require(ggplot2)
        require(pals)
        require(tidyverse)
        require(rlang)
        
        colPatte = colorRampPalette(colors)
        
        if(phylodf){
                dat =phylo
                dat[[taxlevel]] = as.character(dat[[taxlevel]])
        } else {
                if(perc){
                # agglomerate taxa
                glom = tax_glom(phylo, taxrank = taxlevel, NArm = NArm) %>%
                        transform_sample_counts(function(x)100*x/sum(x))
                abund = 100*abund
        } else {
                glom = tax_glom(phylo, taxrank = taxlevel, NArm = NArm) %>%
                        transform_sample_counts(function(x)x/sum(x))
        }
        

        dat = psmelt(glom) 
        # group by tax level than calculate mean
        dat[[taxlevel]] = as.character(dat[[taxlevel]])
        }
        
        # # filter
        # if(perc){
        #         dat[[taxlevel]] = ifelse(dat$Abundance < abund, paste("<", abund, "%"),  dat[[taxlevel]])
        # } else {
        #         dat[[taxlevel]] = ifelse(dat$Abundance < abund, paste("<", abund),  dat[[taxlevel]])
        # }
        # 
        # }
        
        
        facet_arg = paste0("~", facet_)
        
        if(perc){
                dat1 = dat %>%
                        dplyr::group_by(!!!syms(summarizeVariable)) %>%
                        dplyr::summarise(meanRA = mean(Abundance)) %>%
                        dplyr::mutate(new_taxlevel = ifelse(meanRA < 100*abund, paste0("<", 100*abund, "%"), !!sym(taxlevel)))
        } else {
                dat1 = dat %>%
                        dplyr::group_by(!!!syms(summarizeVariable)) %>%
                        dplyr::summarise(meanRA = mean(Abundance)) %>%
                        dplyr::mutate(new_taxlevel = ifelse(meanRA < abund, paste0("<", abund, "%"), !!sym(taxlevel)))
        }
        
        
        n_tax = length(unique(dat1[["new_taxlevel"]]))
        
        if(summarizeby){
                
                outplot = dat1 %>% 
                ggplot(aes_string(x=ggplotx, y="meanRA", fill="new_taxlevel")) +
                geom_bar(stat = "identity", position = "fill") +
                facet_grid(facet_arg, scales = "free_x", space="free_x", ...) +
                scale_fill_manual(values = colPatte(n_tax), na.value = "grey45") +
                scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
                labs(x=xLab, y=yLab, fill=leg_fill) +
                theme_bw()+
                theme(axis.text.x = element_text(angle = 90, hjust=1), legend.position = legpos, text = element_text(size=fontsize)) +
                guides(fill=guide_legend(ncol=guide_col))
        
        }else{
                outplot = 
                        dat1 %>%
                        ggplot(aes_string(x=x, y=y, fill="new_taxlevel")) +
                        geom_bar(stat = "identity", position = "fill") +
                        facet_grid(facet_arg, scales = "free_x", space="free_x", ...) +
                        scale_fill_manual(values = colPatte(n_tax), na.value = "grey45") +
                        scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
                        labs(x=xLab, y=yLab,fill=leg_fill) +
                        theme_bw()+
                        theme(axis.text.x = element_text(angle = 90, hjust=1), legend.position = legpos, text = element_text(size=fontsize)) +
                        guides(fill=guide_legend(ncol=guide_col))
        }
        

        return(list(data=dat, plot=outplot))

}

# --------------------------------------------------------------------------------------------#
# Figure 6 and S7: Results of LDSC Partitioned Heritability
# bar plots of enrichment score and radar plot of enrichment p-value
# requires LDSC data files
# --------------------------------------------------------------------------------------------#
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
# --------------------------------------------------------------------------------------------#
# each file is a GWAS trait (15)
setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")

# /proj/yunligrp/users/lagler/SIP/LDSC/resultsPIRspecific/LDSC_resultsPIR_CTspecific.txt
specific <- fread("LDSC_resultsPIR_CTspecific.txt") 

# /proj/yunligrp/users/lagler/SIP/LDSC/resultsPIR/LDSC_resultsPIR.txt
allsips <- fread("LDSC_resultsPIR.txt")

# --------------------------------------------------------------------------------------------#
# table of trait abbreviations
traits <- data.table(Abbreviation=c("RBC", "HGB", "HCT", "MCV", "MCH", "MCHC", "RDW",
                                    "WBC", "NEU", "MONO", "LYM", "BASO", "EOS",
                                    "PLT", "MPV"),
                     `GWAS Trait`=c("Red Blood Cell Count", "Hemoglobin", "Hematocrit", "Mean Corpuscular Volume",
                                    "Mean Corpuscular Hemoglobin  ", "MCH Concentration", "RBC Distribution Width",
                                    "White Blood Cell Count", "Neutrophil Count", "Monocyte Count", "Lymphocyte Count",
                                    "Basophil Count", "Eosinophil Count",
                                    "Platelet Count", "Mean PLT Volume"))
# --------------------------------------------------------------------------------------------#

CTcolor <- c("#00a450","#00b2de","#25519c","#ec9136","#FF66CC")

# bar plot of enrichment scores
enrichScore <- function(i, data){
  alpha <- ifelse(data[Category==CT[i],]$pstar == "*", 1, .7)
  ggplot(data[Category==CT[i],], aes(x=trait2, y=Enrichment, fill=Category))+
    geom_col(position=position_dodge(width=.8), width=.8, alpha=alpha) +
    geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error,
                      ymax=Enrichment+Enrichment_std_error), alpha=alpha,
                  width=.2, position=position_dodge(.8)) + 
    geom_text(aes(label=pstar, y=(Enrichment+Enrichment_std_error)*1.03),
              position=position_dodge(width=.8),
              color="red", show.legend=F)+
    ggtitle(cellnames[i]) + 
    scale_fill_manual(values = CTcolor[i]) +
    coord_flip() + theme_bw() + 
    theme(legend.position = "none", axis.title = element_blank()) 
}

# radar plot theme
plotTheme <- theme_minimal() +
  theme(legend.position = "top", # change to none if making combined figure
        axis.text.x = element_text(size=12, face="bold"),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        plot.title=element_text(size=16),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = unit(rep(1,4), "mm"), 
        legend.text=element_text(size=10),# comment out if making combined figure
        legend.title=element_blank() # comment out if making combined figure
  )
# --------------------------------------------------------------------------------------------#

#mkplots <- function(data){
  # reorder for plotting (red blood cell, white blood cell, platelet)
  data$trait <- factor(data$trait, 
                       levels=c("HCT", "HGB", "MCH", "MCHC", "MCV", "RBC", "RDW",
                                "BASO", "EOS", "LYM", "MONO", "NEU", "WBC",
                                "PLT", "MPV"))
  # reverse trait order for flipped coordinates
  data$trait2 <- factor(data$trait, levels=rev(levels(data$trait)))
  
  # rename cell types
  data$Category <- factor(data$Category, labels = c("Ery", "MacMon", "MK", "nCD4","Neu"))
  # store cell types
  CT <- levels(data$Category)
  # store unabbreviated names
  cellnames <- c(Ery="Erythrocyte", MacMon="Macrophage/Monocyte", MK="Megakaryocyte",
                 nCD4="Naive CD4 T-cell", Neu="Neutrophil")
  
  # define significance
  setDT(data)[, pstar := ifelse(Enrichment_p < 0.05, "*", "")]
  # -log10 transform
  setDT(data)[, log10p := -log10(Enrichment_p)]
  
  print(summary(data$Enrichment)) # .3 - 8.3
  print(summary(data$log10p)) # 0.007 - 3.75 (1.3 is significant)
  
  # re-order to match radar plot
  traits$Abbreviation <- factor(traits$Abbreviation, levels=levels(data$trait))
  setorder(traits)
  
  # make table
  traitTable <- ggtexttable(traits, rows = NULL,
                            theme=ttheme(base_size=6, 
                                         tbody.style = tbody_style(hjust=0, x=.1, fill="white"),
                                         colnames.style = colnames_style(fill = "white")))
  
  # bar plots of enrichment scores
  scorePlots <- ggarrange(enrichScore(1, data), enrichScore(2, data), enrichScore(3, data),
                          enrichScore(4, data), enrichScore(5, data),
                          traitTable, align="hv") %>%
    annotate_figure(.,
                    left = text_grob("GWAS Trait", color = "black", rot = 90, size=16),
                    bottom = text_grob("Enrichment Score", color = "black", size=16))
  
  # radar plot of enrichment p-values
  pvalPlot <- ggplot(data, aes(x=trait, y=log10p, fill=Category))+
    geom_hline(yintercept = c(1.3, 2.3, 3.3), color="grey")+
    geom_col(position=position_dodge(width=.8), width=.8) +
    geom_hline(yintercept = -.1, color="grey40") +
    # geom_text(aes(label=pstar, y=log10p+.05), position=position_dodge(width=.8),
    #           color="red", show.legend=F)+
    annotate("text", x=rep(0.5, 3), y=c(1.3, 2.3, 3.3), label=c(1.3, 2.3, 3.3), size=4)+
    scale_fill_manual(values = CTcolor) + 
    coord_polar() + ylim(-.75,3.8)+
    ggtitle("Enrichment p-value (-log10)")+ plotTheme
  
#   return(list=c(pval=pvalPlot, score=scorePlots))
# }


# save each separately
# cell type specific SIPs
scorePlots
ggsave("../figures/FigureS7_LDSCbar_specific.png", height=10, width=16)
ggsave("../figures/FigureS7_LDSCbar_specific.pdf", height=10, width=16)

pvalPlot+ggtitle("")
ggsave("../figures/Figure6_LDSCpval_specific.png", height=8, width=8)
ggsave("../figures/Figure6_LDSCpval_specific.pdf", height=8, width=8)

#pval plot for PowerPoint
pvalPlot + ggtitle("") + 
  theme(panel.background = element_rect(fill = "transparent", color=NA), 
  plot.background = element_rect(fill = "transparent", color = NA),
  legend.position = "none")
ggsave("../figures/Figure6_LDSCpval_specific_PowerPoint.png",
       height=8, width=8, bg="transparent")

# all SIPs
scorePlots
ggsave("../figures/FigureS7_LDSCbar_allSIPs.png", height=10, width=16)
ggsave("../figures/FigureS7_LDSCbar_allSIPs.pdf", height=10, width=16)

# pvalPlot+ggtitle("")
# ggsave("../figures/Figure6_LDSCpval_allSIPs.png", height=8, width=8)
# ggsave("../figures/Figure6_LDSCpval_allSIPs.pdf", height=8, width=8)



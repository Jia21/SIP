#----------------------------------------------------------------------------------------------#
# Figure 2(c-d): bar and box plots of number of super-PIRs and PIR score
# requires fraser.rda and SIP_*.txt (from defineSIPlot.R)
#----------------------------------------------------------------------------------------------#
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")
#setwd("/proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/") #fraser.rda
#----------------------------------------------------------------------------------------------#
# Define super PIR and PIR score - get table of results
#----------------------------------------------------------------------------------------------#
# formatted for intra-chromosomal interactions and distance < 2MB
load("fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[31] <- "MacMon"

getData <- function (cell, nInt){
  # select significant interactions
  hic <- fraser %>% filter(!!as.name(cell) >= 5) %>% data.table
  
  # count number of interactions per PIR 
  nbaits <- hic[, .N, by="oeID"]
  # nrow(nbaits[N==1])/nrow(nbaits) ~59% interact with a single promoter fragment
  # nrow(nbaits[N>=4])/nrow(nbaits) ~10% interact with 4 or more promoter fragments
  
  # merge with hic data
  hic2 <- merge(hic, nbaits, by="oeID")
  # find max number of interactions per PIR (PIR score)
  maxN <- hic2[, .("max"=max(N)), by="baitID"]
  # merge back with hic data
  hic3 <- unique(merge(maxN, hic2,
                      by="baitID", all.x=T, all.y=F))
  
  # add SIP indicator (changes data from interaction level to bait level)
  # /proj/yunligrp/users/lagler/SIP/data/
  sip <- fread(paste0("SIP_", cell, ".txt"))
  pir <- merge(sip, hic3[, c("baitID", "max")],
               by="baitID", all.x=T, all.y=F) %>% unique(., by="baitID")
  
  # define super PIR as >=4
  setDT(pir)[, superPIR := ifelse(max>=nInt, 1, 0)]
  
  yy <- nrow(pir[SIP==1 & superPIR==1]) # super PIR and SIP
  yn <- nrow(pir[SIP==1 & superPIR==0]) # regular PIR and SIP
  ny <- nrow(pir[SIP==0 & superPIR==1]) # super PIR and non-SIP
  nn <- nrow(pir[SIP==0 & superPIR==0]) # regular PIR and non-SIP
  
  # 2x2 table of counts plus ratio
  counts <- data.table("CT" = rep(cell, 2),
                       "type" = c("SIP", "non-SIP"),
                       "superPIR" = c(yy, ny),
                       "PIR" = c(yn, nn))
  setDT(counts)[, Ratio := superPIR/(superPIR+PIR)]
  
  # chisq test
  p <- chisq.test(counts[,3:4])$p.value
  setDT(counts)[, chsq.pval := p]
  
  # median PIR score per group
  medPIR <- pir[, .("med" = median(max)), by="SIP"]
  counts$medPIR <- medPIR[c(2,1), med]
  
  # wilcoxon test for PIR score
  w <- wilcox.test(pir[SIP==1, max], pir[SIP==0, max])$p.value
  setDT(counts)[, wilcox.pval := w]
  
  # add cell type
  setDT(pir)[, CT := cell]
  return(list("table"=counts, "data"=pir))

}

# call function for each cell type
pir <- rbind(getData("Ery", 4)$table, getData("MacMon", 4)$table,
             getData("MK", 4)$table, getData("nCD4", 4)$table, getData("Neu", 4)$table)
scores <- rbind(getData("Ery", 4)$data, getData("MacMon", 4)$data,
                getData("MK", 4)$data, getData("nCD4", 4)$data, getData("Neu", 4)$data)

# check other cutoff values for defining super-PIRs
nInt = 10
rbind(getData("Ery", nInt)$table, getData("MacMon", nInt)$table,
      getData("MK", nInt)$table, getData("nCD4", nInt)$table, getData("Neu", nInt)$table)

#----------------------------------------------------------------------------------------------#
# kable of results for supplementary table
#----------------------------------------------------------------------------------------------#

# keep formatting separate incase want to change digits
pirK <- pir
pirK$chsq.pval <- format(pirK$chsq.pval, digits=2) # format p-values
pirK$wilcox.pval <- format(pirK$wilcox.pval, digits=2)
pirK$chsq.pval[c(2,4,6,8,10)] <- NA # remove every other value
pirK$wilcox.pval[c(2,4,6,8,10)] <- NA
pirK$Ratio <- format(pirK$Ratio, digits=2, drop0trailing = F)
pirK$medPIR <- format(pirK$medPIR, digits=2, drop0trailing = F)
pirK[,1] <- NA

options(knitr.kable.NA = "")

kable(pirK, col.names = c("", "Bait Type", "# Super PIR",
                          "# Typical PIR", "Ratio", "Chi-sq p-value",
                          "Median PIR Score", "Wilcoxon p-value"),
      booktabs=T, "latex", align="llccccccc", 
      format.args = list(big.mark=",")) %>%
  pack_rows(index = c("Ery" = 2, "MacMon" = 2, "MK" = 2, "nCD4"=2, "Neu"=2)) %>%
  kable_styling(latex_options = c("scale_down")) %>%
  save_kable("../figures/superPIRTable.pdf")

#----------------------------------------------------------------------------------------------#
# make plots
#----------------------------------------------------------------------------------------------#

# used to add significance marks
maxRatio <- pir[, .("max"=max(Ratio)), by="CT"]
maxRatio$type <- "SIP"

ratioPlot <- ggplot(data=pir, aes(y=Ratio, x=CT, fill=type))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired", direction=1, labels = c("non-SIP", "SIP"))+ 
  geom_text(aes(label=format(Ratio, digits=2, drop0trailing = F)),
            vjust=-.2,
            position=position_dodge(width=0.9),
            color="black", size=3) +
  geom_point(data=maxRatio, aes(y=max+0.05), shape=8, color="red", show.legend=F)+
  ylab("Proportion of baits with super PIRs") + 
  xlab("Cell Type Group") +
  theme_bw() + theme(legend.title = element_blank())

# used to add significance marks
maxScore <- scores[, .("max"=max(max)), by="CT"]
maxScore$SIP <- 1

scorePlot <- ggplot(data=scores, aes(y=max, x=CT, fill=as.factor(SIP)))+
  geom_boxplot(outlier.alpha = 0.5)+
  geom_point(data=maxScore, aes(y=max+3), shape=8, color="red", show.legend=F)+
  scale_fill_brewer(palette="Paired", direction=1, labels = c("non-SIP", "SIP"))+ 
  ylab("PIR Score") + 
  xlab("Cell Type Group") + 
  theme_bw()+theme(legend.title = element_blank())

ggarrange(ratioPlot, scorePlot, 
         labels=letters[1:2], common.legend=T)
ggsave("../figures/superPIR.png", width=8, height=4)
ggsave("../figures/superPIR.pdf", width=8, height=4)

#----------------------------------------------------------------------------------------------#
# run nPIRplot.R to get ninter and medChicago plots
# combine all four plots
annotate_figure(ggarrange(ninter+xlab(element_blank()),
                          medChicago+xlab(element_blank()),
                          ratioPlot+xlab(element_blank()),
                          scorePlot+xlab(element_blank()),
                          nrow=2, ncol=2, common.legend = T,
                          labels=letters[1:4], align="hv"),
                bottom=text_grob("Cell Type"))

ggsave("../figures/Figure2_PIRs.png", width=8, height=8)
ggsave("../figures/Figure2_PIRs.pdf", width=8, height=8)

#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#

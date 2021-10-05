#----------------------------------------------------------------------------------------------#
# Figure 3c: SIP gene overlap with GWAS traits and randomly selected gene overlap with GWAS traits
# creates 
# requires *_specific_genes_allInfo.txt, SIP_*.txt, Chen data (5 files), and Vuckovic data
#----------------------------------------------------------------------------------------------#
library(data.table)
library(tidyr)
library(xlsx)
library(dplyr) # sample_n()
library(scales)

setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")
# setwd("/proj/yunligrp/users/lagler/SIP/data/")

# Vuckovic data (EUR)
vdat <- fread("ST3_ST4_BCX_Vuckovic_et_al.csv", drop=7)
vdat$ancestry <- "EUR"
# match phenotypes to chen
vdat[Phenotype=="RBC#", "Phenotype"] <- "RBC"
vdat[Phenotype=="WBC#", "Phenotype"] <- "WBC"
vdat[Phenotype=="NEUT#", "Phenotype"] <- "NEU"
vdat[Phenotype=="MONO#", "Phenotype"] <- "MONO"
vdat[Phenotype=="LYMPH#", "Phenotype"] <- "LYM"
vdat[Phenotype=="BASO#", "Phenotype"] <- "BASO"
vdat[Phenotype=="EO#", "Phenotype"] <- "EOS"
vdat[Phenotype=="PLT#", "Phenotype"] <- "PLT"

vdat2 <- vdat[Phenotype %in% c("HCT", "HGB", "MCH", "MCHC", "RBC", "RDW",
                           "MONO", "WBC", "PLT", "MPV", "LYM", "NEU") &
                Chr %in% 1:22]
colnames(vdat2)[6] <- "pval"
vdat2$pval <- 10^(-vdat2$pval) # match chen
vdat2$data <- "Vuckovic"

# Chen data has 5 ancestries Get union of SNPs
# /proj/yunligrp/users/jwen/BCX_meta_anal/GWASknownVars/Chen_*_2020.txt
files <- list.files(pattern="^Chen.*\\.txt$")
chen <- lapply(files, fread)
anc <- tstrsplit(files, "_")[[2]] %>% unlist()

# columns and column names not all the same 
chen[[1]] <- chen[[1]][, c(1:4,11,15)]
colnames(chen[[1]])[1] <- "Pheno"
chen[[2]] <- chen[[2]][, c(1:4,11,15)]
chen[[3]] <- chen[[3]][, c(1:4,11,15)]
chen[[4]] <- chen[[4]][, c(1:4,11,14)]
colnames(chen[[4]])[4] <- "POS_(hg19)"
chen[[5]] <- chen[[5]][, c(1:4,9,15)]
# add ancestry
for(i in 1:5){
  chen[[i]]$ancestry <- anc[i]
  chen[[i]]$`P-value` <- as.numeric(chen[[i]]$`P-value`)
}
# merge all CTs together
chenall <- Reduce(function(...) merge(..., all = TRUE), chen)
colnames(chenall) <- c("Phenotype", "VariantID", "Chr", "Position", "pval",
                       "rsID", "ancestry")
chenall$data <- "Chen"

# combine sets 
gwasLong <- rbind(vdat2, chenall)
gwas <- gwasLong
gwas$Chr <- as.numeric(gwas$Chr)
# long to wide, reduce to union 
# gwas <- dcast(unique(gwasLong),
#               Phenotype + VariantID + rsID + Chr + Position +
#                  + data + `pval`~ ancestry ,
#                   value.var = "ancestry", fun.aggregate=NULL) %>% 
#   unite(.,ancestry, 8:13, sep=";", na.rm=T)

# pcHiC data formatted for intra-chromosomal interactions and distance < 2MB
# /proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/
load("fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[[31]] <- "MacMon"
fraser$baitChr <- as.numeric(fraser$baitChr)

#----------------------------------------------------------------------------------------------#
# function to find SNP overlaps by cell type and variant
#----------------------------------------------------------------------------------------------#

getSNPs <- function(celltype, phenotype, windowsize,
                    CTspecific=TRUE, nSample=NA, return.data=FALSE, nSNP=1){
  
  if(nSNP==1){
    gwas_ct <- gwas[Phenotype %in% phenotype,]
  }else{
    gwas_ct <- gwas[Phenotype %in% phenotype & ancestry == "EUR",]
  }
  
  setDT(gwas_ct)[, rsStart := Position - windowsize]
  setDT(gwas_ct)[, rsEnd := Position + windowsize]
  
  setkey(gwas_ct, Chr, rsStart, rsEnd)
  
  # SIP cumulative interaction score
  # /proj/yunligrp/users/lagler/SIP/data/
  sipScore <- fread(paste0("SIP_", celltype, "_specific.txt"), drop=3)
  colnames(sipScore)[[4]] <- "Specific"
  setorder(sipScore, -intscore)
  sipScore$rank <- seq(1:nrow(sipScore))
  
  # combine fraser with sip and filter out non signficant interactions
  fraser_ct <- fraser[eval(as.name(celltype)) >= 5,]
  fraser_ct <- fraser_ct[, c(1:10)]
  sip <- merge(sipScore, 
               fraser_ct,
               by="baitID")
  sip$oeChr <- as.numeric(sip$oeChr)
  
  # sample same number as CT specific SIP genes unless provided
  # if not sampling, look at only SIPs or cell type specific SIPs
  if(is.character(nSample)){
    if(CTspecific==TRUE){
      sipSample <- sip[Specific==1]}
    else{
      if(CTspecific==FALSE){sipSample <- sip[SIP==1]}
    }
  }else{
    if(is.na(nSample)){
      if(CTspecific==TRUE){
        nsample <- uniqueN(sip[Specific==1], "baitID")
      }else{
        if(CTspecific==FALSE)
        nsample <- uniqueN(sip[SIP==1], "baitID")
      }
    }else{
      nsample <- nSample
      }
    # select n baits
    baitSample <- sample(unique(sip$baitID), nsample)
    # sample data
    sipSample <- sip[baitID %in% baitSample]
  }
  
  # find overlaps of PIRs and SNPs
  # some PIRs may have multiple SNPs
  setkey(sipSample, oeChr, oeStart, oeEnd)
  overlaps <- foverlaps(sipSample, gwas_ct, nomatch = NULL)
  
  # count how many baits have at least one PIR overalpping a SNP
  # specific columns meaningless unless looking at SIPs
  counts <- data.table("CT" = celltype,
                       "PT" = paste(phenotype, collapse = ", "),
                       "nSNP" = uniqueN(overlaps, by="baitID"),
                       "nBait" = uniqueN(sipSample, by="baitID"))
  setDT(counts)[, percent := nSNP/nBait]
  #print(counts)
  if(return.data==TRUE){
    setorder(overlaps, rank)
    if(CTspecific==TRUE){
      return(overlaps[,c("baitID", "baitName", "baitChr", "baitStart", "baitEnd",
                         "oeStart", "oeEnd", "oeID", "oeName", 
                         "Phenotype", "VariantID", "rsID", "Position", "pval",
                         "ancestry", "data")])
    }else{
      return(overlaps[,c("baitID", "baitName", "baitChr", "baitStart", "baitEnd",
                         "oeStart", "oeEnd", "oeID", "oeName", "Specific",
                         "Phenotype", "VariantID", "rsID", "Position", "pval",
                         "ancestry", "data")])
    }

  }else{
    return(counts)
  }
}

getTable <- function(celltype, phenotype){
  ss <- getSNPs(celltype, phenotype, 0, TRUE, "SIP") # specific SIPS
  as <- getSNPs(celltype, phenotype, 0, FALSE, "SIP") # all SIPS
  
  # do sample 100 times
  sr100 <- list(); ar100 <-list()
  for(i in 1:100){
    sr100[[i]] <- getSNPs(celltype, phenotype, 0, TRUE, NA) # specific RANDOM
    ar100[[i]] <- getSNPs(celltype, phenotype, 0, FALSE, NA) # all RANDOM
  }
  sr100 <- rbindlist(sr100); ar100 <- rbindlist(ar100)
  sr <- sr100[1,]; ar <- ar100[1,]
  sr$nSNP <- median(sr100$nSNP); ar$nSNP <- median(ar100$nSNP)
  setDT(sr)[, percent := nSNP/nBait]; setDT(ar)[, percent := nSNP/nBait]
  
  out <- rbind(ss, sr, as, ar)
  out$type <- c("spec.SIP", "spec.Random", "all.SIP", "all.Random")
  
  # fisher's exact test
  fisher_s <- fisher.test(out[1:2, 3:4])
  fisher_a <- fisher.test(out[3:4, 3:4])
  out$OR <- c(round(fisher_s$estimate,2), "--", round(fisher_a$estimate,2), "--")
  out$pval <- c(format(fisher_s$p.value,digits=2, scientific=T), "--",
                format(fisher_a$p.value,digits=2, scientific=T), "--")
  out$lower <- c(round(fisher_s$conf.int[1],2), "--",
                 round(fisher_a$conf.int[1],2), "--")
  out$upper <- c(round(fisher_s$conf.int[2],2), "--",
                 round(fisher_a$conf.int[2],2), "--")
  out$percent <- format(out$percent*100, digits=2)
  return(out)
}

# red blood cell related traits
RBC <- c("HCT", "HGB", "MCH", "MCHC", "RBC", "RDW")

set.seed(12162020)
ery <- getTable("Ery", RBC)
macmon <- getTable("MacMon", c("MONO", "WBC"))
mk <- getTable("MK", c("PLT", "MPV"))
ncd4 <- getTable("nCD4", c("LYM", "WBC"))
neu <- getTable("Neu", c("NEU", "WBC"))

# combine cell types
alltables <- rbind(ery, macmon, mk, ncd4, neu)
alltables$OR <- as.numeric(alltables$OR)
alltables$pval <- as.numeric(alltables$pval)
alltables$lower <- as.numeric(alltables$lower)
alltables$upper <- as.numeric(alltables$upper)

# output file for easier plotting later
fwrite(alltables, "GWASannotPIR_tables.txt", sep="\t", col.names=T)

alltables <- fread("GWASannotPIR_tables.txt")
plotdata <- alltables[is.na(pval)==F,]
plotdata$CT <- factor(plotdata$CT, levels=c("Neu", "nCD4", "MK", "MacMon", "Ery"))

allsipOR <- ggplot(plotdata[type=="all.SIP"], aes(x = OR, y = CT)) + 
  #geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = upper, xmin = lower),
                 size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "purple3") +
  # facet_wrap(~type, scales="free", ncol=1,
  #            labeller=labeller(type=c(all.SIP="SIP", spec.SIP="Specific SIP"))) +
  theme_bw()+ theme(panel.grid.minor = element_blank()) +
  ylab("Cell Type Group") + xlab("Odds Ratio for Variant Overlap") 
allsipOR
ggsave("../figures/gwasPIR_ORs_allSIP.png", width=6, height=4, unit="in")

# supplemental figure 4: cell type specific SIPs
ggplot(plotdata[type=="spec.SIP"], aes(x = OR, y = CT)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = upper, xmin = lower),
                 size = .5, height =.2, color = "gray50") +
  geom_point(size = 3.5, color = "purple3") +
  scale_x_continuous(breaks=seq(1,11,2))+
  # facet_wrap(~type, scales="free", ncol=1,
  #            labeller=labeller(type=c(all.SIP="SIP", spec.SIP="Specific SIP"))) +
  theme_bw()+  theme(panel.grid.minor = element_blank()) +
  ylab("Cell Type") + xlab("Odds Ratio") 
ggsave("../figures/FigureS4_gwasPIRspecific.png", width=6, height=4, unit="in")
ggsave("../figures/FigureS4_gwasPIRspecific.pdf", width=6, height=4, unit="in")

#--------------------------------------------------------------------------------------#

# arrange ATAC-seq plots with GWAS plot
# common y axis
annotate_figure(ggarrange(propPlot+xlab(element_blank()),
                          nOvlpPlot+xlab(element_blank()),
                          allsipOR+ylab(element_blank()),
                          align="hv", common.legend = T,
                          nrow=1, labels=letters[1:3],
                          widths = c(.3,.4,.3)),
                left = text_grob("Cell Type", rot=90))
ggsave("../figures/Figure3_ATACor.png", height=4, width=10, unit="in")
ggsave("../figures/Figure3_ATACor.pdf", height=4, width=10, unit="in")

#--------------------------------------------------------------------------------------#
# cell type specific SIP baits with at least one PIR overlapping a SNP  
ery.spec <- getSNPs("Ery", RBC, 0, TRUE, "SIP", TRUE)  
macmon.spec <- getSNPs("MacMon", c("MONO", "WBC"), 0, TRUE, "SIP", TRUE) 
mk.spec <- getSNPs("MK", c("PLT", "MPV"), 0, TRUE, "SIP", TRUE) 
ncd4.spec <- getSNPs("nCD4", c("LYM", "WBC"), 0, TRUE, "SIP", TRUE) 
neu.spec <- getSNPs("Neu", c("NEU", "WBC"), 0, TRUE, "SIP", TRUE) 

# find unique baitIDs and unique genes
getunique <- function(ct.spec, ct){
  ct.unique <- unique(ct.spec, by="baitID")[,c("baitName", "baitID", "rank", "intscore")]
  # one gene per row
  ct.genes <- separate_rows(ct.unique, baitName, sep=";") %>% data.table()
  #NOTE: macmon has no ct specific genes overlapping a SNP
  # cell type specific genes 
   spec.genes <- fread(paste0(ct, "_specific_genes_allInfo.txt"))
  # # subset specific SIPs to only be specific genes
  ct.genes <- ct.genes[baitName %in% spec.genes$baitName]
  ct.genes$CT <- ct
  return(ct.genes)
}

annotgenes <- rbind(getunique(ery.spec, "Ery")[c(1:5),],
                    getunique(macmon.spec, "MacMon")[1:5,],
                    getunique(mk.spec, "MK")[c(2,3,4,5,7,8),],
                    getunique(ncd4.spec, "nCD4")[c(1,3,4,6,8),],
                    getunique(neu.spec, "Neu")[c(1,3,5,6,8),])
fwrite(annotgenes, "annotgenes.txt", sep="\t", col.names = )

#--------------------------------------------------------------------------------------#
# write GWAS PIR overlaps to excel file
# get data
ery.all <- getSNPs("Ery", RBC, 0, FALSE, "SIP", TRUE)  
macmon.all <- getSNPs("MacMon", c("MONO", "WBC"), 0, FALSE, "SIP", TRUE) 
mk.all <- getSNPs("MK", c("PLT", "MPV"), 0, FALSE, "SIP", TRUE) 
ncd4.all <- getSNPs("nCD4", c("LYM", "WBC"), 0, FALSE, "SIP", TRUE) 
neu.all <- getSNPs("Neu", c("NEU", "WBC"), 0, FALSE, "SIP", TRUE) 
# write to file
xlsx::write.xlsx(ery.all, file=paste0("SIP_PIR_GWAS.xlsx"),
                 sheetName="Ery_SIPs", append=T, row.names=F)
xlsx::write.xlsx(macmon.all, file=paste0("SIP_PIR_GWAS.xlsx"),
                 sheetName="MacMon_SIPs", append=T, row.names=F)
xlsx::write.xlsx(mk.all, file=paste0("SIP_PIR_GWAS.xlsx"),
                 sheetName="MK_SIPs", append=T, row.names=F)
xlsx::write.xlsx(ncd4.all, file=paste0("SIP_PIR_GWAS.xlsx"),
                 sheetName="nCD4_SIPs", append=T, row.names=F)
xlsx::write.xlsx(neu.all, file=paste0("SIP_PIR_GWAS.xlsx"),
                 sheetName="Neu_SIPs", append=T, row.names=F)

#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#


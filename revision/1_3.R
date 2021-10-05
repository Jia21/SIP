##########################################################################
### Revised from /proj/yunligrp/users/lagler/SIP/R/GWASannotationPIR.R ###
###                   By Jia Wen 09/2021                               ###
##########################################################################

setwd("/proj/yunligrp/users/jwen/SIP/revise")
# pcHiC data formatted for intra-chromosomal interactions and distance < 2MB
# /proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/fraser.rda
library(data.table)
load("/proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[[31]] <- "MacMon"

summary(fraser$dist)
len_bait <- data.frame(fraser$baitEnd - fraser$baitStart + 1)
colnames(len_bait) <- "length"
len_bait$label <- "bait"
theme <- theme(axis.text.x = element_text(size = 15),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text( size=18, face="bold"),
        legend.title = element_text(color = "black", size = 15),
        legend.text = element_text(color = "black", size = 15),
        axis.title.y = element_text( size=18, face="bold")) 
library(ggplot2)
p <- ggplot(len_bait, aes(x=label,y=length)) + 
  geom_violin() +theme
ggsave("bait_length.png", height=7, width=10.5)
#     bait length                    
#  Min.   :   123     
#  1st Qu.:  3464  
#  Median :  5912   
#  Mean   :  7221                     
#  3rd Qu.:  9456                     
#  Max.   :156451    

fraser$bait_length <- fraser$baitEnd - fraser$baitStart + 1
fraser$oe_length <-  fraser$oeEnd - fraser$oeStart + 1
# y <- distinct(subset(fraser,select=c("baitID","bait_length")))

mat <- matrix(0,nrow=5,ncol=5)
CT <- c("Ery","MK","MacMon","Neu","nCD4")
# mat[,1] <- CT
mat <- NULL
p <- NULL
for (i in 1:length(CT)){
    cat(CT[i],"\n")
    sips <- read.table(paste0("/proj/yunligrp/users/lagler/SIP/data/SIP_",CT[i],".txt"),header = T) %>% filter(SIP==1)   

    # y$label <- ifelse(as.matrix(y$baitID) %in% as.matrix(sips$baitID),1,0)

    # y$kb <- y$bait_length/1000

    # mod1 <- glm(label ~ bait_length , data = y, family = "binomial")
    # mat[i,2] <- formatC(mod1$coefficients["bait_length"],format= "e",digits=2)
    # mat[i,3] <- formatC(summary(mod1)$coefficients["bait_length","Pr(>|z|)"],format= "e",digits=2)

    # mod2 <- glm(label ~ kb , data = y, family = "binomial")
    # mat[i,4] <- formatC(mod2$coefficients["kb"],format= "e",digits=2)
    # mat[i,5] <- formatC(summary(mod2)$coefficients["kb","Pr(>|z|)"],format= "e",digits=2)

    fraser_ct <- fraser[eval(as.name(CT[i])) >= 5,]
    sip_pir <- fraser_ct[fraser_ct$baitID %in% sips$baitID,] %>% select(oeChr,oeStart,oeEnd,oeID) %>% distinct()
    sir_pir <- data.frame(sip_pir)
    sip_pir$length <- sip_pir$oeEnd - sip_pir$oeStart + 1
    sip_pir$ct <- CT[i]
    sip_pir$label <- "SIP PIR"

    nonsip_pir <- fraser_ct[!fraser_ct$baitID %in% sips$baitID,] %>% select(oeChr,oeStart,oeEnd,oeID) %>% distinct()
    nonsip_pir <- data.frame(nonsip_pir)
    nonsip_pir$length <- nonsip_pir$oeEnd - nonsip_pir$oeStart + 1
    nonsip_pir$ct <- CT[i]
    nonsip_pir$label <- "non-SIP PIR"
    mat <- rbind(mat,rbind(sip_pir,nonsip_pir))
    p <- c(p,wilcox.test(abs(sip_pir$length/1000),abs(nonsip_pir$length/1000),alternative = "less")$p.value)
}
# colnames(mat) <- c("CT","beta_bp","pval_bp","beta_kb","pval_kb")

library(ggplot2)
plotd <- ggplot(mat, aes(x=label, y = length/1000, fill=label)) +
        geom_violin() + scale_fill_brewer(palette="Paired") + ylab("Length of SIP/non-SIP PIRs(KB)") + xlab("") +
         geom_point(data=mat[, .("median"=median(abs(length/1000))), by="label"], aes(y=median), show.legend = F,size = 5)+ facet_wrap(.~ct) +
         theme_bw() +
        theme(axis.text=element_text(size=20),axis.title=element_text(size=20)) +
        theme(legend.title = element_blank(),legend.text = element_text(size=20))+
        theme(strip.text.x = element_text(size = 15))
     ggsave(paste0("./1_3/","PIRs_violin.png"),width=16,height=10)
    ggsave(paste0("./1_3/","PIRs_violin.pdf"),width=16,height=10)


CT <- c("Ery","MK","MacMon","Neu","nCD4")
# mat[,1] <- CT
mat <- NULL
p <- NULL
for (i in 1:length(CT)){
    cat(CT[i],"\n")
    sips <- read.table(paste0("/proj/yunligrp/users/lagler/SIP/data/SIP_",CT[i],".txt"),header = T) %>% filter(SIP==1)   

    # y$label <- ifelse(as.matrix(y$baitID) %in% as.matrix(sips$baitID),1,0)

    # y$kb <- y$bait_length/1000

    # mod1 <- glm(label ~ bait_length , data = y, family = "binomial")
    # mat[i,2] <- formatC(mod1$coefficients["bait_length"],format= "e",digits=2)
    # mat[i,3] <- formatC(summary(mod1)$coefficients["bait_length","Pr(>|z|)"],format= "e",digits=2)

    # mod2 <- glm(label ~ kb , data = y, family = "binomial")
    # mat[i,4] <- formatC(mod2$coefficients["kb"],format= "e",digits=2)
    # mat[i,5] <- formatC(summary(mod2)$coefficients["kb","Pr(>|z|)"],format= "e",digits=2)

    fraser_ct <- fraser[eval(as.name(CT[i])) >= 5,]
    sip_pir <- fraser_ct[fraser_ct$baitID %in% sips$baitID,] %>% select(baitChr,baitStart,baitEnd,baitID) %>% distinct()
    sir_pir <- data.frame(sip_pir)
    sip_pir$length <- sip_pir$oeEnd - sip_pir$oeStart + 1
    sip_pir$ct <- CT[i]
    sip_pir$label <- "SIP PIR"

    nonsip_pir <- fraser_ct[!fraser_ct$baitID %in% sips$baitID,] %>% select(oeChr,oeStart,oeEnd,oeID) %>% distinct()
    nonsip_pir <- data.frame(nonsip_pir)
    nonsip_pir$length <- nonsip_pir$oeEnd - nonsip_pir$oeStart + 1
    nonsip_pir$ct <- CT[i]
    nonsip_pir$label <- "non-SIP PIR"
    mat <- rbind(mat,rbind(sip_pir,nonsip_pir))
    p <- c(p,wilcox.test(abs(sip_pir$length/1000),abs(nonsip_pir$length/1000),alternative = "less")$p.value)
}
# colnames(mat) <- c("CT","beta_bp","pval_bp","beta_kb","pval_kb")

library(ggplot2)
plotd <- ggplot(mat, aes(x=label, y = length/1000, fill=label)) +
        geom_violin() + scale_fill_brewer(palette="Paired") + ylab("Length of SIP/non-SIP PIRs(KB)") + xlab("") +
         geom_point(data=mat[, .("median"=median(abs(length/1000))), by="label"], aes(y=median), show.legend = F,size = 5)+ facet_wrap(.~ct) +
         theme_bw() +
        theme(axis.text=element_text(size=20),axis.title=element_text(size=20)) +
        theme(legend.title = element_blank(),legend.text = element_text(size=20))+
        theme(strip.text.x = element_text(size = 15))
     ggsave(paste0("./1_3/","PIRs_violin.png"),width=16,height=10)
    ggsave(paste0("./1_3/","PIRs_violin.pdf"),width=16,height=10)


write.table(mat,"logis_bait_len_sip",quote = F,sep="\t",col.names = T,row.names = F)

# len_PIR <-  data.frame(fraser$oeEnd - fraser$oeStart + 1)
# colnames(len_PIR) <- "length"
# len_PIR$label <- "oe"
# len <- rbind(len_bait,len_PIR)
# p <- ggplot(len, aes(x=label,y=length,color=label)) + 
#   geom_violin() +theme
# ggsave("bait_PIR_length.png", height=7, width=10.5)

# p <- ggplot(len_PIR, aes(x=label,y=length)) + 
#   geom_violin() +theme
# ggsave("PIR_length.png", height=7, width=10.5)
# #   length               
#  Min.   :  151      
#  1st Qu.: 1134   
#  Median : 2471   
#  Mean   : 3609                     
#  3rd Qu.: 4838                     
#  Max.   :39767    


### If using SCZ GWAS study ###
library(data.table)
library(tidyr)
library(xlsx)
library(dplyr) # sample_n()
library(scales)
library(vroom)

sczGWAS <- vroom("/proj/yunligrp/users/jwen/GWAS_sumstat/sullivan-pardinas2018/CLOZUK_PGC2noclo.METAL.assoc.dosage.fix",col_names = T)
sczGWAS_sigVars <- sczGWAS %>% filter(P < 5e-8) %>% dplyr::select(SNP, CHR, BP)
sczGWAS_sigVars$start <- sczGWAS_sigVars$BP-1
sczGWAS_sigVars$end <- sczGWAS_sigVars$BP
sczGWAS_sigVars$Phenotype <- "SCZ"
sczGWAS_sigVars$ancestry <- "EUR"
gwas <- sczGWAS_sigVars
# pcHiC data formatted for intra-chromosomal interactions and distance < 2MB
# /proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/
load("/proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[[31]] <- "MacMon"
fraser$baitChr <- as.numeric(fraser$baitChr)

#----------------------------------------------------------------------------------------------#
# function to find SNP overlaps by cell type and variant
#----------------------------------------------------------------------------------------------#

# nSample denote SIP or non-SIPs or randomly shuffle if nSample = NA #
getSNPs <- function(celltype, phenotype, windowsize,
                    CTspecific=TRUE, nSample=NA, return.data=FALSE, nSNP=1){
  
  if(nSNP==1){
    gwas_ct <- sczGWAS_sigVars
  }else{
    gwas_ct <- sczGWAS_sigVars[Phenotype %in% phenotype & ancestry == "EUR",]
  }
  
  setDT(gwas_ct)[, rsStart := BP - windowsize]
  setDT(gwas_ct)[, rsEnd := BP + windowsize]
  
  setkey(gwas_ct, CHR, rsStart, rsEnd)
  
  # SIP cumulative interaction score
  # /proj/yunligrp/users/lagler/SIP/data/
  sipScore <- fread(paste0("/proj/yunligrp/users/lagler/SIP/data/SIP_", celltype, "_specific.txt"), drop=3)
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
  sr100 <- list() 
  ar100 <-list()
  for(i in 1:100){
    cat(i,"\t")
    sr100[[i]] <- getSNPs(celltype, phenotype, 0, TRUE, NA) # specific RANDOM
    ar100[[i]] <- getSNPs(celltype, phenotype, 0, FALSE, NA) # all RANDOM
  }
  sr100 <- rbindlist(sr100)
  ar100 <- rbindlist(ar100)
  sr <- sr100[1,]
  ar <- ar100[1,]
  sr$nSNP <- median(sr100$nSNP)
  ar$nSNP <- median(ar100$nSNP)
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

# SCZ GWAS
phenptype <- "SCZ"

set.seed(12162020)
ery <- getTable("Ery", phenptype)
macmon <- getTable("MacMon", phenptype)
mk <- getTable("MK", phenptype)
ncd4 <- getTable("nCD4", phenptype)
neu <- getTable("Neu", phenptype)

# combine cell types
alltables <- rbind(ery, macmon, mk, ncd4, neu)
alltables$OR <- as.numeric(alltables$OR)
alltables$pval <- as.numeric(alltables$pval)
alltables$lower <- as.numeric(alltables$lower)
alltables$upper <- as.numeric(alltables$upper)

# output file for easier plotting later
fwrite(alltables, "sczGWASannotPIR_tables.txt", sep="\t", col.names=T)

alltables <- fread("sczGWASannotPIR_tables.txt")
plotdata <- alltables[is.na(pval)==F,]
plotdata$CT <- factor(plotdata$CT, levels=c("Neu", "nCD4", "MK", "MacMon", "Ery"))

allsipOR <- ggplot(plotdata[type=="spec.SIP"], aes(x = OR, y = CT)) + 
  #geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = upper, xmin = lower),
                 size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "purple3") +
  # facet_wrap(~type, scales="free", ncol=1,
  #            labeller=labeller(type=c(all.SIP="SIP", spec.SIP="Specific SIP"))) +
  theme_bw()+ theme(panel.grid.minor = element_blank()) +
  ylab("Cell Type Group") + xlab("Odds Ratio for Variant Overlap") 
allsipOR
ggsave("sczGWASPIR_ORs_specSIP.png", width=6, height=4, unit="in")

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
ggsave("FigureS4_sczGWASPIRspecific.png", width=6, height=4, unit="in")
ggsave("FigureS4_sczGWASPIRspecific.pdf", width=6, height=4, unit="in")

#--------------------------------------------------------------------------------------#
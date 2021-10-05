#----------------------------------------------------------------------------------------------#
# Supplemental File 1: Details of SIP subnetworks
# creates 
# requires SIP_*.txt, ST3_ST4_BCX_Vuckovic_et_al.csv, fraser.rda, and atac-seq data files
#----------------------------------------------------------------------------------------------#
library(data.table)
library(tidyr)
library(xlsx)

setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")

# Vuckovic data (EUR)
vdat <- fread("ST3_ST4_BCX_Vuckovic_et_al.csv") 
# match phenotypes to chen
vdat[Pheno=="RBC#", "Pheno"] <- "RBC"
vdat[Pheno=="WBC#", "Pheno"] <- "WBC"
vdat[Pheno=="NEUT#", "Pheno"] <- "NEU"
vdat[Pheno=="MONO#", "Pheno"] <- "MONO"
vdat[Pheno=="LYMPH#", "Pheno"] <- "LYM"
vdat[Pheno=="BASO#", "Pheno"] <- "BASO"
vdat[Pheno=="EO#", "Pheno"] <- "EOS"
vdat[Pheno=="PLT#", "Pheno"] <- "PLT"

# pcHiC data formatted for intra-chromosomal interactions and distance < 2MB
# /proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/f
load("fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[31] <- "MacMon"

# ATAC-seq data
setwd("/Users/tmlagler/OneDrive/Lab/SIP/ATACseq/")
#setwd("/proj/yunligrp/users/jwen/ATAC-seq_collection/ATAC-seq/2019-08/")
atac_ery <- fread("Ery_peaks.narrowPeak.gz", select=c(1:3))
atac_mon <- fread("Mono_peaks.narrowPeak.gz", select=c(1:3))
atac_ncd4 <- fread("CD4_peaks.narrowPeak.gz", select=c(1:3))
atac_mk <- fread("Mega_peaks.narrowPeak.gz", select=c(1:3))

#----------------------------------------------------------------------------------------------#
# function to find SNP overlaps by cell type and variant
#----------------------------------------------------------------------------------------------#

setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")
# /proj/yunligrp/users/lagler/SIP/data/

getSNPs <- function(celltype, phenotype, windowsize, seqdat){
  
  if(phenotype=="all"){
    gwas <- vdat
  }else{
    gwas <- vdat[Pheno %in% phenotype,]
  }
  
  setDT(gwas)[, rsStart := `POS_(hg19)` - windowsize]
  setDT(gwas)[, rsEnd := `POS_(hg19)` + windowsize]
  
  setkey(gwas, CHR, rsStart, rsEnd)
  
  # SIP cumulative interaction score
  sipScore <- fread(paste0("SIP_", celltype, ".txt"), drop=3)
  setorder(sipScore, -intscore)
  sipScore$rank <- seq(1:nrow(sipScore))
  
  # combine fraser with sip and filter out non signficant interactions
  fraser_ct <- fraser[eval(as.name(celltype)) >= 5,]
  fraser_ct <- fraser_ct[, c(1:10)]
  sip <- merge(sipScore, 
               fraser_ct,
               by="baitID")
  sip$oeChr <- as.numeric(sip$oeChr)
  
  # find overlaps of PIRs and SNPs
  # some PIRs may have multiple SNPs
  setkey(sip, oeChr, oeStart, oeEnd)
  overlaps <- foverlaps(sip, gwas, nomatch = NULL)
  
  if(sum(is.na(seqdat))==1){
    sipOverlaps <- overlaps[SIP==1]
  }else{
    # format ATAC-seq data
    colnames(seqdat) <- c("chr", "start", "end")
    seqdat <- seqdat[chr %in% paste0("chr", 1:22)] 
    seqdat$chrNum <- as.numeric(tstrsplit(seqdat$chr, "chr")[[2]])
    setkey(seqdat, chrNum, start, end)
    
    # PIR overlapping GWAS SNP and ATAC-seq peak
    atac <- foverlaps(overlaps, seqdat, by.x=c("oeChr", "oeStart", "oeEnd"),
                      type="any", mult="first", nomatch=NULL)
    sipOverlaps <- atac[SIP==1]
  }

  
  nPIRs <- unique(sipOverlaps, by=c("baitID", "oeID")) %>%
    count(., baitID, sort=T, name="nPIRs") %>% data.table()
  potential <- merge(nPIRs, sipScore, by="baitID")
  setorder(potential, -nPIRs, rank)
  baits <- potential[nPIRs >= 2]
  
  # find which networks are within 500kb
  baitInfo <- sipOverlaps[baitID %in% baits$baitID]
  baitInfo2 <- baitInfo[, .("minPos"=min(oeStart),
                            "maxPos"=max(oeEnd)), by="baitID"]
  setDT(baitInfo2)[, netDist := maxPos - minPos]
  # remove network constraint
  #baits2 <- baitInfo2[netDist <= 500000]
  baits2 <- baitInfo2[netDist <= Inf]
  
  outData <- merge(sipOverlaps, baits2, by="baitID") %>%
    merge(., baits[,1:2], by="baitID")
  
  # add nSNPs per subnetwork
  nSNPs <- count(data.frame(outData), baitID, sort=T, name="nSNPs") %>% data.table()
  outData2 <- merge(outData, nSNPs, by="baitID", all.x=T)
  

  print(uniqueN(outData2, by="baitID"))
  # reorder and return
  outData2 <- outData2[, c("baitID", "baitName", "baitChr", "baitStart",
                          "baitEnd", "oeStart", "oeEnd", "oeID", "oeName",
                          "MarkerName", "rsID", "POS_(hg19)", "nPIRs", "nSNPs")]
  outData2$baitChr <- as.numeric(outData2$baitChr) #nicer excel printing
  setorder(outData2, -nPIRs, baitID)
  return(outData2)
}

#----------------------------------------------------------------------------------------------#
# call function for each cell type/phenotype combination
ery_hct <- getSNPs("Ery", "HCT", 0, atac_ery) #2
ery_hgb <- getSNPs("Ery", "HGB", 0, atac_ery) #2
ery_mch <- getSNPs("Ery", "MCH", 0, atac_ery) #7
ery_mchc <- getSNPs("Ery", "MCHC", 0, atac_ery) #3
ery_rbc <- getSNPs("Ery", "RBC", 0, atac_ery) #4
ery_rdw <- getSNPs("Ery", "RDW", 0, atac_ery) #11

macmon_mono <- getSNPs("MacMon", "MONO", 0, atac_mon) #5
macmon_wbc <- getSNPs("MacMon", "WBC", 0, atac_mon) #1

mk_plt <- getSNPs("MK", "PLT", 0, atac_mk) #14
mk_mpv <- getSNPs("MK", "MPV", 0, atac_mk) #10

ncd4_lym <- getSNPs("nCD4", "LYM", 0, atac_ncd4) #15
ncd4_wbc <- getSNPs("nCD4", "WBC", 0, atac_ncd4) #2

# no Neu ATAC-seq data
neu_neu <- getSNPs("Neu", "NEU", 0, NA) #16
neu_wbc <- getSNPs("Neu", "WBC", 0, NA) #22

#----------------------------------------------------------------------------------------------#
# write to excel file
library(xlsx)
writeSheet <- function(ct_trait){
  xlsx::write.xlsx(ct_trait, file="SIPsubnetworks.xlsx",
                   sheetName=toupper(deparse(substitute(ct_trait))),
                   append=T, row.names=F)
}

writeSheet(ery_hct); writeSheet(ery_hgb); writeSheet(ery_mch)
writeSheet(ery_mchc); writeSheet(ery_rbc); writeSheet(ery_rdw)
writeSheet(macmon_mono); writeSheet(macmon_wbc)
writeSheet(mk_plt); writeSheet(mk_mpv)
writeSheet(ncd4_lym); writeSheet(ncd4_wbc)
writeSheet(neu_neu); writeSheet(neu_wbc)
#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#

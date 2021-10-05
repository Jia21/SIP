#----------------------------------------------------------------------------------------------#
# Get PIR regions of cell type specific SIP genes for LDSC Regression
# outputs bed file of PIR locations for LDSC regression 
# requires SIP_*_specific.txt, fraser.rda
#----------------------------------------------------------------------------------------------#

library(data.table)
library(dplyr)
setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")

# pcHiC data formatted for intra-chromosomal interactions and distance < 2MB
# /proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/fraser.rda
load("fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[31] <- "MacMon"

# SIP call files
# /proj/yunligrp/users/lagler/SIP/data/
files <- list.files(pattern="^SIP.*\\.txt$")
files <- files[grep("specific", files)] # only include specific files
data <- lapply(files, fread) 
cells <- tstrsplit(files, "_")[[2]] %>% tstrsplit(., ".txt") %>% unlist()

#----------------------------------------------------------------------------------------------#

# function to get PIR positions for each cell type
getPIRs <- function(sipdata, ct){
  # merge sip data with pcHiC to get PIRs
  sipdata <- merge(sipdata, fraser, by="baitID",
                   all.x=T, all.y=F) %>% data.table()
 
   # only keep significant interactions
  sipdataPIRs <- sipdata[eval(as.name(ct)) > 5] %>%
    unique(., by="oeID") %>%
    dplyr::select(.,c("oeChr", "oeStart", "oeEnd")) %>% data.table()
  
  setDT(sipdataPIRs)[, chr := paste0("chr", oeChr)]
  return(sipdataPIRs[,c("chr", "oeStart", "oeEnd")])
}

for(i in 1:5){
  # store sip call data and cell type
  ctData <- data[[i]]
  ct <- cells[i]
  
  #  only SIPs
  sips <- ctData[SIP==1]
  
  # only cell type-specific SIPs
  specCol <- colnames(ctData)[5]
  spec.sips <- ctData[eval(as.name(specCol))==1]
  
  # get positions of unique PIRs
  sipPIRs <- getPIRs(sips , ct)
  specPIRs <- getPIRs(spec.sips, ct)
  
  # write to file
  # /proj/yunligrp/users/lagler/SIP/LDSC/PIRregions/
  write.table(sipPIRs, paste0("all_sipPIRs_", ct, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)
  write.table(specPIRs, paste0("specific_sipPIRs_", ct, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)
}


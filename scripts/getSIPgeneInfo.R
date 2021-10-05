#----------------------------------------------------------------------------------------------#
# Supplemental File 2
# creates lists of cell type group-specific genes and shared genes
# outputs *_specific_genes.txt for LDSC regression (only ENSEMBL IDs)
#     *_specific_genes_allInfo.txt (gene name, ID, baitID, position, and gene expression)
#     allSIPs_geneInfo.txt (all genes, SIP or non-SIP)
# requires SIP_*.txt, fraser.rda, GeneExpressionMatrix_combine_Mac_Mon.txt
# input for distancePlot.R and LDSC regression
# Updated 7/31/21 to include SIPs without baitNames
#----------------------------------------------------------------------------------------------#
library(data.table)
library(dplyr)
library(tidyr) # separate rows
library(Biobase) # get gene names from ensembl ids
library(biomaRt)
library(xlsx)
setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/data/")
#----------------------------------------------------------------------------------------------#
# read in data 
#----------------------------------------------------------------------------------------------#

# pcHiC data formatted for intra-chromosomal interactions and distance < 2MB
# /proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/fraser.rda
load("fraser.rda")
fraser <- data.table(fraser)
fraser <- fraser[baitChr %in% 1:22]
colnames(fraser)[31] <- "MacMon"

# SIP call files
# /proj/yunligrp/users/lagler/SIP/data/SIP_*.txt
files <- list.files(pattern="^SIP.*\\.txt$")
files <- files[-grep("specific", files)] # don't include specific files
data <- lapply(files, fread) 
cells <- tstrsplit(files, "_")[[2]] %>% tstrsplit(., ".txt") %>% unlist()

# change SIP column to by CT, remove unnecessary columns
for(i in 1:5){
  colnames(data[[i]])[4] <- paste0("SIP.", cells[i])
  data[[i]]$intscore <- NULL
  data[[i]]$rank <- NULL
}

# merge all CTs together
mergedData <- Reduce(function(...) merge(..., all = TRUE), data)
# replace NAs with 0
mergedData[is.na(mergedData)] <- 0 
# indicator if any sip
mergedData$anySIP <- rowSums(mergedData[,2:6])

# add bait positions and gene names one row per gene name
mergedData2 <- merge(mergedData, fraser[,c(1:5)], by="baitID", all.x=T, all.y=F) %>%
  unique(., by="baitID") %>% 
  separate_rows(., "baitName", sep=";") %>% data.table(.)
# some baits don't have associated genes, remove those
# want to keep this for excel file - updated 7/31/21
# mergedData2 <- mergedData2[is.na(baitName)==F]

# get ensemblIDs (GRCh37)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl",
                   host="grch37.ensembl.org", path="/biomart/martservice")
gene_names2 <- getBM(attributes=c('hgnc_symbol', "ensembl_gene_id",
                                  "start_position", "end_position"), 
                     filters = 'hgnc_symbol', 
                     values = mergedData2$baitName,
                     mart = ensembl)

# add to sip data (some genes don't have ensembl IDs)
mergedData3 <- merge(mergedData2, data.table(gene_names2),
                     by.x="baitName", by.y="hgnc_symbol",
                     all.x=T, all.y=F)
colnames(mergedData3)[12] <- "ENSEMBL_GENEID"

# fill in missing ensembleIDs (gains 723/7073)
# /proj/yunligrp/users/ycyang/Gencode/gene_info.gencode.v28lift37.txt
grch37 <- fread("/Users/tmlagler/OneDrive/Lab/HiC/data/Final/gene_info.gencode.v28lift37.txt",
                select=c("Gene_ID", "Gene_NAME", "StartSite", "EndSite"))

# missing genes
missing <- mergedData3[is.na(ENSEMBL_GENEID)]
# find which missing genes have info
info <- grch37[Gene_NAME %in% missing$baitName]
missing2 <- missing[baitName %in% info$Gene_NAME]
missingInfo <- merge(missing2, info, by.x="baitName", by.y="Gene_NAME")
missingInfo$ENSEMBL_GENEID <- NULL;
missingInfo$start_position <- NULL; missingInfo$end_position <- NULL
colnames(missingInfo)[12:14] <- c("ENSEMBL_GENEID", "start_position", "end_position")
# add back to full data
mergedData4 <- rbind(mergedData3[!(baitID %in% missingInfo$baitID &
              baitName %in% missingInfo$baitName)], missingInfo)


# BLUEPRINT expression data (mapps baitIDs to ENSEMBL ids)
# /proj/yunligrp/users/jwen/SIP/Fig.3a_combine_mac_mon/GeneExpressionMatrix_combine_Mac_Mon.txt
expr <- fread("GeneExpressionMatrix_combine_Mac_Mon.txt") # Jia's file
# remove genes without baitIDs
expr <- expr[is.na(BaitID)==F]
# exponentiate to get expression comparable to RPKM
expr2 <- data.table(expr[,1:2], apply(expr[,3:7], 2, exp))

# add expression data to sip data 
# genes without ensembl ids will be missing expression data
expr3 <- merge(mergedData4, expr2,
               by.x=c("baitID", "ENSEMBL_GENEID"),
               by.y=c("BaitID", "ENSEMBL_GENEID"), all.x=T, all.y=F)
colnames(expr3)[15] <- "MacMon"

# write to file. info on all genes, SIP and non-SIP
# this file will now include SIPs without gene names
fwrite(expr3, "allSIPs_geneInfo.txt", sep="\t", col.names=T)

# no longer needed
# write file of all gene positions for LDSC
# positions <- expr3[!is.na(start_position),
#                    c("ENSEMBL_GENEID", "baitChr", "start_position", "end_position")]
# colnames(positions) <- c("GENE", "CHR", "START", "END")
# positions$CHR <- paste0("chr", positions$CHR)
# fwrite(positions, "ENSG_coord_v2.txt", sep="\t", col.names=T)

#----------------------------------------------------------------------------------#
# count shared and cell type specific SIP genes
# write to file
#----------------------------------------------------------------------------------#
#expr3 <- fread("allSIPs_geneInfo.txt")
expr3 <- expr3[anySIP > 0]
# accounts for the same gene at multiple SIP baits
expr3 <- merge(expr3, expr3[, .("anySIP2"=sum(anySIP)/sum(!!anySIP)), by="baitName"],
               by="baitName",all.x=T)

# SIP expression data (cell type specific)
EryExpr <- expr3[SIP.Ery==1 & anySIP2==1, ]
MacMonExpr <- expr3[SIP.MacMon==1 & anySIP2==1, ]
MKExpr <- expr3[SIP.MK==1 & anySIP2==1, ]
nCD4Expr <- expr3[SIP.nCD4==1 & anySIP2==1, ]
NeuExpr <- expr3[SIP.Neu==1 & anySIP2==1, ]

# cell type-group specific genes and shared genes
ErySpecific <- unique(EryExpr, by="baitName")
MacMonSpecific <- unique(MacMonExpr, by="baitName")
MKSpecific <- unique(MKExpr, by="baitName")
nCD4Specific <- unique(nCD4Expr, by="baitName")
NeuSpecific <- unique(NeuExpr, by="baitName")
Shared <- unique(expr3[anySIP==5, ], by="baitName")

# all genes
allGenes <- rbind(ErySpecific, MacMonSpecific, MKSpecific,
                  nCD4Specific, NeuSpecific, Shared)
allGenes$CT <- c(rep("Ery", nrow(ErySpecific)), rep("MacMon", nrow(MacMonSpecific)),
                 rep("MK", nrow(MKSpecific)), rep("nCD4", nrow(nCD4Specific)),
                 rep("Neu", nrow(NeuSpecific)), rep("Shared", nrow(Shared)))

# find shared genes to remove
rmEry <- intersect(ErySpecific$baitName, allGenes[CT != "Ery",baitName])
rmMacMon <- intersect(MacMonSpecific$baitName, allGenes[CT != "MacMon",baitName])
rmMK <- intersect(MKSpecific$baitName, allGenes[CT != "MK",baitName])
rmnCD4 <- intersect(nCD4Specific$baitName, allGenes[CT != "nCD4",baitName])
rmNeu <- intersect(NeuSpecific$baitName, allGenes[CT != "Neu",baitName])

# remove shared genes
ErySpecific <- ErySpecific[!(baitName %in% rmEry)]
MacMonSpecific <- MacMonSpecific[!(baitName %in% rmMacMon)]
MKSpecific <- MKSpecific[!(baitName %in% rmMK)]
nCD4Specific <- nCD4Specific[!(baitName %in% rmnCD4)]
NeuSpecific <- NeuSpecific[!(baitName %in% rmNeu)]

# no longer needed
# write only ensemblIDs to file
# fwrite(list(ErySpecific$ENSEMBL_GENEID[ErySpecific$ENSEMBL_GENEID!=""]),
#        "Ery_specific_genes.txt", col.names=F, sep="\t")
# fwrite(list(MacMonSpecific$ENSEMBL_GENEID[MacMonSpecific$ENSEMBL_GENEID!=""]),
#        "MacMon_specific_genes.txt", col.names=F, sep="\t")
# fwrite(list(MKSpecific$ENSEMBL_GENEID[MKSpecific$ENSEMBL_GENEID!=""]),
#        "MK_specific_genes.txt", col.names=F, sep="\t")
# fwrite(list(nCD4Specific$ENSEMBL_GENEID[nCD4Specific$ENSEMBL_GENEID!=""]),
#        "nCD4_specific_genes.txt", col.names=F, sep="\t")
# fwrite(list(NeuSpecific$ENSEMBL_GENEID[NeuSpecific$ENSEMBL_GENEID!=""]),
#        "Neu_specific_genes.txt", col.names=F, sep="\t")
# fwrite(list(Shared$ENSEMBL_GENEID[Shared$ENSEMBL_GENEID!=""]),
#        "Shared_genes.txt", col.names=F, sep="\t")

# write entire specific gene info to file 
fwrite(ErySpecific[, c(1:3, 10:12, 15:19)], "Ery_specific_genes_allInfo.txt", col.names=T, sep="\t")
fwrite(MacMonSpecific[, c(1:3, 10:12, 15:19)], "MacMon_specific_genes_allInfo.txt", col.names=T, sep="\t")
fwrite(MKSpecific[, c(1:3, 10:12, 15:19)], "MK_specific_genes_allInfo.txt", col.names=T, sep="\t")
fwrite(nCD4Specific[, c(1:3, 10:12, 15:19)], "nCD4_specific_genes_allInfo.txt", col.names=T, sep="\t")
fwrite(NeuSpecific[, c(1:3, 10:12, 15:19)], "Neu_specific_genes_allInfo.txt", col.names=T, sep="\t")
fwrite(Shared[, c(1:3, 10:12, 15:19)], "Shared_genes_allInfo.txt", col.names=T, sep="\t")

#----------------------------------------------------------------------------------#
# write SIP gene info to Excel files for supplemental tables
# add back in SIP score and rank

expr3 <- fread("allSIPs_geneInfo.txt")
all <- expr3[anySIP > 0] # only SIPs
all$start_position <- NULL; all$end_position <- NULL

getCT <- function(CT){
  # select SIPs in cell type of interest  
  allCT <- all[eval(as.name(paste0("SIP.", CT))) > 0 ]
  # colKeep <- c("baitID", "baitName", "baitChr", "baitStart", "baitEnd",
  #             "ENSEMBL_GENEID", "Ery", "MacMon", "MK", "nCD4", "Neu")
  # allCT <- allCT[, ..colKeep]
  # replace missing with NA
  allCT[ENSEMBL_GENEID=="", ]$ENSEMBL_GENEID <- NA
  
  # SIP cumulative interaction score
  sipScore <- fread(paste0("SIP_", CT, ".txt"), drop=3:4)
  # rank highest to lowest 
  setorder(sipScore, -intscore)
  sipScore$rank <- seq(1:nrow(sipScore))
  
  # combine
  sipCT <- merge(allCT, sipScore, by="baitID", all.x=T, all.y=F)
  setorderv(sipCT, CT, order=-1, na.last=T)
  # remove duplicates
  sipCT <- unique(sipCT, by=c("baitID", "baitName"))
  
  # specific/shared
  # note that some genes may correspond to multiple baits and therefore the sum
  # of specific/shared genes in the table will be greater than the unique total
  specific <- fread(paste0(CT, "_specific_genes_allInfo.txt"), select="baitName")
  setDT(sipCT)[, specific := ifelse(baitName %in% specific$baitName, 1, 0)]
  shared <- fread("shared_genes_allInfo.txt", select = "baitName")
  setDT(sipCT)[, shared := ifelse(baitName %in% shared$baitName, 1, 0)]
  
  # reorder and write file
  sipCT <- sipCT[, c("baitID", "baitName", "baitChr", "baitStart", "baitEnd",
                     "intscore", "rank", "specific", "shared",
                     "SIP.Ery", "SIP.MacMon", "SIP.MK", "SIP.nCD4", "SIP.Neu",
                     "ENSEMBL_GENEID", "Ery", "MacMon", "MK", "nCD4", "Neu")]
  setorder(sipCT, rank)
  xlsx::write.xlsx(sipCT, file=paste0("SIPgeneInfo.xlsx"),
                   sheetName=paste0(CT, "_SIPs"), append=T, row.names=F)
}

getCT("Ery")
getCT("MacMon")
getCT("MK")
getCT("nCD4")
getCT("Neu")

#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#






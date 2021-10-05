#----------------------------------------------------------------------------------------------#
# Figure S2: Upset Plot of SIP calls
# requires SIP_*.txt
# outputs upset plot and txt files of specific SIPs
#----------------------------------------------------------------------------------------------#
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(UpSetR)
setwd("/Users/tmlagler/OneDrive/Lab/SIP/figureCode/")
#----------------------------------------------------------------------------------------------#
# read in SIP call files
#----------------------------------------------------------------------------------------------#
# /proj/yunligrp/users/lagler/SIP/data/
files <- list.files(pattern="^SIP.*\\.txt$")
files <- files[-grep("specific", files)] # don't include specific files
data <- lapply(files, fread) 
cells <- tstrsplit(files, "_")[[2]] %>% tstrsplit(., ".txt") %>% unlist()
# change SIP column to by CT, remove unnecessary columns
for(i in 1:5){
  colnames(data[[i]])[4] <- cells[i]
  data[[i]]$intscore <- NULL
  data[[i]]$rank <- NULL
}

# merge all CTs together
mergedData <- Reduce(function(...) merge(..., all = TRUE), data)
# replace NAs with 0
mergedData[is.na(mergedData)] <- 0 


#----------------------------------------------------------------------------------------------#
# make upset plot
#----------------------------------------------------------------------------------------------#
upsetPlot <- upset(mergedData, 
      sets = rev(cells), # reverse order so Ery is on top
      sets.bar.color = "black",
      main.bar.color="grey28",
      sets.x.label = "# of SIPs per Cell Type",
      #c(int. size title, int. size tick labels, set size title,
      # set size tick labels, set names, numbers above bars)
      text.scale = c(1.5, 1.25, 1.25, 1, 1.3, 1.75),
      order.by = "freq", keep.order = T,
      # color shared and specific groups
      queries = list(list(query = intersects, 
                          params = list(cells), color = "green4", active = T),
                     list(query = intersects, 
                          params = list(cells[1]), color = "deepskyblue3", active = T),
                     list(query = intersects, 
                          params = list(cells[2]), color = "deepskyblue3", active = T),
                     list(query = intersects, 
                          params = list(cells[3]), color = "deepskyblue3", active = T),
                     list(query = intersects, 
                          params = list(cells[4]), color = "deepskyblue3", active = T),
                     list(query = intersects, 
                          params = list(cells[5]), color = "deepskyblue3", active = T)))
# save plot
#tiff("../figures/FigureS2_upset.tiff", res=300, height=7, width=10.5, units="in")
png("../figures/FigureS2_upset.png", res=300, height=7, width=10.5, units="in")
upsetPlot
dev.off()

#----------------------------------------------------------------------------------------------#
# save cell type group specific and shared SIP
#----------------------------------------------------------------------------------------------#
# probably a more efficient way to do this...
setDT(mergedData)[,Shared := ifelse(Ery==1 & MacMon==1 & MK==1 & nCD4==1 & Neu==1, 1, 0)]
setDT(mergedData)[,ErySpecific := ifelse(Ery==1 & MacMon==0 & MK==0 & nCD4==0 & Neu==0, 1, 0)]
setDT(mergedData)[,MacMonSpecific := ifelse(Ery==0 & MacMon==1 & MK==0 & nCD4==0 & Neu==0, 1, 0)]
setDT(mergedData)[,MKSpecific := ifelse(Ery==0 & MacMon==0 & MK==1 & nCD4==0 & Neu==0, 1, 0)]
setDT(mergedData)[,nCD4Specific := ifelse(Ery==0 & MacMon==0 & MK==0 & nCD4==1 & Neu==0, 1, 0)]
setDT(mergedData)[,NeuSpecific := ifelse(Ery==0 & MacMon==0 & MK==0 & nCD4==0 & Neu==1, 1, 0)]

# add specific/shared indicator to SIP call
# recall files so they're unmodified
files <- list.files(pattern="^SIP.*\\.txt$")
data <- lapply(files, fread)

# save to file
# /proj/yunligrp/users/lagler/SIP/data/
fwrite(merge(data[[1]], mergedData[, c("baitID", "ErySpecific")],
           by="baitID", all.x = T, all.y = F),
       "SIP_Ery_specific.txt", col.names=T, sep="\t")
fwrite(merge(data[[2]], mergedData[, c("baitID", "MacMonSpecific")],
             by="baitID", all.x = T, all.y = F),
       "SIP_MacMon_specific.txt", col.names=T, sep="\t")
fwrite(merge(data[[3]], mergedData[, c("baitID", "MKSpecific")],
             by="baitID", all.x = T, all.y = F),
       "SIP_MK_specific.txt", col.names=T, sep="\t")
fwrite(merge(data[[4]], mergedData[, c("baitID", "nCD4Specific")],
             by="baitID", all.x = T, all.y = F),
       "SIP_nCD4_specific.txt", col.names=T, sep="\t")
fwrite(merge(data[[5]], mergedData[, c("baitID", "NeuSpecific")],
             by="baitID", all.x = T, all.y = F),
       "SIP_Neu_specific.txt", col.names=T, sep="\t")

#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#

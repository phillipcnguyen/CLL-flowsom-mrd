#Merging FCS files
#Multiple FCS files are sometimes required to achieve sufficient WBC counts
#This script merges files with the same name and edits channel names

library(FlowSOM) 
library(flowCore)

#Set working directory 
setwd("/Users/phil/Documents/FlowSOM project/CLLproject/Duraclone/")

#Loop to merge and rename
dir_split <- "FCS files/Split validation 2/" 
dir_merged <- "FCS files/Merged validation 2/"
split_str <- sub(" .*","",list.files(path=dir_split)) 
split_uniq <- split_str[!duplicated(split_str)]
for (i in split_uniq){
  split <- list.files(path=dir_split, pattern=i)
  agg <- AggregateFlowFrames(paste0(dir_split, split),
                             dataset=2,
                             cTotal = 1000000000000,
                             writeOutput = FALSE,
                             outputFile = paste0(dir_merged, i,".fcs"))
  colnames(agg)[6:15] <- c("CD81 FITC","ROR1 PE","Nil1 ECD","CD79b PerCP-Cy5.5","CD19 PE-Cy7","CD5 APC","Nil2 AF700","CD43 AF750","CD20 Pac Blue","CD45 KO")
  parameters(agg)@data$desc[5:15] <- c("SS-A","CD81","ROR1","Nil1","CD79b","CD19","CD5","Nil2","CD43","CD20","CD45")
  write.FCS(agg,file = paste0(dir_merged, i,".fcs"))
}
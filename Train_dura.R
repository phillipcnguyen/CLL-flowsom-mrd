#Training the Reference SOM for Duraclone samples
#Code adaptated From Quintelier et al, Nature Protocols, 2021

library(FlowSOM) 
library(ggplot2)
library(openCyto) 
library(CytoExploreR) 
library(flowWorkspace)

#Set working directory 
setwd("/Users/phil/Documents/FlowSOM project/CLLproject/Duraclone/")

#Establishing directories
#Put merged (see Merging_files.R) training FCS files into Training/Data/Raw
dir_data <- "Training/Data/"
dir_raw <- "Training/Data/Raw/" 
dir_prepr <- "Training/Data/Preprocessed/" 
dir_QC <- "Training/Data/Preprocessed/QC/" 
dir_RDS <- "Training/RDS/" 
dir_results <- "Training/Results/" 
path_comp <- "Compensation1.csv" 
for (path in c(dir_data, dir_raw, dir_prepr, dir_QC, dir_RDS, dir_results)){
  dir.create(path)}

#Reading in reference sample to generate transformation list
#Markers of interest searches for the text in the channel names
#Depending on FCS file output, may need to rename channels. See Merging_files.R
reference_sample <- read.FCS("Training/Data/Raw/12882429.fcs")
reference_marker <- "ROR1 PE" 
markers_of_interest <- c("SS-A","CD81","ROR1 PE","Nil1","CD79b","CD19","CD5","Nil2","CD43","CD20","CD45")
channels_of_interest <- GetChannels(object = reference_sample, markers = markers_of_interest, exact = FALSE)
compensation_matrix <- read.csv(path_comp, check.names = FALSE, row.names = 1)
colnames(compensation_matrix) <- channels_of_interest[2:11]
rownames(compensation_matrix) <- channels_of_interest[2:11]

#Compute transformation list
ff_m <- PeacoQC::RemoveMargins(reference_sample, channels_of_interest)
ff_c <- flowCore::compensate(ff_m, compensation_matrix)
translist <- estimateLogicle(ff_c, colnames(compensation_matrix), m=5) 
ff_t <- flowCore::transform(ff_c, translist)
q5_goal <- quantile(exprs(ff_t)[,reference_marker], 0.05) 
q95_goal <- quantile(exprs(ff_t)[,reference_marker], 0.95) 
q5_SSCA <- quantile(exprs(ff_t)[,"SS-A"], 0.05)
q95_SSCA <- quantile(exprs(ff_t)[,"SS-A"], 0.95)
SSCA_a <- (q95_goal - q5_goal) / (q95_SSCA - q5_SSCA)
SSCA_b <- q5_goal - q5_SSCA * (q95_goal - q5_goal) / (q95_SSCA - q5_SSCA) 
translist <- c(translist, transformList("SS-A", flowCore::linearTransform(a = SSCA_a, b = SSCA_b)))

#Generating preprocessed files 
files <- list.files(path = dir_raw)
for(i in files){
  ff <- read.FCS(paste0(dir_raw,i), truncate_max_range = FALSE) 
  ff_c <- flowCore::compensate(ff, compensation_matrix) 
  ff_t <- flowCore::transform(ff_c, translist) 
  ff_s <- PeacoQC::RemoveDoublets(ff_t, channel1 = "FS-A", channel2 = "FS-H") 
  gate <-cyto_gate_draw(ff_s, parent = "root", alias = "Lymphocytes", channels = c("SS-A","CD45 KO"), xlim=c(0.25,3))
  gate_matrix <- gate[["Lymphocytes"]]@boundaries
  lymph_gate <- flowCore::polygonGate(.gate = gate_matrix)
  lymph_filt <- filter(ff_s, lymph_gate)
  ff_lymph <- ff_s[lymph_filt@subSet,]
  write.FCS(ff_lymph,file = paste0(dir_prepr, i))
}

#Creating aggregate file from preprocessed files
agg <- AggregateFlowFrames(paste0(dir_prepr, files), 
                           cTotal = 10000000000,
                           writeOutput = TRUE,
                           outputFile = paste0(dir_prepr, "aggregate.fcs"))
agg <- read.FCS(paste0(dir_prepr, "aggregate.fcs")) 

#Hierarchical clustering to generate Reference SOM
#level 1 clustering to identify B cells
fsom1 <- FlowSOM(input = agg,
                 scale = FALSE,
                 colsToUse = c("SS-A", "CD19 PE-Cy7"),
                 seed = 1,
                 nClus = 5,
                 xdim = 1, ydim = 10)
PlotStars(fsom = fsom1,  
          backgroundValues = fsom1$metaclustering,
          equalNodeSize = FALSE,
          maxNodeSize = 1,
          list_insteadof_ggarrange = FALSE)
ggsave(paste0(dir_results, "fsom1tree.pdf"), height = 8.5, width = 11)

#Hierarachical clustering to generate Reference SOM
#Filtering B cells to pass to level 2
#In this example, the B cells are in metacluster 1,2. 
fsom_tmp <- NewData(fsom=fsom1, input=agg) 
agg_bcells <- agg[GetMetaclusters(fsom_tmp) %in% c(1,2),] #Adapt numbers to metaclusters containing B cells
write.FCS(agg_bcells,file = paste0(dir_prepr, "agg_bcells.fcs"))

#Hierarchical clustering to generate Reference SOM: 
#Level 2 clustering
fsom2 <- FlowSOM(input = agg_bcells,
                 scale = FALSE,
                 colsToUse = markers_of_interest[c(2,3,5,6,7,9,10,11)],
                 seed = 1,
                 nClus = 10,
                 xdim = 10, ydim = 10)
PlotStars(fsom = fsom2,  
          backgroundValues = fsom2$metaclustering,
          equalNodeSize = FALSE,
          maxNodeSize = 1.5,
          list_insteadof_ggarrange = FALSE)
ggsave(paste0(dir_results, "fsom2tree.pdf"), height = 8.5, width = 11)

#Save fsom1 and fsom2 for use in validation sets
saveRDS(fsom1, paste(dir_RDS,"fsom1dura.rds"))
saveRDS(fsom2, paste(dir_RDS,"fsom2dura.rds"))

#Overview of flowsom results
FlowSOMmary(fsom = fsom1, plotFile = paste0(dir_results, "fsom1dura_summ.pdf"))
FlowSOMmary(fsom = fsom2, plotFile = paste0(dir_results, "fsom2dura_summ.pdf"))

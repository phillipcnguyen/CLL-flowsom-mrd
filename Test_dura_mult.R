#Determining MRD as a loop for multiple files

library(FlowSOM) 
library(ggplot2)
library(openCyto) #Uses Raphael Gottardos Rcytosuite for gating/compensation
library(CytoExploreR) 
library(flowWorkspace)

#Set working directory 
setwd("/Users/phil/Documents/FlowSOM project/CLLproject/Duraclone/")

#Load the Reference SOM
#These fsom maps made in the Train_dura.R script
fsom1 <- readRDS("Training/RDS/ fsom1dura.rds") 
fsom2 <- readRDS("Training/RDS/ fsom2dura.rds")  

#Establishing other directories.
#Put merged validation FCS files into Validation/Data/Raw
dir_data <- "Validation/Data/"
dir_raw <- "Validation/Data/Raw/" 
dir_prepr <- "Validation/Data/Preprocessed/" 
dir_QC <- "Validation/Data/Preprocessed/QC/" 
dir_RDS <- "Validation/RDS/" 
dir_results <- "Validation/Results/" 
path_comp <- "Compensation1.csv" 
for (path in c(dir_data, dir_raw, dir_prepr, dir_QC, dir_RDS, dir_results)){
  dir.create(path)}
files <- list.files(path = dir_raw)

#Reading in Duraclone FCS files, channels, and compensation matrix
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

# Generating preprocessed files from raw with QC checks
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

#Mapping test files to Reference SOM
validation <- vector(mode= "list", length = length(files))
validation <- setNames(validation,files)
for(i in files){
  ff_prepr <- read.FCS(paste0(dir_prepr, i))
  fsom_lymphs <- NewData(fsom=fsom1, input=ff_prepr) 
  clustering <- GetMetaclusters(fsom_lymphs)
  ff_bcells <- ff_prepr[clustering %in% c(1,2),]
  fsom_bcells <- NewData(fsom = fsom2,input = ff_bcells)
  validation[[i]] <- fsom_bcells
}

#Generating Results Table
CLL_metacluster <- 4
results <- matrix(NA, nrow = length(validation), 
                  ncol=2, 
                  byrow=FALSE, 
                  dimnames=list(NULL,c("File", "Events")))
results[,1] <- files
for (i in 1:length(files)){ 
  fsom_temp <- validation[[i]]
  counts <- GetCounts(fsom=fsom_temp, level = "metaclusters")[CLL_metacluster]
  counts[is.na(counts)] <-0 
  results[i,2] <- counts
}
print(results)
write.csv(results, file = paste0(dir_results,"results.csv"))
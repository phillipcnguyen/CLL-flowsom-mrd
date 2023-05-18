#Used for determining MRD in the validation cohort(s)
#Script below was adapted for each of the different cohorts (ie typical or atypical)
#Code adapted From Quintelier et al, Nature Protocols, 2021

library(FlowSOM) 
library(ggplot2)
library(openCyto)
library(CytoExploreR)
library(dplyr)

# Reading in the reference map
fsom <- readRDS("Documents/FlowSOM project/CLLproject/Training/RDS/ fsom.RDS") #Loading reference flowsom
fsom_level2 <- readRDS("Documents/FlowSOM project/CLLproject/Training/RDS/ fsom_level2.RDS") #Loading reference flowsom

#Define variables
file_pattern <- ".fcs" 
reference_file <- read.FCS("Documents/FlowSOM project/CLLproject/Validation/Data/Raw/12349141_VINE_CLLMRDPRE_COMBINATION.fcs")
reference_marker <- "PE-A" 
markers_of_interest <- c("SSC-A", "CD81FITC", "CD79b PE", "CD22 PE-Cy5-5", "CD19 PE-Cy7", "CD43 APC", "CD20 APC-H7", "CD5 BV-421", "CD3 V510")

#Establishing directories
dir_prepr <- "Documents/FlowSOM project/CLLproject/Validation/Data/Preprocessed/" #where the preprocessed data will be stored
dir_QC <- "Documents/FlowSOM project/CLLproject/Validation/Data/Preprocessed/QC/" #where the data QC results will be stored
dir_RDS <- "Documents/FlowSOM project/CLLproject/Validation/RDS/" #where the R objects will be stored
dir_results <- "Documents/FlowSOM project/CLLproject/Validation/Results/" #where the results will be stored
dir_raw <- "Documents/FlowSOM project/CLLproject/Validation/Data/Raw/" #where the raw data is located
path_comp <- "Documents/FlowSOM project/CLLproject/Compensation.csv"
for (path in c(dir_prepr, dir_QC, dir_RDS, dir_results)){
  dir.create(path)}

#Additional information for preprocessing
files <- list.files(path = dir_raw, pattern = file_pattern)
channels_of_interest <- GetChannels(object = reference_file, markers = markers_of_interest, exact = FALSE)
compensation_matrix <- read.csv(path_comp, check.names = FALSE, row.names = 1)
colnames(compensation_matrix) <- sub(" :: .*", "", colnames(compensation_matrix))

#Compute transformation list
ff_m <- PeacoQC::RemoveMargins(reference_file, channels_of_interest)
ff_c <- flowCore::compensate(ff_m, compensation_matrix)
translist <- estimateLogicle(ff_c, colnames(compensation_matrix)) 
ff_t <- flowCore::transform(ff_c, translist)
q5_goal <- quantile(exprs(ff_t)[,reference_marker], 0.05) 
q95_goal <- quantile(exprs(ff_t)[,reference_marker], 0.95) 
q5_SSCA <- quantile(exprs(ff_t)[,"SSC-A"], 0.05)
q95_SSCA <- quantile(exprs(ff_t)[,"SSC-A"], 0.95)
SSCA_a <- (q95_goal - q5_goal) / (q95_SSCA - q5_SSCA)
SSCA_b <- q5_goal - q5_SSCA * (q95_goal - q5_goal) / (q95_SSCA - q5_SSCA) 
translist <- c(translist, transformList("SSC-A", flowCore::linearTransform(a = SSCA_a, b = SSCA_b)))

#Preprocessing pipeline for all files using a for loop
for (file in files){
  ff <- read.FCS(paste0(dir_raw, file), truncate_max_range = FALSE)
  ff_m <- PeacoQC::RemoveMargins(ff, channels_of_interest)
  ff_c <- flowCore::compensate(ff_m, compensation_matrix)
  ff_t <- flowCore::transform(ff_c, translist)
  ff_s <- PeacoQC::RemoveDoublets(ff_t)
  gate <-cyto_gate_draw(ff_s, parent = "root", alias = "Lymphocytes", channels = c("FSC-A", "SSC-A"), ylim=c(0,2))
  gate_matrix <- gate[["Lymphocytes"]]@boundaries
  lymph_gate <- flowCore::polygonGate(filterId = "Lymphocytes", .gate = gate_matrix)
  lymph_filter <- filter(ff_s, lymph_gate)
  ff_l <- ff_s[lymph_filter@subSet,]
  PQC <- PeacoQC::PeacoQC(ff = ff_l,
                          channels = channels_of_interest,
                          plot = TRUE, save_fcs = FALSE,
                          output_directory = dir_QC)
  write.FCS(PQC$FinalFF,
            file = paste0(dir_prepr, file))
  to_plot <- list(list(ff_pre = ff,
                       ff_post = ff_m,
                       title = "Removed margin events",
                       channel_x = "PE-Cy7-A",
                       channel_y = "APC-A"),
                  list(ff_pre = ff_t,
                       ff_post = ff_s,
                       title = "Removed doublets",
                       channel_x = "FSC-A",
                       channel_y = "FSC-H"),
                  list(ff_pre = ff_s,
                       ff_post = ff_l,
                       title = "Removed debris and dead cells",
                       channel_x = "FSC-A",
                       channel_y = "SSC-A"),
                  list(ff_pre = ff_l,
                       ff_post = PQC$FinalFF,
                       title = "Removed low quality events",
                       channel_x = "Time",
                       channel_y = "PE-Cy7-A"))
  plot_list <- list()
  for (plot in to_plot) {
    plot_list[[length(plot_list) + 1]] <- filter_plot(ff_pre = plot$ff_pre,
                                                      ff_post = plot$ff_post,
                                                      title = plot$title,
                                                      channel_x = plot$channel_x,
                                                      channel_y = plot$channel_y)
  }
  png(paste0(dir_QC, sub("fcs", "png", file)), width = 1920)
  print(ggpubr::ggarrange(plotlist = plot_list, nrow = 1))
  dev.off()
}

#Generate a vector with multiple mapped FCS files
validation <- vector(mode= "list", length = length(files))
validation <- setNames(validation,files)
for(i in 1:length(files)){
  file <- files[i]
  ff_prepr <- read.FCS(paste0(dir_prepr, file))
  ff_level1 <- NewData(fsom=fsom, input=ff_prepr) #Not convinced I need this line of code
  clustering <- GetMetaclusters(ff_level1)
  ff_prepr_level1 <- ff_prepr[clustering %in% c(1,3),]
  ff_level2 <- NewData(fsom = fsom_level2,input = ff_prepr_level1)
  validation[[i]] <- ff_level2
}

# Generating results matrix from the validation vector
# Define which clusters or metaclusters are disese
CLL_typical_clusters <- c(193,208,209,224,225,203,217, 218,219,223,207,222)
CLL_typical_metaclusters <- c(28, 29, 30)
results <- matrix(NA, nrow = length(validation), 
                  ncol=3, 
                  byrow=FALSE, 
                  dimnames=list(NULL,c("File", "Events", "Percent_bcells")))
results[,1] <- files
for (i in 1:length(files)){ 
  fsom_temp <- validation[[i]]
  counts <- GetCounts(fsom=fsom_temp, level = "clusters")[CLL_typical_clusters] 
  percentages <- GetPercentages(fsom=fsom_temp, level = "clusters")[CLL_typical_clusters] 
  counts[is.na(counts)] <-0 
  percentages[is.na(percentages)] <-0
  counts <- sum(counts)
  percentages <- sum(percentages)
  results[i,2] <- counts
  results[i,3] <- percentages
}
print(results)
write.csv(results, file = paste0(dir_results,"results.csv"))

#Generating an FCS file with FlowSOM MST for interrogation in Kaluza
#Adapated from https://github.com/bc2zb/cyttools/blob/master/FlowSOM.R#L81-L120

outlier_IDnumber <- c(5,35) #Numbers correspond to files where Kaluza interrogation desired
for(i in outlier_IDnumber){
test <- validation[[i]]
file.interest <- files[i] 
ResultsTable <- as.data.frame(test$data)
fileLabels <- test$metaData %>% unlist() %>% matrix(ncol = 2, byrow = T)
fileNames <- vector(length = nrow(ResultsTable))
fileNames <- c(rep(files[i], nrow(ResultsTable)))
ResultsTable$FileNames <- fileNames
ResultsTable$Mapping <- test$map$mapping[,1]
ResultsTable$DistToNode <- test$map$mapping[,2]
ResultsTable <- ResultsTable %>% left_join(test$MST$l %>%
                                             as.data.frame() %>%
                                             setNames(c("Flowsom_x", "Flowsom_y")) %>%
                                             tibble::rownames_to_column("Mapping") %>% ##added tibble to get rownames_to_column recognised
                                             mutate(Mapping = as.numeric(Mapping)))
Prepr.FCS <- read.FCS(paste0(dir_prepr,file.interest))
Raw.FCS <- read.FCS(paste0(dir_raw,file.interest))
Level1.FSOM <- NewData(fsom = fsom, input = Prepr.FCS) #Maps Prepre.FCS onto the FlowSOM 1
Level1.Metaclusters <- GetMetaclusters(Level1.FSOM) # Generates a factor with cells and the associated metacluster numbers of FlowSOM 1
Filtered.Prepr.FCS <- Prepr.FCS[Level1.Metaclusters %in% c(1,3),] #Flowframe with only selected metaclusters
Level2.FSOM <- NewData(fsom = fsom_level2,input = Filtered.Prepr.FCS) #Maps this filtered flowframe to FlowSOM 2
Level2.Metaclusters <- GetMetaclusters(Level2.FSOM) # Generates a factor with cells and the associated metacluster numbers of FlowSOM 2

clusterData <- ResultsTable %>%
  select(Mapping, DistToNode, Flowsom_x, Flowsom_y)
clustering_raw <- matrix(0, nrow(exprs(Raw.FCS)), 4)
colnames(clustering_raw) <- c("Mapping", "DistToNode","Flowsom_x","Flowsom_y") 
clustering_raw[c(exprs(Filtered.Prepr.FCS)[,"Original_ID"]),] <- as.matrix(clusterData)
clusterFCS <- flowCore::fr_append_cols(Raw.FCS, clustering_raw)
write.FCS(clusterFCS, paste0(paste(dir_results,"Outliers/", sep=""), gsub(".fcs", "",file.interest), "_flowsom_MST.fcs"))
}


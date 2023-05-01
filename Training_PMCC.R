#Training the Reference SOM for PMCC diagnostic samples
#Code adapted From Quintelier et al, Nature Protocols, 2021

library(FlowSOM) 
library(ggplot2)
library(openCyto) 
library(CytoExploreR) 
library(flowWorkspace)

#Define variables
file_pattern <- ".fcs" 
reference_file <- read.FCS("Documents/FlowSOM project/CLLproject/Training/Data/Raw/0_12349141_VINE_CLLMRDPRE_COMBINATION.fcs")
reference_marker <- "PE-A" 
markers_of_interest <- c("SSC-A", "CD81FITC", "CD79b PE", "CD22 PE-Cy5-5", "CD19 PE-Cy7", "CD43 APC", "CD20 APC-H7", "CD5 BV-421", "CD3 V510")

#Establishing directories
dir_prepr <- "Documents/FlowSOM project/CLLproject/Training/Data/Preprocessed/" #where the preprocessed data will be stored
dir_QC <- "Documents/FlowSOM project/CLLproject/Training/Data/Preprocessed/QC/" #where the data QC results will be stored
dir_RDS <- "Documents/FlowSOM project/CLLproject/Training/RDS/" #where the R objects will be stored
dir_results <- "Documents/FlowSOM project/CLLproject/Training/Results/" #where the results will be stored
dir_raw <- "Documents/FlowSOM project/CLLproject/Training/Data/Raw/" #where the raw data is located
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

#Preprocessing loop for all files
for (file in files[1]){
  ff <- read.FCS(paste0(dir_raw, file), truncate_max_range = FALSE)
  ff_m <- PeacoQC::RemoveMargins(ff, channels_of_interest)
  ff_c <- flowCore::compensate(ff_m, compensation_matrix)
  ff_t <- flowCore::transform(ff_c, translist)
  ff_s <- PeacoQC::RemoveDoublets(ff_t)
  # Following lines of code allows gating of lymphocytes for each FCS file
  gate <-cyto_gate_draw(ff_s, parent = "root", alias = "Lymphocytes", channels = c("FSC-A", "SSC-A"), ylim=c(0,2))
  gate_matrix <- gate[["Lymphocytes"]]@boundaries
  lymph_gate <- flowCore::polygonGate(filterId = "Lymphocytes", .gate = gate_matrix)
  lymph_filter <- filter(ff_s, lymph_gate)
  ff_l <- ff_s[lymph_filter@subSet,] #New flow frame of lymphocytes
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

#Creating aggregate file.
n <- 50000000 #High number to capture all events
agg <- AggregateFlowFrames(paste0(dir_prepr, files), 
                           cTotal = n,
                           writeOutput = TRUE,
                           outputFile = paste0(dir_prepr, "aggregate.fcs"))

#Generating the reference SOM
SOM_x <- 10
SOM_y <- 10
n_meta <- 10
seed <- 1
scaling <- FALSE
fsom <- FlowSOM(input = agg,
                scale = scaling,
                colsToUse = c("CD3 V510", "CD19 PE-Cy7"),
                seed = seed,
                nClus = n_meta,
                xdim = SOM_x, ydim = SOM_y)
PlotStars(fsom = fsom,  
          backgroundValues = fsom$metaclustering,
          equalNodeSize = TRUE,
          list_insteadof_ggarrange = FALSE)
ggsave(paste0(dir_results, "fsom_level1_tree.pdf"), height = 8.5, width = 11)
fsom_tmp <- NewData(fsom=fsom, input=agg) 
clustering <- GetMetaclusters(fsom_tmp)
ff_tmp <- agg[clustering %in% c(1,3),] 
SOM_x_level2 <- 15
SOM_y_level2 <- 15
n_meta_level2 <- 30
fsom_level2 <- FlowSOM(input = ff_tmp,
                       scale = scaling,
                       colsToUse = markers_of_interest[2:9],
                       seed = seed,
                       nClus = n_meta_level2,
                       xdim = SOM_x_level2, ydim = SOM_y_level2)
PlotStars(fsom = fsom_level2,  
          backgroundValues = fsom_level2$metaclustering,
          equalNodeSize = TRUE,
          maxNodeSize = 1.5,
          list_insteadof_ggarrange = FALSE)
ggsave(paste0(dir_results, "fsom_level2_tree.pdf"), height = 8.5, width = 11)

#Save fsom1 and fsom2 for use in validation sets
saveRDS(fsom, paste(dir_RDS,"fsom.rds"))
saveRDS(fsom_level2, paste(dir_RDS,"fsom_level2.rds"))

#Overview of flowsom results
FlowSOMmary(fsom = fsom, plotFile = paste0(dir_results, "fsom_summary.pdf"))
FlowSOMmary(fsom = fsom_level2, plotFile = paste0(dir_results, "fsom_level2_summary.pdf"))



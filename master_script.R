# Master script for lncRNA analysis
# Load external scripts
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("cluster_deseq.R")
source("get_neighbours.R")
source("analyse_neighbours.R")
library(tidyverse)
library(data.table)
library(hash)
library(pals)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(tibble)
# <- list.files(path= "/home/harry/Documents/lncRNA_21_10/counts/split_files", pattern="*.csv", full.names=TRUE, recursive=FALSE)


prefix="itra"

count_folder<- "/home/harry/Documents/lncRNA_21_10/counts/itra/"
split_file <- paste("/home/harry/Documents/lncRNA_21_10/counts/split_files/",prefix,"_split.csv", sep = "", collapse="")
gtf_file <- "/home/harry/Documents/lncRNA_21_10/annotation/lncRNA_v11.4_merged.gtf"
image_folder <- "/home/harry/Documents/lncRNA_21_10/images/"

#### Count presence in GTF file ####
# Extract transcript lines
gtf <- readLines(gtf_file)
transcript_lines <- grep('\ttranscript\t', gtf)
transcript_lines <- gtf[transcript_lines]
# Count all lncRNA
all_lncRNA <- grep("MSTRG", transcript_lines)
# Count K 
k_id <- grep("class_code \"k\"", transcript_lines)
# Count intergenic
intergenic_id <- grep("class_code \"u\"", transcript_lines)
intergenic_names <- c()
for(i in intergenic_id){
  intergenic <- str_match(transcript_lines[i], "transcript_id \"(.*?)\"")[1,2]
  intergenic_names <- c(intergenic_names, intergenic)
}
# Count antisense
antisense_id <- grep("class_code \"x\"", transcript_lines)
antisense_names <- c()
# GGenerate a false DF for s/as pairs
sense_transcripts <- c()
sense_ids <- c()
false_distance <- c()
for(i in antisense_id){
  antisense <- str_match(transcript_lines[i], "transcript_id \"(.*?)\"")[1,2]
  antisense_names <- c(antisense_names, antisense)
  sense_transcript <- str_match(transcript_lines[i], "cmp_ref \"(.*?)\"")[1,2]
  sense_transcripts <- c(sense_transcripts, sense_transcript)
  sense_id <- str_match(transcript_lines[i], "cmp_ref_gene \"(.*?)\"")[1,2]  
  sense_ids <- c(sense_ids, sense_id)
  false_distance <- c(false_distance, 0)
 
}
sas_out <- do.call(rbind, Map(data.frame, lncRNA=antisense_names, PCG_id=sense_ids, PCG_transcript=sense_transcripts, distance=false_distance))

# Get names for these
lncRNA_lines <- transcript_lines[all_lncRNA]
PCG_lines <- transcript_lines[grep("MSTRG", transcript_lines, invert = TRUE)]

# Generate logfold change
source("repo_deseq.R")
deseq_results <- deseq_func(count_directory = count_folder, split_file_path = split_file,
                        p_val = 0.05 , min_base_mean = 30 , min_LFC = 1, signif_base_mean = 5000)




setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
logchange <- deseq_results[[1]]
signif_bmean <- deseq_results[[2]]
raw_bmean <- deseq_results[[3]]
colnames(raw_bmean) <- colnames(logchange)
raw_counts <- deseq_results[[4]]

write.csv(logchange,paste(c(prefix,"_lfc.csv"),sep="", collapse=""), row.names = T)
# Only generate heatmaps for itra data
if(prefix == "itra"){
  # Generate clusters and heatmaps
  ### GENERATE HEATMAP - Harry Chown ###
  # Convert data into a matrix which can be used as input for heatmap
  set.seed(2)
  gene_names <- logchange[,1]
  clust_num <- 30
  out<- pheatmap(logchange, kmeans_k = clust_num, scale = "row", 
                 cluster_cols = F)
  out2<- pheatmap(logchange,
                  show_rownames=F, cluster_cols=F, cluster_rows=T, scale="row",
                  cex=1, clustering_distance_rows="euclidean", cex=1,
                  clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE)
  cluster_labels <- out$kmeans$cluster
  
  tiff(file = "/home/harry/Documents/lncRNA_21_10/images/log_cluster_heatmap.tiff",
       unit = "in",
       res=600,
       height = 7,
       width = 4.5)
  
  print(out)
  
  dev.off()
  
  tiff(file = "/home/harry/Documents/lncRNA_21_10/images/log_hierarchical_heatmap.tiff",
       unit = "in",
       res=600,
       height = 7,
       width = 4.5)
  
  print(out2)
  

  
  dev.off()
  logchange_w_cluster <-  add_column(logchange, cluster = cluster_labels, .before = 1)
  write.csv(logchange_w_cluster, "/home/harry/Documents/lncRNA_21_10/clustering/lfc_w_cluster.csv", row.names = T)
  
}








# Identify and annotate lncRNA with neighbouring genes
source("get_neighbours.R")
lncRNA_pcg_neighbours <- get_neighbours(gtf, lncRNA_lines, PCG_lines, 5000, 10)



antisense_pcg_neighbours <- get_neighbours(gtf, transcript_lines[antisense_id], PCG_lines, 5000, 10)
intergenic_pcg_neighbours <- get_neighbours(gtf, transcript_lines[intergenic_id], PCG_lines, 5000, 10)
k_pcg_neighbours <- get_neighbours(gtf, transcript_lines[k_id], PCG_lines, 5000, 10)
# Identify which neighbour pairs are correlated
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("analyse_neighbours.R")
a_n_stats <- analyse_neighbours("antisense_neighbours", "/home/harry/Documents/lncRNA_21_10/neighbours/", antisense_pcg_neighbours, logchange, raw_counts, bmean_df=raw_bmean)
i_n_stats <- analyse_neighbours("intergenic_neighbours", "/home/harry/Documents/lncRNA_21_10/neighbours/", intergenic_pcg_neighbours, logchange, raw_counts, bmean_df=raw_bmean)
k_n_stats <- analyse_neighbours("k_neighbours", "/home/harry/Documents/lncRNA_21_10/neighbours/", k_pcg_neighbours, logchange, raw_counts, bmean_df=raw_bmean)

DE_a_n <- a_n_stats[[1]]
# Identify points with the most extreme values

#p<-ggplot(DE_A_N, aes(x=r_val.cor, y=lncRNA_r_val.cor)) +
#  geom_point() +
#  lims(x=c(-1,1),y=c(-1,1)) +
#  theme_minimal() +
#  coord_fixed() +  
#  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) 
#p


# Identify which sense/antisense pairs are correlated
# Iterate through the pairs and extract the corresponding info
#lncRNA, PCG_id, PCG_transcript, distance (0), ordered_desc


sas_df <- annotate_pcg(sas_out)
sas_stats <- analyse_neighbours("sas", "/home/harry/Documents/lncRNA_21_10/s_as_pairs/", sas_df, logchange, raw_counts, bmean_df=raw_bmean)



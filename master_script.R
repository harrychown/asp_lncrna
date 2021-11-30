# Master script for lncRNA analysis
# Load external scripts
library(tidyverse)
library(data.table)
library(hash)
library(pals)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(tibble)
library("biomaRt")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#source("cluster_deseq.R")
source("get_neighbours.R")
source("analyse_neighbours.R")

### SETUP INPUT/OUTPUT ###
configuration_file <- readLines("analysis.config")
inputdir <- configuration_file[1]
outputdir <- configuration_file[2]
# Save split and count files/folders
split_dir <- paste(c(inputdir, "/split_files"), sep="", collapse="")
all_split_files <- list.files(path=split_dir)
count_dir <- paste(c(inputdir, "/counts"), sep="", collapse="")
gtf_file <- paste(c(inputdir, "/lncRNA_v11.4_merged.gtf"), sep="", collapse="")

### EXTRACT LNCRNA GROUPS FROM GTF ###
# Extract transcript lines
gtf <- readLines(gtf_file)
raw_transcripts <- grep('\ttranscript\t', gtf)
raw_transcripts <- gtf[raw_transcripts]

# Extract all PCG
raw_pcg <- grep("MSTRG", raw_transcripts, invert = T)
raw_pcg <- raw_transcripts[raw_pcg]
pcg_names <- c()
for(i in raw_pcg){
  name <- str_match(i, "gene_id \"(.*?)\"")[1,2]
  pcg_names <- c(pcg_names, name)
}

# Extract all lncRNA
raw_lncrna <- grep("MSTRG", raw_transcripts)
raw_lncrna <- raw_transcripts[raw_lncrna]


# Extract intergenic
raw_intergenic <- grep("class_code \"u\"", raw_lncrna)
raw_intergenic <- raw_lncrna[raw_intergenic]
intergenic_names <- c()
for(i in raw_intergenic){
  name <- str_match(i, "transcript_id \"(.*?)\"")[1,2]
  intergenic_names <- c(intergenic_names, name)
}

# Extract antisense and sense partner
raw_antisense <- grep("class_code \"x\"", raw_lncrna)
raw_antisense <- raw_lncrna[raw_antisense]
antisense_names <- c()
sense_names <- c()
sense_ids <- c()
for(i in raw_antisense){
  a_name <- str_match(i, "transcript_id \"(.*?)\"")[1,2]
  antisense_names <- c(antisense_names, a_name)
  s_name <- str_match(i, "cmp_ref \"(.*?)\"")[1,2]
  sense_names <- c(sense_names, s_name)
  sense_id <- str_match(i, "cmp_ref_gene \"(.*?)\"")[1,2]
  sense_ids <- c(sense_ids, sense_id)
}

# Extract "k"
raw_k <- grep("class_code \"k\"", raw_lncrna)
raw_k <- raw_lncrna[raw_k]
k_names <- c()
for(i in raw_k){
  name <- str_match(i, "transcript_id \"(.*?)\"")[1,2]
  k_names <- c(k_names, name)
}


### OBTAIN DESCRIPTIONS FOR ALL PCGs ###
# Biomart for ensembl fungi
ensembl_fungi <- useMart(host="https://fungi.ensembl.org", 
                         biomart="fungi_mart", 
                         port = 443, dataset = "afumigatusa1163_eg_gene")
# Obtain the attributes of the sense genes
gene_attr <- getBM(attributes=c('ensembl_gene_id', 'description'),
                    filters = "ensembl_gene_id",
                    values = pcg_names, 
                    mart = ensembl_fungi) 
gene_attr <- apply(gene_attr, 2, function(x) gsub("^$|^ $", NA, x))
rownames(gene_attr) <- gene_attr[,1]
gene_attr <- gene_attr[,2]




### IDENTIFY NEIGHBOURING GENES ###
# Use Bedtools closest to identify whether there are lncRNA directly next to PCGs
lncrna_gtf <- raw_lncrna
pcg_gtf <- raw_pcg
# Number of nearest neighbours = 1
num_neighbours  <-  1
# Convert GTF to BED
gtf2bed <- function(gtf_in){
  bed_tmp <- tempfile(pattern = "bed_tmp")
  bed_cmd <- paste(c("gtf2bed | sort -k1,1 -k2,2n > ",bed_tmp), sep = "", collapse = "")
  system(bed_cmd, input = gtf_in, intern = FALSE)
  return(bed_tmp)
}
# Convert lncRNA and PCG GTFs to BED format
lncrna_bed <- gtf2bed(raw_lncrna)
pcg_bed <- gtf2bed(raw_pcg)
# Run closest to find the closest transcript to a PCG
# Perform downstream and upstream analysis so that we can identify what flanks each transcript
closest_downstream <- paste(c("bedtools closest -io -iu -k ", as.character(num_neighbours)," -t all -mdb all -D a -a ", pcg_bed, " -b ", lncrna_bed, " ", pcg_bed), 
                     sep = "", collapse ="")
closest_upstream <- paste(c("bedtools closest -io -id -k ", as.character(num_neighbours)," -t all -mdb all -D a -a ", pcg_bed, " -b ", lncrna_bed, " ", pcg_bed), 
                            sep = "", collapse ="")
downstream_results <- system(closest_downstream, intern = TRUE)
upstream_results <- system(closest_upstream, intern = TRUE)
# Combine upstream and downstream results
closest_results <- c(downstream_results, upstream_results)
# Read closest results and make it human-readable
pcg_transcripts <- c()
pcg_genes <- c()
transcript_neighbours <- c()
distances <- c()
sqrd_distances <- c()
for(r in closest_results){
  split <- strsplit(r, "\t")
  distance <-  as.integer(split[[1]][22])
  distances <- c(distances, distance)
  sqrd_distances <- c(sqrd_distances, distance ** 2)
  pcg_info <- split[[1]][10]
  pcg_transcript <- strsplit(pcg_info, "\"")[[1]][2]
  pcg_transcripts <- c(pcg_transcripts, pcg_transcript)
  pcg_gene <- strsplit(pcg_info, "\"")[[1]][4]
  pcg_genes <- c(pcg_genes, pcg_gene)
  neighbour_info <-  split[[1]][21]
  transcript_neighbour <- strsplit(neighbour_info, "\"")[[1]][2]
  transcript_neighbours <- c(transcript_neighbours, transcript_neighbour)
}
neighbour_out <- do.call(rbind, Map(data.frame, pcg_transcript = pcg_transcripts,
                                    pcg_gene = pcg_genes,
                                    transcript_neighbour = transcript_neighbours,
                                    distance = distances,
                                    sqrd_distance = sqrd_distances))

# Filter neighbour output to analyse only transcripts within 5kb of each other
neighbour_size_filter <- neighbour_out[neighbour_out$sqrd_distance<5000**2,]
# Filter to only include lncRNA neighbours
neighbour_lncrna_filter <- neighbour_size_filter[grep("MSTRG", neighbour_size_filter$transcript_neighbour),]
# Identify which lncRNA are intergenic/antisense/k and add that as a column
lncrna_type <- rep("k", nrow(neighbour_lncrna_filter))
lncrna_type[which(neighbour_lncrna_filter$transcript_neighbour %in% intergenic_names)] <- "intergenic"
lncrna_type[which(neighbour_lncrna_filter$transcript_neighbour %in% antisense_names)] <- "antisense"
neighbour_lncrna <- cbind(neighbour_lncrna_filter, lncrna_type)
# Add gene descriptions to neighbour results
gene_descriptions <- c()
for(gene in neighbour_lncrna$pcg_gene){
  desc <- gene_attr[gene]
  gene_descriptions <- c(gene_descriptions, desc)
}
neighbour_lncrna <-  cbind(neighbour_lncrna, gene_descriptions)


### MAKE A SIMILAR DATAFRAME FOR SENSE/ANTISENSE PAIRS ###
sas_out <- do.call(rbind, Map(data.frame, pcg_transcript = sense_names,
                                    pcg_gene = sense_ids,
                                    antisense_transcript = antisense_names,
                                    distance = rep(0, length(sense_names)),
                                    sqrd_distance = rep(0, length(sense_names)),
                                    lncrna_type = rep("antisense", length(sense_names))))
gene_descriptions <- c()
for(gene in sas_out$pcg_gene){
  desc <- gene_attr[gene]
  gene_descriptions <- c(gene_descriptions, desc)
}
sense_antisense_pairs <-  cbind(sas_out, gene_descriptions)


### PERFORM DE ###


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
  
  tiff(file = paste(c(output_dir,"/log_cluster_heatmap.tiff"),sep="", collapse=""),
       unit = "in",
       res=600,
       height = 7,
       width = 4.5)
  
  print(out)
  
  dev.off()
  
  tiff(file = paste(c(output_dir,"/log_hierarchical_heatmap.tiff"),sep="", collapse=""),
       unit = "in",
       res=600,
       height = 7,
       width = 4.5)
  
  print(out2)
  

  
  dev.off()
  logchange_w_cluster <-  add_column(logchange, cluster = cluster_labels, .before = 1)
  write.csv(logchange_w_cluster, paste(c(output_dir,"/lfc_w_cluster.csv"),sep="", collapse=""), row.names = T)
  
}








# Identify and annotate lncRNA with neighbouring genes
source("get_neighbours.R")
lncRNA_pcg_neighbours <- get_neighbours(gtf, lncRNA_lines, PCG_lines, 5000, 10)

nextdoor_lncrna <- get_neighbours(gtf, lncRNA_lines, PCG_lines, 5000, 1)


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



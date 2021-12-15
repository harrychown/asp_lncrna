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
library("DESeq2")
library(factoextra)
library(tidyr)
library(mvnormtest)
library(scales)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#source("cluster_deseq.R")
source("get_neighbours.R")
source("analyse_neighbours.R")

### DESEQ2 FOR CONDITIONAL DATA ###
getdeseq <- function(indata, countsdir){
  # Build input for DESeq
  sampleTable<-data.frame(sampleName=indata$file, fileName=indata$file, condition=factor(indata$condition))
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=countsdir, design=~condition)
  ddsHTSeq$condition <- factor(ddsHTSeq$condition, levels=unique(ddsHTSeq$condition))
  all_conditions <- levels(ddsHTSeq$condition)
  
  # Perform DESeq
  dds <- DESeq(ddsHTSeq)
  # Number of transcripts
  num_transcripts <- length(dds@rowRanges)
  
  # Generate an empty results table for LFC, p-value and baseMean
  lfc_matrix <- as.data.frame(matrix(data=NA, nrow=num_transcripts, ncol=length(all_conditions) - 1))
  pval_matrix <- as.data.frame(matrix(data=NA, nrow=num_transcripts, ncol=length(all_conditions) - 1))
  basemean_matrix <- as.data.frame(matrix(data=NA, nrow=num_transcripts, ncol=length(all_conditions) - 1))
  sqrd_lfc_matrix <- as.data.frame(matrix(data=NA, nrow=num_transcripts, ncol=length(all_conditions) - 1))
  
  colnames(lfc_matrix) <- colnames(pval_matrix) <- colnames(basemean_matrix) <-  colnames(sqrd_lfc_matrix) <- all_conditions[2:length(all_conditions)]
  
  # Perform DESeq
  dds <- DESeq(ddsHTSeq)
  
  # Extract results from each condition
  for(i in 2:length(all_conditions)){
    res <- results(dds, c("condition", all_conditions[i], all_conditions[1]))
    res_df <- as.data.frame(res)
    rownames(lfc_matrix) <- rownames(pval_matrix) <- rownames(basemean_matrix) <-  rownames(sqrd_lfc_matrix) <- rownames(res_df)
    
    # Store results in previously generated matrices
    lfc_matrix[,i-1] <- res_df$log2FoldChange
    pval_matrix[,i-1] <- res_df$padj
    basemean_matrix[,i-1] <- res_df$baseMean
    sqrd_lfc_matrix[,i-1] <- res_df$log2FoldChange ** 2
  }
  
  # Generate overview dataframe for filtering
  # Min Pval, Max basemean, Max sqrdLFC
  min_pval <- apply(pval_matrix, 1, min)
  max_basemean <- apply(basemean_matrix, 1, max)
  max_sqrd_lfc <- apply(sqrd_lfc_matrix, 1, max)
  # Build dataframe
  best_dataframe <- data.frame(pval=min_pval, basemean=max_basemean, sqrdlfc=max_sqrd_lfc)
  # Store results in a list
  output <- list("lfc"=lfc_matrix, "pval"=pval_matrix, "basemean"=basemean_matrix, "sqrdlfc"=sqrd_lfc_matrix)
  # Filter
  keep_index <- which((best_dataframe$basemean>30&best_dataframe$pval<0.05&best_dataframe$sqrdlfc>1)|best_dataframe$basemean>5000)
  keep_index <- keep_index[! keep_index %in% which(is.na(best_dataframe))]
  keep_rownames <- rownames(lfc_matrix)[keep_index]
  for(i in 1:length(output)){
    dm <- output[[i]]
    print(nrow(dm))
    print(max(keep_index))
    filtered_dm <- dm[keep_index, ]
    rownames(filtered_dm) <- keep_rownames
    output[[i]] <- filtered_dm
  }
  
  return(output)
  
}


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
lncrna_names <- c()
for(i in raw_lncrna){
  l_name <- str_match(i, "transcript_id \"(.*?)\"")[1,2]
  lncrna_names <- c(lncrna_names, l_name)
}


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

### PERFORM DESEQ ON ALL DATA ###
if(file.exists("lncnrna_deseq.rds")){
 all_deseq <- readRDS("lncnrna_deseq.rds")
} else{
  data_prefixes <- c()
  all_deseq <- list()
  counter <- 1
  for(split_file in all_split_files){
    prefix <- str_split(split_file, "_")[[1]][1]
    split_file_path <- paste(c(inputdir, "/split_files/", split_file), sep="", collapse="")
    print(split_file_path)
    count_data <- read.csv(split_file_path, header=T)
    deseq_data <- getdeseq(count_data, count_dir)
    all_deseq[[counter]] <- deseq_data
    names(all_deseq)[counter] <- prefix
    counter <- counter + 1
  }
  saveRDS(all_deseq, file = "lncnrna_deseq.rds")
}


### PERFORM CLUSTERING OF ITRA DATA###
itra_data <- all_deseq$itra
itra_lfc <- itra_data$lfc
# Scale each transcript relative to itself
norm_wide_data <-  t(scale(t(itra_lfc)))
# Generate inverse data
lncrna_id <- grep("MSTRG", rownames(norm_wide_data))
inv_wide_data <- norm_wide_data
inv_wide_data[lncrna_id, ] <- inv_wide_data[lncrna_id, ] * -1
# Perform clustering
set.seed(0)
km <- pheatmap(norm_wide_data, kmeans_k = 20, scale = "row",
               cluster_cols = F, cluster_rows = T)
inv_km <- pheatmap(inv_wide_data, kmeans_k = 20, scale = "row",
               cluster_cols = F, cluster_rows = T)
hierarch <- pheatmap(norm_wide_data, show_rownames=F, cluster_cols=F, cluster_rows=T, 
                clustering_distance_rows="euclidean",
                cex=1, clustering_distance_cols="euclidean", 
                clustering_method="complete", border_color=FALSE)
inv_hierarch <- pheatmap(inv_wide_data, show_rownames=F, cluster_cols=F, cluster_rows=T, 
                    clustering_distance_rows="euclidean",
                    cex=1, clustering_distance_cols="euclidean", 
                    clustering_method="complete", border_color=FALSE)
print(km)
dev.off()
print(hierarch)
dev.off()

### GENERATE CLUSTER GRAPHS ###
# Convert wide data to long
wide_data <-  as.data.frame(norm_wide_data)
# Add transcript ID, cluster number and transcript type (lncrna/pcg)
wide_data$transcript <- rownames(wide_data)
wide_data$cluster <- km$kmeans$cluster
wide_data$inv_cluster <- inv_km$kmeans$cluster
wide_data$type <- rep("pcg", nrow(wide_data))
wide_data$type[grep("MSTRG", wide_data$transcript)] <- "lncrna"

# Generate cluster medians
medianscores <- as.data.frame(matrix(nrow=length(unique(wide_data$cluster)), ncol=5))
for(c in 1:length(unique(wide_data$cluster))){
  medianscores[c, ] <- apply(wide_data[wide_data$cluster==c,1:5] ,2, median)
}
medianscores$cluster <- 1:length(unique(wide_data$cluster))
colnames(medianscores) <- c(0.25, 0.5, 1, 2, 4, "cluster")
medianscores$transcript <-  rep("Median", nrow(medianscores))
median_long <- melt(setDT(medianscores), id.vars = c("cluster", "transcript"), variable.name = "condition")
median_long$condition <- as.numeric(levels(median_long$condition))[median_long$condition]

long_data <- melt(setDT(wide_data), id.vars = c("transcript", "cluster", "inv_cluster", "type"), variable.name = "condition")
# Set the MIC to numeric format

long_data$condition <- as.numeric(levels(long_data$condition))[long_data$condition]

# Plot clusters
p <- ggplot(data=long_data, aes(x=condition, y=value, group=transcript)) +
  geom_line(color="#00BFC4") + 
  geom_line(data=median_long, aes(x=condition, y=value, group=transcript), color="red") +
  facet_wrap(~ cluster) +
  xlab("MIC") + ylab("Normalized Log2Fold Change")
print(p)

### CHECK WHETHER THE NEIGHBOURS ARE IN THE SAME CLUSTER ###
if(file.exists("neighbour_data.csv")){
  neighbour_lncrna_w_cluster <- read.csv("neighbour_data.csv")
}else{
  # Generate a new dataframe with cluster columns
neighbour_lncrna_w_cluster <- neighbour_lncrna
neighbour_lncrna_w_cluster$t_clust <- rep(NA, nrow(neighbour_lncrna_w_cluster))
neighbour_lncrna_w_cluster$n_clust <- rep(NA, nrow(neighbour_lncrna_w_cluster))
neighbour_lncrna_w_cluster$inv_t_clust <- rep(NA, nrow(neighbour_lncrna_w_cluster))
neighbour_lncrna_w_cluster$inv_n_clust <- rep(NA, nrow(neighbour_lncrna_w_cluster))

for(i in 1:nrow(neighbour_lncrna)){
  # Extract transcript names
  t_name <- neighbour_lncrna$pcg_transcript[i]
  n_name <- neighbour_lncrna$transcript_neighbour[i]
  
  # Extract relevant IDs in wide data
  t_wide_id <- which(t_name == wide_data$transcript)
  n_wide_id <- which(n_name == wide_data$transcript)
  # Check transcript exists in DESEQ results, if so add cluster details
  if(length(t_wide_id) > 0){
    t_data <- wide_data[t_wide_id, ]
    neighbour_lncrna_w_cluster$t_clust[i] <- t_data$cluster
    neighbour_lncrna_w_cluster$inv_t_clust[i] <- t_data$inv_cluster
  }
  if(length(n_wide_id) > 0){
    n_data <- wide_data[n_wide_id, ]
    neighbour_lncrna_w_cluster$n_clust[i] <- n_data$cluster
    neighbour_lncrna_w_cluster$inv_n_clust[i] <- n_data$inv_cluster
  }
  
}
write.csv(neighbour_lncrna_w_cluster,"neighbour_data.csv", row.names = FALSE)
}

if(!file.exists("clustered_neighbours.csv")){
  same_clusters <- which(neighbour_lncrna_w_cluster$t_clust == neighbour_lncrna_w_cluster$n_clust)
  same_clust_df <- neighbour_lncrna_w_cluster[same_clusters, ]
  write.csv(same_clust_df, "clustered_neighbours.csv", row.names = F)
  inv_clusters <- which(neighbour_lncrna_w_cluster$inv_t_clust == neighbour_lncrna_w_cluster$inv_n_clust)
  inv_clust_df <- neighbour_lncrna_w_cluster[inv_clusters, ]
  write.csv(inv_clust_df, "inv_clustered_neighbours.csv", row.names = F)
  
}


### ANALYSE SENSE/ANTISENSE PAIRS ###
if(file.exists("sas_data.csv")){
  sas_w_cluster <- read.csv("sas_data.csv")
}else {
sas_w_cluster <- sense_antisense_pairs
sas_w_cluster$t_clust <- rep(NA, nrow(sas_w_cluster))
sas_w_cluster$n_clust <- rep(NA, nrow(sas_w_cluster))
sas_w_cluster$inv_t_clust <- rep(NA, nrow(sas_w_cluster))
sas_w_cluster$inv_n_clust <- rep(NA, nrow(sas_w_cluster))

for(i in 1:nrow(sense_antisense_pairs)){
  # Extract transcript names
  t_name <- sense_antisense_pairs$pcg_transcript[i]
  n_name <- sense_antisense_pairs$antisense_transcript[i]
  
  # Extract relevant IDs in wide data
  t_wide_id <- which(t_name == wide_data$transcript)
  n_wide_id <- which(n_name == wide_data$transcript)
  # Check transcript exists in DESEQ results, if so add cluster details
  if(length(t_wide_id) > 0){
    t_data <- wide_data[t_wide_id, ]
    sas_w_cluster$t_clust[i] <- t_data$cluster
    sas_w_cluster$inv_t_clust[i] <- t_data$inv_cluster
  }
  if(length(n_wide_id) > 0){
    n_data <- wide_data[n_wide_id, ]
    sas_w_cluster$n_clust[i] <- n_data$cluster
    sas_w_cluster$inv_n_clust[i] <- n_data$inv_cluster
  }
  
}
write.csv(sas_w_cluster,"sas_data.csv", row.names = FALSE)
}

if(!file.exists("clustered_sas.csv")){
  same_clusters <- which(sas_w_cluster$t_clust == sas_w_cluster$n_clust)
  same_clust_df <- sas_w_cluster[same_clusters, ]
  write.csv(same_clust_df, "clustered_sas.csv", row.names = F)
  inv_clusters <- which(sas_w_cluster$inv_t_clust == sas_w_cluster$inv_n_clust)
  inv_clust_df <- sas_w_cluster[inv_clusters, ]
  write.csv(inv_clust_df, "inv_clustered_sas.csv", row.names = F)
}



### PERFORM STATISTICAL ANALYSES ###
# For each gene is there a significant difference in it's expression between drugs affecting ergosterol biosynthetic pathway
# MICs from 0.5 to 4 were analysed
#Itra, sim, terb, = ergosterol 
# Extract all datasets and add a condition column
itra_lfc <-  all_deseq$itra$lfc
itra_lfc$condition <- "itra"
itra_lfc$transcript <- rownames(itra_lfc)
itra_lfc <- itra_lfc[,-1]

simv_lfc <- all_deseq$simv$lfc
simv_lfc$condition <- "simv"
simv_lfc$transcript <- rownames(simv_lfc)

terb_lfc <- all_deseq$terb$lfc
terb_lfc$condition <- "terb"
terb_lfc$transcript <- rownames(terb_lfc)

fc_lfc <- all_deseq$`5fc`$lfc
fc_lfc$condition <- "5fc"
fc_lfc$transcript <- rownames(fc_lfc)

dodin_lfc <- all_deseq$dodin$lfc
dodin_lfc$condition <- "dodin"
dodin_lfc$transcript <- rownames(dodin_lfc)

hyg_lfc <- all_deseq$hyg$lfc
hyg_lfc$condition <- "hyg"
hyg_lfc$transcript <- rownames(hyg_lfc)

milt_lfc <- all_deseq$milt$lfc
milt_lfc$condition <- "milt"
milt_lfc$transcript <- rownames(milt_lfc)

all_drug_lfc <- rbind(itra_lfc, simv_lfc, terb_lfc, fc_lfc, dodin_lfc, hyg_lfc, milt_lfc)

# Check correlation between each condition
conditions <- unique(all_drug_lfc$condition)
# Per transcript, check correlation
long_results <- c()
for(transcript in lncrna_names){
  transcript <- "MSTRG.45.2"
  transcript_lfc <- all_drug_lfc[all_drug_lfc$transcript == transcript, ]
  transcript_pcc <- matrix(nrow=length(conditions), ncol=length(conditions))
  row_count <- 0
  # Extract condition 1 
  for(drug_i in conditions){
    row_count <- row_count + 1
    drug_i_lfc <- transcript_lfc[transcript_lfc$condition == drug_i, 1:4]
    # If there is no result for drug i keep values as NA
    
    col_count <- 0
    # Extract condition 2
    for(drug_j in conditions){
      col_count <- col_count + 1
      if(drug_i == drug_j | drug_i_lfc == 0){
        P <- NA
        R2 <- NA
      }else{ 
        drug_j_lfc <- transcript_lfc[transcript_lfc$condition == drug_j, 1:4]
        
        # Perform PCC 
        cor.res <- cor.test(as.numeric(drug_i_lfc), as.numeric(drug_j_lfc))
        P <- cor.res$p.value
        R2 <- cor.res$estimate ** 2
        stop()
        }
      #print(cor(drug_i_lfc, drug_j_lfc, method = c("pearson")))
      
    }
  }
  
}

# We would want to adjust the p-values once we have them all
# This could be done after ALL results for ALL drugs and ALL lncRNA are obtained
# Or, once ALL results for ALL drugs for a single lncRNA is obtained
# I'm going to try both and see what happens
# We can use p.adjust to assess this






# Extract each lncRNA name for subsetting
# Subset each lncRNA
lncrna <- lncrna_names[1]
itra_res <- itra_lfc[lncrna, ]
simv_res <- simv_lfc[lncrna, ]
terb_res <- terb_lfc[lncrna, ]
fc_res <- fc_lfc[lncrna, ]
dodin_res <- dodin_lfc[lncrna, ]
hyg_res <- hyg_lfc[lncrna, ]
milt_res <- milt_lfc[lncrna, ]
all_drug_res <- rbind(itra_res, simv_res, terb_res, fc_res, dodin_res, hyg_res, milt_res)

# Subset to remove drugs which do not appear in the data
sub_drug <- all_drug_res[-which(is.na(all_drug_res$condition)), -ncol(all_drug_res)]
sub_conditions <- sub_drug$condition
drug_t <- t(sub_drug)
colnames(drug_t) <- sub_conditions
# Perform correlation
res <- cor(all_drug_res[,-5])
colnames(res) <- rownames(res) <- sub_conditions


# res.man <- manova(. ~ condition, data = erg_drugs)
# 
# model <- lm(. ~ condition, data = erg_drugs)
# cor(model$model$dd, model$fitted.values)
# 
# sqrt(summary(model)$r.squared)
# 
# erg_drugs <- rbind(itra_lfc, simv_lfc, terb_lfc)#, fc_lfc, dodin_lfc, hyg_lfc, milt_lfc)
# pca_res <- prcomp(erg_drugs[,-5], scale. = TRUE)
# autoplot(pca_res, data = erg_drugs, colour = 'condition')














# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ### PERFORM DE ###
# 
# 
# # Generate logfold change
# source("repo_deseq.R")
# deseq_results <- deseq_func(count_directory = count_folder, split_file_path = split_file,
#                         p_val = 0.05 , min_base_mean = 30 , min_LFC = 1, signif_base_mean = 5000)
# 
# 
# 
# 
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# logchange <- deseq_results[[1]]
# signif_bmean <- deseq_results[[2]]
# raw_bmean <- deseq_results[[3]]
# colnames(raw_bmean) <- colnames(logchange)
# raw_counts <- deseq_results[[4]]
# 
# write.csv(logchange,paste(c(prefix,"_lfc.csv"),sep="", collapse=""), row.names = T)
# # Only generate heatmaps for itra data
# if(prefix == "itra"){
#   # Generate clusters and heatmaps
#   ### GENERATE HEATMAP - Harry Chown ###
#   # Convert data into a matrix which can be used as input for heatmap
#   set.seed(2)
#   gene_names <- logchange[,1]
#   clust_num <- 30
#   out<- pheatmap(logchange, kmeans_k = clust_num, scale = "row", 
#                  cluster_cols = F)
#   out2<- pheatmap(logchange,
#                   show_rownames=F, cluster_cols=F, cluster_rows=T, scale="row",
#                   cex=1, clustering_distance_rows="euclidean", cex=1,
#                   clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE)
#   cluster_labels <- out$kmeans$cluster
#   
#   tiff(file = paste(c(output_dir,"/log_cluster_heatmap.tiff"),sep="", collapse=""),
#        unit = "in",
#        res=600,
#        height = 7,
#        width = 4.5)
#   
#   print(out)
#   
#   dev.off()
#   
#   tiff(file = paste(c(output_dir,"/log_hierarchical_heatmap.tiff"),sep="", collapse=""),
#        unit = "in",
#        res=600,
#        height = 7,
#        width = 4.5)
#   
#   print(out2)
#   
# 
#   
#   dev.off()
#   logchange_w_cluster <-  add_column(logchange, cluster = cluster_labels, .before = 1)
#   write.csv(logchange_w_cluster, paste(c(output_dir,"/lfc_w_cluster.csv"),sep="", collapse=""), row.names = T)
#   
# }
# 
# 
# 
# 
# 
# 
# 
# 
# # Identify and annotate lncRNA with neighbouring genes
# source("get_neighbours.R")
# lncRNA_pcg_neighbours <- get_neighbours(gtf, lncRNA_lines, PCG_lines, 5000, 10)
# 
# nextdoor_lncrna <- get_neighbours(gtf, lncRNA_lines, PCG_lines, 5000, 1)
# 
# 
# antisense_pcg_neighbours <- get_neighbours(gtf, transcript_lines[antisense_id], PCG_lines, 5000, 10)
# intergenic_pcg_neighbours <- get_neighbours(gtf, transcript_lines[intergenic_id], PCG_lines, 5000, 10)
# k_pcg_neighbours <- get_neighbours(gtf, transcript_lines[k_id], PCG_lines, 5000, 10)
# # Identify which neighbour pairs are correlated
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# source("analyse_neighbours.R")
# a_n_stats <- analyse_neighbours("antisense_neighbours", "/home/harry/Documents/lncRNA_21_10/neighbours/", antisense_pcg_neighbours, logchange, raw_counts, bmean_df=raw_bmean)
# i_n_stats <- analyse_neighbours("intergenic_neighbours", "/home/harry/Documents/lncRNA_21_10/neighbours/", intergenic_pcg_neighbours, logchange, raw_counts, bmean_df=raw_bmean)
# k_n_stats <- analyse_neighbours("k_neighbours", "/home/harry/Documents/lncRNA_21_10/neighbours/", k_pcg_neighbours, logchange, raw_counts, bmean_df=raw_bmean)
# 
# DE_a_n <- a_n_stats[[1]]
# # Identify points with the most extreme values
# 
# #p<-ggplot(DE_A_N, aes(x=r_val.cor, y=lncRNA_r_val.cor)) +
# #  geom_point() +
# #  lims(x=c(-1,1),y=c(-1,1)) +
# #  theme_minimal() +
# #  coord_fixed() +  
# #  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) 
# #p
# 
# 
# # Identify which sense/antisense pairs are correlated
# # Iterate through the pairs and extract the corresponding info
# #lncRNA, PCG_id, PCG_transcript, distance (0), ordered_desc
# 
# 
# sas_df <- annotate_pcg(sas_out)
# sas_stats <- analyse_neighbours("sas", "/home/harry/Documents/lncRNA_21_10/s_as_pairs/", sas_df, logchange, raw_counts, bmean_df=raw_bmean)
# 


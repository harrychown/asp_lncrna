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
library(data.table)
library(venn)
library(UpSetR)

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
  norm_lfc_matrix <-  t(scale(t(lfc_matrix)))
  # Store results in a list
  output <- list("lfc"=lfc_matrix, "pval"=pval_matrix, "basemean"=basemean_matrix, "sqrdlfc"=sqrd_lfc_matrix, "normlfc"=norm_lfc_matrix)
  # Filter
  keep_index <- which((best_dataframe$basemean>30&best_dataframe$pval<0.05&best_dataframe$sqrdlfc>1)|best_dataframe$basemean>5000)
  keep_index <- keep_index[! keep_index %in% which(is.na(best_dataframe))]
  keep_rownames <- rownames(lfc_matrix)[keep_index]
  for(i in 1:length(output)){
    dm <- output[[i]]
    filtered_dm <- dm[keep_index, ]
    rownames(filtered_dm) <- keep_rownames
    output[[i]] <- filtered_dm
  }
  
  return(output)
  
}

# Convert GTF to BED
gtf2bed <- function(gtf_in){
  bed_tmp <- tempfile(pattern = "bed_tmp")
  bed_cmd <- paste(c("gtf2bed | sort -k1,1 -k2,2n > ",bed_tmp), sep = "", collapse = "")
  system(bed_cmd, input = gtf_in, intern = FALSE)
  return(bed_tmp)
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

### IDENTIFY MULTI-INTERSECTING ANTISENSE lncRNA ###
# Use Bedtools intersect to validate antisense lncRNA
antisense_bed <- gtf2bed(raw_antisense)
pcg_bed <- gtf2bed(raw_pcg)
bedtools_interesect <- paste(c("bedtools intersect -a ", antisense_bed, " -b ", pcg_bed," -wa -wb -S"), 
                            sep = "", collapse ="")
bedtools_results <- system(bedtools_interesect, intern = TRUE)
a_id <- c()
s_id <- c()
i_bp <- c()
i_percent <- c()
s_t <- c()
for(r in bedtools_results){
  split <- strsplit(r, "\t")
  antisense_desc <-  split[[1]][10]
  antisense <- strsplit(antisense_desc, "\"")[[1]][2]
  sense <-  split[[1]][14]
  sense_desc <- split[[1]][20]
  transcript <- strsplit(sense_desc, "\"")[[1]][2]
  
  a_start <- as.integer(split[[1]][2])
  a_end <- as.integer(split[[1]][3])
  a_size = a_end - a_start
  
  s_start <- as.integer(split[[1]][12])
  s_end <- as.integer(split[[1]][13])
  s_size = s_end - s_start
  if(a_start < s_end & a_start > s_start & a_end > s_end & a_end > s_start){
    intersection = s_end - a_start
  }else if(a_start < s_start & a_start < s_end & a_end > s_start & a_end < s_end){
    intersection = a_end - s_start
  }else if(a_start < s_start & a_start < s_end & s_end > a_start & s_end > a_end){
    intersection = s_size
  }else if(a_start > s_start & a_start < s_end & a_end > s_start & a_end < s_end){
    intersection = a_size

  }

  percentage_i <- (intersection / a_size) * 100 
  i_percent <- c(i_percent, percentage_i)
  i_bp <- c(i_bp, intersection)
  a_id <- c(a_id, antisense)
  s_id <- c(s_id, sense)
  s_t <- c(s_t, transcript)
}
sas_out <- do.call(rbind, Map(data.frame, pcg_transcript = s_t,
                                    pcg_gene = s_id, antisense_transcript = a_id,
                                    intersection_bp = i_bp, gene_percent_intersection = i_percent,
                                    lncrna_type = rep("antisense", length(s_id))))
# Count how many have multiple overlaps
overlap_freq <- as.data.frame(table(sas_out$antisense))
overlap_overview <- as.data.frame(table(overlap_freq$Freq))
                                   

### IDENTIFY NEIGHBOURING GENES ###
# Use Bedtools closest to identify whether there are lncRNA directly next to PCGs

# Number of nearest neighbours = 1
num_neighbours  <-  1

# Convert lncRNA and PCG GTFs to BED format
lncrna_bed <- gtf2bed(raw_lncrna)
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


### ADD GENE DESCRIPTIONS TO S/AS DATA ###
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
  multiclust_list <- list()
  multipresence_list <- list()
  for(split_file in all_split_files){
    prefix <- str_split(split_file, "_")[[1]][1]
    split_file_path <- paste(c(inputdir, "/split_files/", split_file), sep="", collapse="")
    count_data <- read.csv(split_file_path, header=T)
    deseq_data <- getdeseq(count_data, count_dir)
    write.csv(deseq_data$lfc,paste(c(prefix, "_lfc_data.csv"), sep="", collapse=""), row.names = TRUE)
    write.csv(deseq_data$normlfc,paste(c(prefix, "_normlfc_data.csv"), sep="", collapse=""), row.names = TRUE)
    # Add the prefix to the normalised LFC ready for multidrug clustering
    if(prefix %in% c("5fc","dodin","hyg","itra","milt","simv","terb")){
      norm_df <- deseq_data$normlfc
      # Remove 0.25 x MIC from itra
      if(prefix == "itra"){
        norm_df <- norm_df[,-1]
      }
      multiclust <- as.data.frame(cbind(norm_df, rep(prefix, nrow(norm_df)), rownames(norm_df)))
      multiclust_list[[prefix]] <- multiclust
    }
    # Store all LFC results
    multipresence <- as.data.frame(cbind(deseq_data$lfc, rep(prefix, nrow(deseq_data$lfc)), rownames(deseq_data$lfc)))
    multipresence_list[[prefix]] <- multipresence
    all_deseq[[prefix]] <- deseq_data
  }
  multiclust_df <- as.data.frame(rbindlist(multiclust_list))
  all_deseq[["multiclust"]] <- multiclust_df
  all_deseq[["multipresence"]] <- multipresence_list
  saveRDS(all_deseq, file = "lncnrna_deseq.rds")
}

### PERFORM PRESENCE ANALYSIS OF ALL LNCRNA
all_lncrna <- list()
for(i in names(all_deseq$multipresence)){
  condition <- i
  all_lncrna[[condition]] <- rownames(all_deseq$multipresence[[i]])[grep("MSTRG", rownames(all_deseq$multipresence[[i]]))]
}

upset_p <- upset(fromList(all_lncrna), sets = names(all_lncrna), point.size = 3.5, line.size = 2, order.by = "freq", text.scale = c(1.3, 1.3, 1, 1, 2, 1), mainbar.y.label = "Number of lncRNA", sets.x.label = "DE lncRNA Per Condition")

png("present_upset.png", width = 3000, height = 3000, pointsize = 20, res = 300)
print(upset_p)
dev.off()

### PERFORM CLUSTERING OF MULTIDRUG LNCRNA DATA ###
k_all <- pheatmap(all_deseq$multiclust[,1:4], kmeans_k = 20,
               cluster_cols = F, cluster_rows = T)
# Extract only lncRNA
multi_lncrna <-  all_deseq$multiclust[grep("MSTRG", all_deseq$multiclust$V6),]
multi_lncrna_clusters <- k_all$kmeans$cluster[grep("MSTRG", all_deseq$multiclust$V6)]
multidrug_df <- cbind(multi_lncrna, multi_lncrna_clusters)


### GENERATE CLUSTER GRAPHS ###
# Convert wide data to long
MULTIwide_data <-  as.data.frame(multidrug_df[,1:4])
colnames(MULTIwide_data) <- c(0.5, 1, 2, 4)
MULTIwide_data <- data.frame(apply(MULTIwide_data, 2, function(x) as.numeric(as.character(x))))
# Add transcript ID, cluster number and transcript type (lncrna/pcg)
MULTIwide_data$transcript <- paste(multidrug_df$V5, multidrug_df$V6,  sep="_")
MULTIwide_data$cluster <- multi_lncrna_clusters
MULTIwide_data$type <- rep("lncrna", nrow(MULTIwide_data))

# Generate cluster medians
MULTImedianscores <- as.data.frame(matrix(nrow=length(unique(MULTIwide_data$cluster)), ncol=4))
for(c in 1:length(unique(MULTIwide_data$cluster))){
  MULTImedianscores[c, ] <- apply(MULTIwide_data[MULTIwide_data$cluster==c,1:4] ,2, median)
}
MULTImedianscores$cluster <- 1:length(unique(MULTIwide_data$cluster))
colnames(MULTImedianscores) <- c(0.5, 1, 2, 4, "cluster")
MULTImedianscores$transcript <-  rep("Median", nrow(MULTImedianscores))
MULTImedian_long <- melt(setDT(MULTImedianscores), id.vars = c("cluster", "transcript"), variable.name = "condition")

MULTIlong_data <- melt(setDT(MULTIwide_data), id.vars = c("transcript", "cluster", "type"), variable.name = "condition")


# Set the MIC to numeric format
MULTIlong_data$condition <- gsub("X", "", MULTIlong_data$condition) 
#MULTIlong_data$condition <- as.numeric(levels(MULTIlong_data$condition))[MULTIlong_data$condition]

# Plot clusters
p <- ggplot(data=MULTIlong_data, aes(x=condition, y=value, group=transcript)) +
  geom_line(color="#00BFC4") + 
  geom_line(data=MULTImedian_long, aes(x=condition, y=value, group=transcript), color="red") +
  geom_hline(yintercept=0) +
  facet_wrap(~ cluster) +
  xlab("MIC") + ylab("Normalized Log2Fold Change")

png("multidrug_cluster.png", width = 2000, height = 2000, pointsize = 20, res = 300)
print(p)
dev.off()




# Attempt a Venn Diagram
venn_data <- multidrug_df
# Combine the cluster ID with the transcript ID
venn_data$x <- paste(venn_data$V6,venn_data$multi_lncrna_clusters, sep="_")
# Save each transcript per drug
venn_list <- list()
for(i in unique(venn_data$V5)){
  venn_list[[i]] <- venn_data$x[venn_data$V5 == i]
}
venn_out <- venn(venn_list)
venn_overview <- as.data.frame(cbind(rownames(venn_out),venn_out$counts))[-1,]
colnames(venn_overview) <- c("combination", "count")

upset_p <- upset(fromList(venn_list), sets = names(venn_list),
                 queries = list(list(query = intersects, params = list("simv", "terb"), color = "orange", active = T),
                                list(query = intersects, params = list("simv", "terb", "milt"), color = "green", active = T)),
                 point.size = 3.5, line.size = 2, 
                 order.by = "freq", text.scale = c(1.3, 1.3, 1, 1, 2, 1), 
                 mainbar.y.label = "Number of lncRNA", sets.x.label = "DE lncRNA Per Condition")

# Look at what is shared between the drugs
x1 <- unlist(venn_list, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
# Generate an empty matrix for each combination
cluster_frequency <- matrix(nrow = length(venn_overview$combination), ncol = 20, data = 0)
exact_cluster_frequency <- matrix(nrow = length(venn_overview$combination), ncol = 20, data = 0)
for(combo_i in 1:length(venn_overview$combination)){
  combo <- venn_overview$combination[combo_i]
  split_combo <- str_split(combo, ":")
  combo_binary <- upset_p$New_data[split_combo[[1]]]
  # Compounding value transcripts
  combo_transcripts <- x1[(rowSums(combo_binary) == ncol(combo_binary))]
  freq_table <- table(unlist(lapply(combo_transcripts, function(x) as.numeric(str_split(x, "_")[[1]][2]))))
  cluster_ids <- names(freq_table)
  for(cluster_i in 1:length(cluster_ids)){
    cluster_frequency[combo_i, as.numeric(cluster_ids[cluster_i])] <- as.numeric(freq_table[cluster_i])
  }
  
  # Exact value transcripts
  all_transcripts <- x1[(rowSums(upset_p$New_data) == ncol(combo_binary))]
  exact_transcripts <- combo_transcripts[combo_transcripts %in% all_transcripts]
  exfreq_table <- table(unlist(lapply(exact_transcripts, function(x) as.numeric(str_split(x, "_")[[1]][2]))))
  excluster_ids <- names(exfreq_table)
  for(cluster_i in 1:length(excluster_ids)){
    exact_cluster_frequency[combo_i, as.numeric(excluster_ids[cluster_i])] <- as.numeric(exfreq_table[cluster_i])
  }
  
}
cluster_freq_df <- as.data.frame(cluster_frequency, row.names = venn_overview$combination)
colnames(cluster_freq_df) <- as.character(c(1:20))
cluster_freq_df$total <- rowSums(cluster_freq_df)

excluster_freq_df <- as.data.frame(exact_cluster_frequency, row.names = venn_overview$combination)
colnames(excluster_freq_df) <- as.character(c(1:20))
excluster_freq_df$total <- rowSums(excluster_freq_df)


# Save cluster frequency matrix
write.csv(cluster_freq_df, "multidrug_cluster_frequency.csv")
write.csv(excluster_freq_df, "multidrug_exact_cluster_frequency.csv")
png("upset.png", width = 3000, height = 3000, pointsize = 20, res = 300)
print(upset_p)
dev.off()

### PERFORM CLUSTERING OF ITRA DATA ###
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
km <- pheatmap(norm_wide_data, kmeans_k = 20,
               cluster_cols = F, cluster_rows = T)
inv_km <- pheatmap(inv_wide_data, kmeans_k = 20,
               cluster_cols = F, cluster_rows = T)
hierarch <- pheatmap(norm_wide_data, show_rownames=F, cluster_cols=F, cluster_rows=T, 
                clustering_distance_rows="euclidean",
                cex=1, clustering_distance_cols="euclidean", 
                clustering_method="complete", border_color=FALSE)
inv_hierarch <- pheatmap(inv_wide_data, show_rownames=F, cluster_cols=F, cluster_rows=T, 
                    clustering_distance_rows="euclidean",
                    cex=1, clustering_distance_cols="euclidean", 
                    clustering_method="complete", border_color=FALSE)
png("k_means_itra.png", width = 1600, height = 2000, pointsize = 20, res = 300)
print(km)
dev.off()
png("hierarch_itra.png", width = 1600, height = 2000, pointsize = 20, res = 300)
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
  geom_hline(yintercept=0) +
  facet_wrap(~ cluster) +
  xlab("MIC") + ylab("Normalized Log2Fold Change")

png("cluster.png", width = 2000, height = 2000, pointsize = 20, res = 300)
print(p)
dev.off()

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


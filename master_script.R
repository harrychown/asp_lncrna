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
library(wrapr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#set.seed(0)

#source("get_neighbours.R")
#source("analyse_neighbours.R")

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
  all_basemean <- output$basemean
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

  output[["all_basemean"]] <- all_basemean
  
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
  bmean_list <- list()
  for(split_file in all_split_files){
    prefix <- str_split(split_file, "_")[[1]][1]
    split_file_path <- paste(c(inputdir, "/split_files/", split_file), sep="", collapse="")
    count_data <- read.csv(split_file_path, header=T)
    deseq_data <- getdeseq(count_data, count_dir)
    write.csv(deseq_data$lfc,paste(c(prefix, "_lfc_data.csv"), sep="", collapse=""), row.names = TRUE)
    write.csv(deseq_data$normlfc,paste(c(prefix, "_normlfc_data.csv"), sep="", collapse=""), row.names = TRUE)
    # Save all basemeans for each condition
    bmean_raw <- deseq_data$all_basemean
    bmean_df <- as.data.frame(cbind(bmean_raw[,1],rownames(bmean_raw), rep(prefix, nrow(bmean_raw))))
    write.csv(bmean_df,paste(c(prefix, "_basemean_data.csv"), sep="", collapse=""), row.names = TRUE)
    bmean_list[[prefix]] <- bmean_df
    
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
  multibmean_df <- as.data.frame(rbindlist(bmean_list))
  multiclust_df <- as.data.frame(rbindlist(multiclust_list))
  all_deseq[["multibmean"]] <- multibmean_df
  all_deseq[["multiclust"]] <- multiclust_df
  all_deseq[["multipresence"]] <- multipresence_list
  saveRDS(all_deseq, file = "lncnrna_deseq.rds")
}

### PERFORM DE ANALYSIS OF ALL LNCRNA
all_lncrna <- list()
for(i in names(all_deseq$multipresence)){
  condition <- i
  all_lncrna[[condition]] <- rownames(all_deseq$multipresence[[i]])[grep("MSTRG", rownames(all_deseq$multipresence[[i]]))]
}

upset_de <- upset(fromList(all_lncrna), sets = names(all_lncrna), point.size = 3.5, line.size = 2, order.by = "freq", text.scale = c(1.3, 1.3, 1, 1, 2, 1), mainbar.y.label = "Number of lncRNA", sets.x.label = "DE lncRNA Per Condition")

png("de_upset.png", width = 3000, height = 3000, pointsize = 20, res = 300)
print(upset_de)
dev.off()
#set.seed(0)



### PERFORM PRESENCE ANALYSIS OF ALL LNCRNA
present_lncrna <- list()
lncrna_bmean <- all_deseq$multibmean[grep("MSTRG", (all_deseq$multibmean$V2)), ]
for(i in unique(lncrna_bmean$V3)){
    condition_bmean <- lncrna_bmean[lncrna_bmean$V3 == i,]
  condition_bmean$V1 <- as.numeric(condition_bmean$V1)
  condition_present <- condition_bmean[condition_bmean$V1>30,]
  present_lncrna[[i]] <- condition_present[,2]
}
upset_present <- upset(fromList(present_lncrna), sets = names(present_lncrna), point.size = 3.5, line.size = 2, order.by = "freq", text.scale = c(1.3, 1.3, 1, 1, 2, 1), mainbar.y.label = "Number of lncRNA", sets.x.label = "Transcribed lncRNA Per Condition")

# Look at what lncRNA are present between conditions
up_present_lncrna <- unlist(present_lncrna, use.names = FALSE)
up_present_lncrna <- up_present_lncrna[ !duplicated(up_present_lncrna) ]
rownames(upset_present$New_data) <- up_present_lncrna

png("present_upset.png", width = 3000, height = 3000, pointsize = 20, res = 300)
print(upset_present)
dev.off()



### PERFORM CLUSTERING OF MULTIDRUG LNCRNA DATA ###
multi_lncrna <- all_deseq$multiclust[grep("MSTRG", all_deseq$multiclust$V6),]
multi_data <- apply(multi_lncrna[,1:4], 2, as.numeric)
multi_data <-  t(scale(t(multi_data)))


source("clustering.R")

st_r <- list("strict"=40, "relaxed"=20)
multi_str <- list()
# Perform clustering with strict/relaxed K-value
for(i in 1:length(st_r)){
  output_name <- paste(c("multidrug_", names(st_r[i])), sep="", collapse="")
  multidrug_clust <- cluster_data(multi_data, as.numeric(unname(st_r[i])), 1, output_name)
  multidrug_df <- cbind(multi_lncrna, multidrug_clust)
  # Generate line graphs from cluster data
  multi_graph(multidrug_df, output_name)
  # Generate upset plots of drug combinations
  multi_p <- multi_shared(multidrug_df, output_name)
  multi_str[[names(st_r[i])]] <- multi_p
}

### IDENTIFY WHICH LNCRNA ARE FOUND IN THE STRICT/RELAXED MULTIDRUG CLUSTERING ###
strict_multidrug <- extract_multi_clust_info(multi_str$strict$New_data, multi_str$strict$x1)
relaxed_multidrug <- extract_multi_clust_info(multi_str$relaxed$New_data, multi_str$relaxed$x1)


### PERFORM CLUSTERING OF ITRA DATA ###
itra_data <- all_deseq$itra
itra_lfc <- itra_data$lfc
# Scale each transcript relative to itself
norm_wide_data <-  t(scale(t(itra_lfc)))
# Generate inverse data
lncrna_id <- grep("MSTRG", rownames(norm_wide_data))
inv_wide_data <- norm_wide_data
inv_wide_data[lncrna_id, ] <- inv_wide_data[lncrna_id, ] * -1


source("clustering.R")


# Perform clustering with strict/relaxed K-value
itra_str <- list()
for(i in 1:length(st_r)){
  output_name <- paste(c("itra_", names(st_r[i])), sep="", collapse="")
  itra_clust <- cluster_data(norm_wide_data, as.numeric(unname(st_r[i])), 1, output_name)
  itra_inv <- cluster_data(inv_wide_data, as.numeric(unname(st_r[i])), 1, output_name)
  itra_df <- cbind(norm_wide_data, itra_clust)
  itra_w_cluster <- itra_graph(itra_df, itra_inv, output_name)
  itra_w_cluster <- as.data.frame(itra_w_cluster)
  # Check whether neighbours/sas are in the same cluster
  # And check similarities in the multidrug data
  matched_neighbours <- match_clusters(neighbour_lncrna, itra_w_cluster)
  matched_sas <- match_clusters(sense_antisense_pairs, itra_w_cluster)

}
  






  
  
  



# ### CHECK TO SEE IF NEIGHBOURS ARE IN THE SAME CLUSTER ###
# strict_neighbours <- match_clusters(neighbour_lncrna, itra_str$strict)
# relaxed_neighbours <- match_clusters(neighbour_lncrna, itra_str$relaxed)
# 
# ### CHECK TO SEE IF SAS PAIRS ARE IN THE SAME CLUSTER ###
# strict_sas <- match_clusters(sense_antisense_pairs, itra_str$strict)
# relaxed_sas <- match_clusters(sense_antisense_pairs, itra_str$relaxed)







# ### CHECK IF LNCRNA OF INTEREST ARE PRESENT IN CLUSTERS OR CONDITIONS ###
# add_presence <- function(input_cluster_df, input_matrix){
#   condition_combinations <- c()
# condition_numbers <- c()
# 
# for(lncrna in input_cluster_df[,3]){
#   if(lncrna %in% rownames(input_matrix)){
#     lncrna_presence <- input_matrix[lncrna,]
#     combo <- paste(colnames(lncrna_presence)[which(lncrna_presence == 1)], collapse=":")
#     number_of_conditions <- unname(rowSums(lncrna_presence))
# 
#   }else{
#     combo <- 0
#     number_of_conditions <- 0
#      }
#  condition_combinations <- c(condition_combinations, combo)
#   condition_numbers <- c(condition_numbers, number_of_conditions)
# }
# output_cluster_df <- cbind(input_cluster_df, condition_numbers, condition_combinations)
# return(output_cluster_df)
# }
# 
# 
# 
# present_same_neigh_df <- add_presence(same_neigh_cluster_df, upset_present$New_data)
# present_inv_neigh_df <- add_presence(inv_neigh_cluster_df, upset_present$New_data)
# 
# present_same_sas_df <- add_presence(same_sas_cluster_df, upset_present$New_data)
# present_inv_sas_df <- add_presence(inv_sas_cluster_df, upset_present$New_data)
# 
# write.csv(present_same_neigh_df, "clustered_neighbours.csv", row.names = F)
# write.csv(present_inv_neigh_df, "inv_clustered_neighbours.csv", row.names = F)
# write.csv(present_same_sas_df, "clustered_sas.csv", row.names = F)
# write.csv(present_inv_sas_df, "inv_clustered_sas.csv", row.names = F)


# ### IDENTIFY WHICH LNCRNA ARE FOUND IN SIMILAR MULTIDRUG CLUSTERS ###
# multidrug_cluster_df <- multi_p$New_data
# rownames(multidrug_cluster_df) <- multi_p$x1
# multidrug_num <- as.character(unname(rowSums(multidrug_cluster_df)))
# multidrug_clust <- unlist(lapply(multi_p$x1, function(x) unlist(as.character(str_split(x, "_")[[1]][2]))))
# multidrug_names <- unlist(lapply(multi_p$x1, function(x) unlist(as.character(str_split(x, "_")[[1]][1]))))
# multidrug_combo <- c()
# for(i in 1:nrow(multidrug_cluster_df)){
#   lncrna_presence <- multidrug_cluster_df[i, ]
#   combo <- paste(colnames(lncrna_presence)[which(lncrna_presence == 1)], collapse=":")
#   multidrug_combo <- c(multidrug_combo, combo)
# }
# multidrug_lncrna_df <- as.data.frame(cbind(multidrug_names, multidrug_clust, multidrug_num, multidrug_combo))
# write.csv(multidrug_lncrna_df, "multidrug_clustered_lncrna.csv", row.names = F)


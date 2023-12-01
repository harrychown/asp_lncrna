# Master Script version 2
# DESeq analysis on RNA-seq mutlidrug experiments
library(argparser, quietly = T)
library(stringr, quietly = T)
library(biomaRt, quietly = T) # version 2,50.3
library(DESeq2, quietly = T) # version 1.34.0
library(data.table, quietly = T) # version 1.14.2
library(ggplot2, quietly = T) 
#### DESEQ2 ####
getdeseq <- function(indata, countsdir, deprefix, deout){
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
  res_df <- data.frame()
  # Extract results from each condition
  for(i in 2:length(all_conditions)){
    res <- results(dds, c("condition", all_conditions[i], all_conditions[1]))
    res_df <- as.data.frame(res)
    rownames(lfc_matrix) <- rownames(pval_matrix) <- rownames(basemean_matrix) <-  rownames(sqrd_lfc_matrix) <- rownames(res_df)
    write.csv(res_df, paste(c(argv$outdir, "/", deprefix, "_raw_", all_conditions[i], ".csv"), sep="", collapse=""), row.names = T)
    # Store results in previously generated matrices
    pval_matrix[,i-1] <- res_df$padj
    basemean_matrix[,i-1] <- res_df$baseMean
    sqrd_lfc_matrix[,i-1] <- res_df$log2FoldChange ** 2
    lfc_matrix[,i-1] <- res_df$log2FoldChange
  }
  write.csv(pval_matrix, paste(c(argv$outdir, "/", deprefix, "_raw_pval", ".csv"), sep="", collapse=""), row.names = T)
  # Generate overview dataframe for filtering
  # Min Pval, Max basemean, Max sqrdLFC
  min_pval <- apply(pval_matrix, 1, min)
  max_basemean <- apply(basemean_matrix, 1, max)
  max_sqrd_lfc <- apply(sqrd_lfc_matrix, 1, max)
  # Build dataframe
  best_dataframe <- data.frame(pval=min_pval, basemean=max_basemean, sqrdlfc=max_sqrd_lfc)
  norm_lfc_matrix <-  t(scale(t(lfc_matrix)))
  write.csv(pval_matrix, paste(c(argv$outdir, "/", deprefix, "_raw_pval", ".csv"), sep="", collapse=""), row.names = T)
  write.csv(lfc_matrix, paste(c(argv$outdir, "/", deprefix, "_raw_lfcl", ".csv"), sep="", collapse=""), row.names = T)
  write.csv(norm_lfc_matrix, paste(c(argv$outdir, "/", deprefix, "_raw_stdlfc", ".csv"), sep="", collapse=""), row.names = T)
  
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

#### GTF2BED ####
gtf2bed <- function(gtf_in){
  bed_tmp <- tempfile(pattern = "bed_tmp")
  bed_cmd <- paste(c("gtf2bed | sort -k1,1 -k2,2n > ",bed_tmp), sep = "", collapse = "")
  system(bed_cmd, input = gtf_in, intern = FALSE)
  return(bed_tmp)
}

#### INPUT/OUTPUT ####
parser <- arg_parser('Analyse DESeq data')
parser <- add_argument(parser, 'splitdir',  type = "character",
                    help='Directory containing split files')
parser <- add_argument(parser, 'countdir', type = "character",
                    help='Directory containing count files')
parser <- add_argument(parser, 'gtf', type = "character",
                    help='GTF file')
parser <- add_argument(parser, 'outdir', type = "character",
                       help='Output directory')
parser <- add_argument(parser, '--rds', type = "character",
                       help='Path to previously stored DESeq .RDS file')
parser <- add_argument(parser, '--kclust', type = "character",
                       help='(K) Number of k-means clusters for hierarchical-kmeans',
                       default = "40")
parser <- add_argument(parser, '--iterk', type = "character",
                       help='Number of iterations for k-means',
                       default = "10")
argv <- parse_args(parser)
# Add forward slash to folders if none are present
print(argv$gtf)


#### EXTRACT GTF INFO ####
### EXTRACT LNCRNA GROUPS FROM GTF ###
# Extract transcript lines

gtf <- readLines(argv$gtf)
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

#### CHECKPOINT - STATS ####
print(paste(c("Total number of intergenic: ",as.character(length(raw_intergenic))), sep="", collapse=""), quote=F)
print(paste(c("Total number of antisense: ",as.character(length(raw_antisense))), sep="", collapse=""), quote=F)

#### OBTAIN DESCRIPTIONS FOR ALL PCGs ####
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

#### IDENTIFY MULTI-INTERSECTING ANTISENSE lncRNA ####
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
colnames(overlap_overview) <- c("Overlaps", "Freq")
print("Number of antisense with multiple sense partners:", quote=F)
print(overlap_overview)
write.csv(overlap_overview, file=paste(c(argv$outdir, "/multiantisense.csv"), sep="", collapse=""), row.names=F)


#### IDENTIFY NEIGHBOURING GENES ####
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
write(closest_results, "closest.debug")
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
print(class(neighbour_out))
write.csv(neighbour_out, "neighbour_out.debug")
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
write.csv(neighbour_lncrna, "neighbour_lncrna.debug")

print(paste(c("Total number of PCGs with a lncRNA neighbour: ",as.character(nrow(neighbour_lncrna))), sep="", collapse=""), quote=F)

#### ADD GENE DESCRIPTIONS TO S/AS DATA ####
gene_descriptions <- c()
for(gene in sas_out$pcg_gene){
  desc <- gene_attr[gene]
  gene_descriptions <- c(gene_descriptions, desc)
}
sense_antisense_pairs <-  cbind(sas_out, gene_descriptions)
print(paste(c("Total number of sense/antisense pairs: ",as.character(nrow(sense_antisense_pairs))), sep="", collapse=""), quote=F)


#### PERFORM DESEQ ANALYSIS ####
if(is.na(argv$rds)){
  data_prefixes <- c()
  # Extract split filenames
  all_split_files <- list.files(path=argv$splitdir)
  all_deseq <- list()
  multiclust_list <- list()
  multipresence_list <- list()
  bmean_list <- list()
  for(split_file in all_split_files){
    prefix <- str_split(split_file, "_")[[1]][1]
    split_file_path <- paste(c(argv$splitdir, "/", split_file), sep="", collapse="")
    count_data <- read.csv(split_file_path, header=T)
    deseq_data <- getdeseq(count_data, argv$countdir, prefix, argv$outdir)
    write.csv(deseq_data$lfc,paste(c(argv$outdir, "/", prefix, "_lfc_data.csv"), sep="", collapse=""), row.names = TRUE)
    write.csv(deseq_data$normlfc,paste(c(argv$outdir, "/", prefix, "_normlfc_data.csv"), sep="", collapse=""), row.names = TRUE)
    # Save all basemeans for each condition
    bmean_raw <- deseq_data$all_basemean
    bmean_df <- as.data.frame(cbind(bmean_raw[,1],rownames(bmean_raw), rep(prefix, nrow(bmean_raw))))
    write.csv(bmean_df,paste(c(argv$outdir, "/", prefix, "_basemean_data.csv"), sep="", collapse=""), row.names = TRUE)
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
  saveRDS(all_deseq, file = paste(c(as.character(Sys.Date()), ".rds"), sep="", collapse=""))
} else{
    all_deseq <- readRDS(argv$rds)
}

#### PERFORM CLUSTERING OF MULTIDRUG LNCRNA DATA ####
source("clustering.R")
# Extract lncRNA data for all conditions
multi_lncrna <- all_deseq$multiclust[grep("MSTRG", all_deseq$multiclust$V6),]
multi_data <- apply(multi_lncrna[,1:4], 2, as.numeric)
# Set output filename
output_name <- paste(c(argv$outdir, "/multidrug_k", as.character(argv$k)) , sep="", collapse="")
# Cluster multidrug data
multidrug_clust <- cluster_data(multi_data, as.numeric(argv$k), 1, output_name, as.numeric(argv$iterk))
multidrug_df <- cbind(multi_lncrna, multidrug_clust)
# Generate line graphs from cluster data
multi_graph(multidrug_df, output_name)
# Generate upset plots of drug combinations
multi_p <- multi_shared(multidrug_df, output_name)
multidrug_data_full <- extract_multi_clust_info(multi_p$New_data, multi_p$x1)


#### PERFORM CLUSTERING OF ITRA DATA ####
itra_data <- all_deseq$itra
itra_normlfc <- itra_data$normlfc
# Generate inverse data
lncrna_id <- grep("MSTRG", rownames(itra_normlfc))
inv_wide_data <- itra_normlfc
inv_wide_data[lncrna_id, ] <- inv_wide_data[lncrna_id, ] * -1
# Perform clustering
output_name <- paste(c(argv$outdir, "/itra_k", as.character(argv$k)) , sep="", collapse="")
output_inv <- paste(c(argv$outdir, "/itra_inv_k", as.character(argv$k)) , sep="", collapse="")
itra_clust <- cluster_data(itra_normlfc, as.numeric(argv$k), 1, output_name, as.numeric(argv$iterk))
itra_inv <- cluster_data(inv_wide_data, as.numeric(argv$k), 1, output_inv, as.numeric(argv$iterk))
itra_df <- cbind(itra_normlfc, itra_clust)
itra_w_cluster <- itra_graph(itra_df, itra_inv, output_name)
itra_w_cluster <- as.data.frame(itra_w_cluster)
# Check whether neighbours/sas are in the same cluster

match_clusters <- function(desc_data, clust_data){
  clust_info <- list("t_clust"=c(),
                     "n_clust"=c(),
                     "inv_t_clust"=c(),
                     "inv_n_clust"=c())
  for(i in 1:nrow(desc_data)){
    # Extract pcg/lnc transcript ID
    pcg_transcript <- desc_data[i,1] 
    lnc_transcript <- desc_data[i,3]
    # Search for clusters
    # If the sequence is not DE then replace values with NA else get clusters
    if(length(which(clust_data$transcript == pcg_transcript)) == 0){
      clust_info$t_clust <- c(clust_info$t_clust, NA)
      clust_info$inv_t_clust <- c(clust_info$inv_t_clust, NA)
    } else{
      clust_info$t_clust  <- c(clust_info$t_clust, clust_data[which(clust_data$transcript == pcg_transcript), "cluster"])
      clust_info$inv_t_clust <- c( clust_info$inv_t_clust, clust_data[which(clust_data$transcript == pcg_transcript), "inv_cluster"])
    }
    
    # Do the same for the lncRNA
    if(length(which(clust_data$transcript == lnc_transcript)) == 0){
      clust_info$n_clust <- c(clust_info$n_clust, NA)
      clust_info$inv_n_clust <- c(clust_info$inv_n_clust, NA)
    } else{
      clust_info$n_clust  <- c(clust_info$n_clust, clust_data[which(clust_data$transcript == lnc_transcript), "cluster"])
      clust_info$inv_n_clust <- c(clust_info$inv_n_clust, clust_data[which(clust_data$transcript == lnc_transcript), "inv_cluster"])
    }
    
  }
  clust_info_df <- as.data.frame(do.call(cbind, clust_info))
  input_data_w_cluster <- cbind(desc_data, clust_info_df)
  
  
  same_clusters <- which(input_data_w_cluster$t_clust == input_data_w_cluster$n_clust)
  same_cluster_df <- input_data_w_cluster[same_clusters, ]
  
  inv_clusters <- which(input_data_w_cluster$inv_t_clust == input_data_w_cluster$inv_n_clust)
  inv_cluster_df <- input_data_w_cluster[inv_clusters, ]
  
  outlist <- list(raw_cluster_data=input_data_w_cluster,
                  same_clust_data=same_cluster_df,
                  inv_clust=inv_cluster_df)
  
  
  return(outlist)
  
}

matched_neighbours <- match_clusters(neighbour_lncrna, itra_w_cluster)
write.csv(matched_neighbours$same_clust_data, paste(c(output_name, "_same_neighbours.csv"), sep="", collapse=""))
write.csv(matched_neighbours$inv_clust, paste(c(output_name, "_inv_neighbours.csv"), sep="", collapse=""))
write.csv(matched_neighbours$raw_cluster_data, paste(c(output_name, "_neighbours.csv"), sep="", collapse=""))


matched_sas <- match_clusters(sense_antisense_pairs, itra_w_cluster)
write.csv(matched_sas$same_clust_data, paste(c(output_name, "_same_sas.csv"), sep="", collapse=""))
write.csv(matched_sas$inv_clust, paste(c(output_name, "_inv_sas.csv"), sep="", collapse=""))
write.csv(matched_sas$raw_cluster_data, paste(c(output_name, "_sas.csv"), sep="", collapse=""))


#### IDENTIFY WHICH LNCRNA ARE FOUND IN SIMILAR MULTIDRUG CLUSTERS ####
get_lncrna_multidrug_clust <- function(upset_data, upset_names, outname){
  rownames(upset_data) <- upset_names
  multidrug_num <- as.character(unname(rowSums(upset_data)))
  multidrug_clust <- unlist(lapply(upset_names, function(x) unlist(as.character(str_split(x, "_")[[1]][2]))))
  multidrug_names <- unlist(lapply(upset_names, function(x) unlist(as.character(str_split(x, "_")[[1]][1]))))
  multidrug_combo <- c()
  for(i in 1:nrow(upset_data)){
    lncrna_presence <- upset_data[i, ]
    combo <- paste(colnames(lncrna_presence)[which(lncrna_presence == 1)], collapse=":")
    multidrug_combo <- c(multidrug_combo, combo)
  }
  multidrug_lncrna_df <- as.data.frame(cbind(multidrug_names, multidrug_clust, multidrug_num, multidrug_combo))
  write.csv(multidrug_lncrna_df, outname, row.names = F)
  return(multidrug_lncrna_df)
}

multidrug_k_lncrna <- get_lncrna_multidrug_clust(multi_p$New_data, multi_p$x1, paste(c(argv$outdir, "/multidrug_k", as.character(argv$k), "_lncrna.csv"), sep="", collapse=""))


#### INVESTIGATE THE CO-OCCURENCE OF SIMILAR S/AS AND NEIGHBOURS ####
matched_list_names <- c("same_clust_data","inv_clust" )
 for(i in matched_list_names){
   # Extract neighbour data
   N_res <- matched_neighbours[[i]]
   N_field <- strsplit(i, "_")[[1]][1]
   for(j in matched_list_names){
     SN_res_vec <- c()
     # Extract s/as data
     S_res <- matched_sas[[j]]
     S_field <- strsplit(j, "_")[[1]][1]
     # Identify shared lncRNA
     shared_lnc <- intersect(N_res$transcript_neighbour, S_res$antisense_transcript)
     # Extract pairs and neighbours
     if(length(shared_lnc)>0){
       for(lnc in shared_lnc){
         sense <- S_res$pcg_transcript[S_res$antisense_transcript==lnc]
         neighbours <- N_res$pcg_transcript[N_res$transcript_neighbour==lnc]
         SN_res <- c(lnc, sense, neighbours)
         SN_transcript_line <- paste(SN_res, sep="", collapse=",")
          SN_res_vec <- c(SN_res_vec, SN_transcript_line)
       }
       filename <- paste(c(argv$outdir, "/itra_k", as.character(argv$k), "_N-N_", i, "_and_SAS_", j, ".txt"), sep = "", collapse = "")
       write(SN_res_vec, filename, sep="")
     }

   }
}


# matched_list_names <- c("same_clust_data","inv_clust" )
# for(i in matched_list_names){
#   # Extract neighbour data
#   N_res <- matched_neighbours[[i]]
#   # i prefix
#   N_field <- strsplit(i, "_")[[1]][1]
#   for(j in matched_list_names){
#     # Extract s/as data
#     S_res <- matched_sas[[j]]
#     S_field <- strsplit(j, "_")[[1]][1]
#     # Identify antisense lncRNA in neighbours
# 
#   # S_res_matchedN <-  S_res[which(S_res[,"antisense_transcript"] %in% N_res[,"transcript_neighbour"]),]
#   S_res_matchedN <- S_res[which(S_res[,"antisense_transcript"] %in% N_res[,"transcript_neighbour"]),]
#   N_res_matchedS <- N_res[which(N_res[,"transcript_neighbour"] %in% S_res[,"antisense_transcript"]),]
#   filename <- paste(c(argv$outdir, "/itra_k", as.character(argv$k), "_N-N_", i, "_and_SAS_", j, ".csv"), sep = "", collapse = "")
#   filename2 <- paste(c(argv$outdir, "/itra_k", as.character(argv$k), "_N-N_", i, "_and_SAS_", j, "Nres.csv"), sep = "", collapse = "")
#   # Since a lncRNA may have multiple neighbours we can must show each of these
#   # Output will be the following format: lncRNA, sense, neighbour1, neighbour2
#   if(nrow(S_res_matchedN)>0){
#     sas_nn_results <- list("lncRNA"=c(),
#                            "sense"=c(),
#                            "neighbour1"=c(),
#                            "neighbour2"=c())
#     for(sas_line in 1:nrow(S_res_matchedN)){
#       lncRNA <- S_res_matchedN[sas_line, "antisense_transcript"]
#       sas_nn_results$lncRNA <- c(sas_nn_results$lncRNA, lncRNA)
#       sense <- S_res_matchedN[sas_line, "pcg_transcript"]
#       sas_nn_results$sense <- c(sas_nn_results$sense, sense)
#       n_n_data <- N_res_matchedS[which(lncRNA == N_res_matchedS$transcript_neighbour),]
#       if(nrow(n_n_data) > 1){
#         for(i_n in 1:nrow(n_n_data)){
#           neighbourX <- n_n_data$pcg_transcript[i_n]
#           print(sas_nn_results[i_n + 2])
#           print(neighbourX)
#           print(c(sas_nn_results[i_n + 2], neighbourX))
#           print(class(sas_nn_results[i_n + 2]))
#           print(class(c(sas_nn_results[i_n + 2], neighbourX)))
#           sas_nn_results[[i_n + 2]] <- c(sas_nn_results[[i_n + 2]], neighbourX)
#         }
#       }else{
#         sas_nn_results$neighbour1 <- c(sas_nn_results$neighbour1, n_n_data$pcg_transcript)
#         sas_nn_results$neighbour2 <- c(sas_nn_results$neighbour2, NA)
#       }
#     }
#     
#   }
# 
#   sas_nn_results_df <- as.data.frame(do.call(cbind, sas_nn_results))
#   print(sas_nn_results_df)
#   write.csv(sas_nn_results_df, filename)
#   }
# }
### LINE GRAPH FUNCTION ###
# Supply transcript IDs
#c("transcript:EDP48075", "MSTRG.6767.A")


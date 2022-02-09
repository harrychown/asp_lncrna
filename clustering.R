### PERFORM CLUSTERING ###
iterative_kmeans <- function(input_data, k_val, iter_num, k_alg="LLoyd"){
  kmeans_results <- c()
  kmeans_withinss <- c()
  for(i in 1:iter_num){
    set.seed(i)
    km_out <-  kmeans(input_data, k_val, iter.max = 1000, nstart = 5*k_val, algorithm=k_alg)
    km_wss <- km_out$tot.withinss
    kmeans_withinss <- c(kmeans_withinss, km_wss)
  }
  seed_val <- which(kmeans_withinss==min(kmeans_withinss))
  set.seed(seed_val)
  km_out <-  kmeans(input_data, k_val, iter.max = 1000, nstart = 5*k_val, algorithm=k_alg)
  return(c(km_out, seed_val))
}

cluster_data <- function(input_data, k_val, treeheight, fileprefix, k_alg="LLoyd"){
  library(pheatmap)
  #set.seed(0)
  # Perform initial hierarchical clustering
  hierarch <- pheatmap(input_data, show_rownames=F, cluster_cols=F, cluster_rows=T, 
                       clustering_distance_rows="euclidean",
                       cex=1, clustering_distance_cols="euclidean", 
                       clustering_method="complete", border_color=FALSE)
  png(paste(c(fileprefix, ".hierarch.png"), sep="", collapse=""), width = 1600, height = 2000, pointsize = 20, res = 300)
  print(hierarch)
  dev.off()
  
  
  kmeans_res <- iterative_kmeans(input_data, k_val, 1, k_alg)

 
  # Perform first round of k-means clustering
  km <- pheatmap(kmeans_res$centers,
                 cluster_cols = F, cluster_rows = T)
  png(paste(c(fileprefix, ".kmeans.png"), sep="", collapse=""), width = 1600, height = 2000, pointsize = 20, res = 300)
  print(km)
  dev.off()
  # Visualise cluster distances
  png(paste(c(fileprefix, ".clustertree.png"), sep="", collapse=""), width = 2000, height = 1600, pointsize = 10, res = 300)
  plot(km$tree_row)
  abline(h=treeheight, col="red", lty=2, lwd=2)
  dev.off()
  # Generate merged clusters based on hierarchical clustering of centroids
  new_labels <- sort(cutree(km$tree_row, h=treeheight))
  km_names <- names(kmeans_res$cluster)
    new_cluster <- c()
  for(i in kmeans_res$cluster){
    label <- new_labels[which(names(new_labels) == i)]
    new_cluster <- c(new_cluster, label)
  }
  names(new_cluster) <- km_names
  cluster_table <- table(new_cluster)
  # Update the centroid values for input back into pheatmap
  new_centers <- matrix(data = NA, nrow = length(unique(new_labels)), ncol=ncol(kmeans_res$centers))
  for(i in unique(new_labels)){
    old_centers <- names(which(new_labels == i))
    if(length(old_centers)>1){
      new_centers[i,] <- colMeans(kmeans_res$centers[old_centers,])
    }else{
      new_centers[i,] <- kmeans_res$centers[old_centers,]
    }
  }
  # Update axis
  rownames(new_centers) <- lapply(unique(new_labels), function(x) paste(c("Cluster: ", 
                                                                          as.character(x), 
                                                                          " Size: ", 
                                                                          as.character(cluster_table[x])),
                                                                        collapse=""))
  colnames(new_centers) <- colnames(kmeans_res$centers)
  # Generate a merged heatmap
  km_u <- pheatmap(new_centers,
                   cluster_cols = F, cluster_rows = T)
  png(paste(c(fileprefix, ".kmeans_hierarch.png"), sep="", collapse=""), width = 1600, height = 2000, pointsize = 20, res = 300)
  print(km_u)
  dev.off()
  
  # Show new hierarchical clustering of the merged
  km_u <- pheatmap(new_centers,
                   cluster_cols = F, cluster_rows = T)
  png(paste(c(fileprefix, ".kmeans_hierarch.png"), sep="", collapse=""), width = 1600, height = 2000, pointsize = 20, res = 300)
  print(km_u)
  dev.off()
  return(new_cluster)
}

### GENERATE CLUSTER GRAPHS FOR MULTIDRUG ###
multi_graph <- function(input_data, fileprefix){
  # Convert wide data to long
  MULTIwide_data <-  as.data.frame(input_data[,1:4])
  colnames(MULTIwide_data) <- c(0.5, 1, 2, 4)
  MULTIwide_data <- data.frame(apply(MULTIwide_data, 2, function(x) as.numeric(as.character(x))))
  # Add transcript ID, cluster number and transcript type (lncrna/pcg)
  MULTIwide_data$transcript <- paste(multidrug_df$V5, multidrug_df$V6,  sep="_")
  MULTIwide_data$cluster <-  input_data[, 7]
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
  
  # Plot clusters
  p <- ggplot(data=MULTIlong_data, aes(x=condition, y=value, group=transcript)) +
    geom_line(color="#00BFC4") + 
    geom_line(data=MULTImedian_long, aes(x=condition, y=value, group=transcript), color="red") +
    geom_hline(yintercept=0) +
    facet_wrap(~ cluster) +
    xlab("MIC") + ylab("Normalized Log2Fold Change")
  
  png(paste(c(fileprefix, ".clustgraph.png"), sep="", collapse=""), width = 2000, height = 2000, pointsize = 20, res = 300)
  print(p)
  dev.off()
}

### ANALYSE GROUPINGS OF DRUG RESPONSES ###
multi_shared <- function(input_data, fileprefix){
  input_data <- multidrug_df
  # Attempt a Venn Diagram
  venn_data <- input_data
  # Combine the cluster ID with the transcript ID
  venn_data$x <- paste(venn_data[,6],venn_data[,7], sep="_")
  # Save each transcript per drug
  venn_list <- list()
  for(i in unique(venn_data[,5])){
    venn_list[[i]] <- venn_data$x[venn_data[,5] == i]
  }
  venn_out <- venn(venn_list)
  venn_overview <- as.data.frame(cbind(rownames(venn_out),venn_out$counts))[-1,]
  colnames(venn_overview) <- c("combination", "count")
  
  upset_p <- upset(fromList(venn_list), sets = names(venn_list),
                   #queries = list(list(query = intersects, params = list("simv", "terb"), color = "orange", active = T),
                   #               list(query = intersects, params = list("simv", "terb", "milt"), color = "green", active = T)),
                   point.size = 3.5, line.size = 2, 
                   order.by = "freq", text.scale = c(1.3, 1.3, 1, 1, 2, 1), 
                   mainbar.y.label = "Number of lncRNA", sets.x.label = "DE lncRNA Per Condition")
  
  # Look at what is shared between the drugs
  x1 <- unlist(venn_list, use.names = FALSE)
  x1 <- x1[ !duplicated(x1) ]
  # Generate an empty matrix for each combination
  cluster_frequency <- matrix(nrow = length(venn_overview$combination), ncol = length(unique(input_data[, 7])), data = 0)
  exact_cluster_frequency <- matrix(nrow = length(venn_overview$combination), ncol = length(unique(input_data[, 7])), data = 0)
  exact_combo_list <- list()
  combo_transcript_list <- list()
  for(combo_i in 1:length(venn_overview$combination)){
    combo <- venn_overview$combination[combo_i]
    single_combo_list <- list()
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
    exfreq_raw_cluster <- unlist(lapply(exact_transcripts, function(x) as.numeric(str_split(x, "_")[[1]][2])))
    exfreq_raw_name <- unlist(lapply(exact_transcripts, function(x) str_split(x, "_")[[1]][1]))
    excluster_ids <- names(exfreq_table)
    
    for(cluster_i in 1:length(excluster_ids)){
      cluster_transcripts <- exfreq_raw_name[which(exfreq_raw_cluster == excluster_ids[cluster_i])]
      single_combo_list[as.character(excluster_ids[cluster_i])] <- list(cluster_transcripts)
      exact_cluster_frequency[combo_i, as.numeric(excluster_ids[cluster_i])] <- as.numeric(exfreq_table[cluster_i])
      
    }
    if(length(single_combo_list)>0){
      combo_transcript_list[combo] <- list(single_combo_list)
    }
    
  }
  cluster_freq_df <- as.data.frame(cluster_frequency, row.names = venn_overview$combination)
  colnames(cluster_freq_df) <- as.character(c(1:length(unique(input_data[,7]))))
  cluster_freq_df$total <- rowSums(cluster_freq_df)
  
  excluster_freq_df <- as.data.frame(exact_cluster_frequency, row.names = venn_overview$combination)
  colnames(excluster_freq_df) <- as.character(c(1:length(unique(input_data[,7]))))
  excluster_freq_df$total <- rowSums(excluster_freq_df)
  
  # Save transcript data
  saveRDS(combo_transcript_list, file = paste(c(fileprefix,".clustered_lncrna.rds"), sep ="", collapse=""))
  
  # Save cluster frequency matrix
  write.csv(cluster_freq_df, paste(c(fileprefix,".cluster_freq.csv"), sep ="", collapse=""))
  write.csv(excluster_freq_df, paste(c(fileprefix,".exact_cluster_freq.csv"), sep ="", collapse=""))
  png(paste(c(fileprefix,".upset.png"), sep ="", collapse=""), width = 3000, height = 3000, pointsize = 20, res = 300)
  print(upset_p)
  dev.off()
  upset_p$x1 <- x1
  return(upset_p)
}


### GENERATE CLUSTER GRAPHS FOR ITRA ###
itra_graph <- function(input_data, inverse_clusters, fileprefix){
  # Convert wide data to long
  wide_data <-  as.data.frame(input_data[,-ncol(input_data)])
  # Add transcript ID, cluster number and transcript type (lncrna/pcg)
  wide_data$transcript <- rownames(wide_data)
  wide_data$cluster <- input_data[,ncol(input_data)]
  wide_data$inv_cluster <- inverse_clusters
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
  
  png(paste(c(fileprefix, ".clustgraph.png"), sep="", collapse=""), width = 2000, height = 2000, pointsize = 20, res = 300)
  print(p)
  dev.off()
  return(wide_data)
}

check_multi <- function(input_cluster_df, input_matrix){
  condition_combinations <- c()
  condition_numbers <- c()
  
  for(lncrna in input_cluster_df[,3]){
    if(lncrna %in% rownames(input_matrix)){
      lncrna_presence <- input_matrix[lncrna,]
      combo <- paste(colnames(lncrna_presence)[which(lncrna_presence == 1)], collapse=":")
      number_of_conditions <- unname(rowSums(lncrna_presence))
  
    }else{
      combo <- 0
      number_of_conditions <- 0
       }
   condition_combinations <- c(condition_combinations, combo)
    condition_numbers <- c(condition_numbers, number_of_conditions)
  }
  output_cluster_df <- cbind(input_cluster_df, condition_numbers, condition_combinations)
  return(output_cluster_df)
}

## IDENTIFY WHETHER TRANSCRIPTS ARE IN THE SAME CLUSTER
match_clusters <- function(input_data, cluster_data){
  # Generate a new dataframe with cluster columns
  input_data_w_cluster <- input_data
  input_data_w_cluster$t_clust <- rep(NA, nrow(input_data_w_cluster))
  input_data_w_cluster$n_clust <- rep(NA, nrow(input_data_w_cluster))
  input_data_w_cluster$inv_t_clust <- rep(NA, nrow(input_data_w_cluster))
  input_data_w_cluster$inv_n_clust <- rep(NA, nrow(input_data_w_cluster))
  
  for(i in 1:nrow(input_data)){
    # Extract transcript names
    t_name <- input_data[i,1]
    n_name <- input_data[i,3]
    
    # Extract relevant IDs in wide data
    t_wide_id <- which(t_name == cluster_data$transcript)
    n_wide_id <- which(n_name == cluster_data$transcript)
    # Check transcript exists in DESEQ results, if so add cluster details
    if(length(t_wide_id) > 0){
      t_data <- cluster_data[t_wide_id, ]
      input_data_w_cluster$t_clust[i] <- t_data$cluster
      input_data_w_cluster$inv_t_clust[i] <- t_data$inv_cluster
    }
    if(length(n_wide_id) > 0){
      n_data <- cluster_data[n_wide_id, ]
      input_data_w_cluster$n_clust[i] <- n_data$cluster
      input_data_w_cluster$inv_n_clust[i] <- n_data$inv_cluster
    }
    
  }
  
  
  same_clusters <- which(input_data_w_cluster$t_clust == input_data_w_cluster$n_clust)
  same_cluster_df <- input_data_w_cluster[same_clusters, ]
  # if(length(upset_data)>0){
  #   same_cluster_df <- check_multi(same_cluster_df, upset_data)
  # }
  # 
  inv_clusters <- which(input_data_w_cluster$inv_t_clust == input_data_w_cluster$inv_n_clust)
  inv_cluster_df <- input_data_w_cluster[inv_clusters, ]
  # if(length(upset_data)>0){
  #   inv_cluster_df <- check_multi(inv_cluster_df, upset_data)
  # }
  return(list(raw_cluster_data=input_data_w_cluster,
              same_clust_data=same_cluster_df,
              inv_clust=inv_cluster_df))
}


extract_multi_clust_info <- function(input_data, input_names){
  multidrug_cluster_df <- input_data
  rownames(multidrug_cluster_df) <- input_names
  multidrug_num <- as.character(unname(rowSums(multidrug_cluster_df)))
  multidrug_clust <- unlist(lapply(input_names, function(x) unlist(as.character(str_split(x, "_")[[1]][2]))))
  multidrug_names <- unlist(lapply(input_names, function(x) unlist(as.character(str_split(x, "_")[[1]][1]))))
  multidrug_combo <- c()
  for(i in 1:nrow(multidrug_cluster_df)){
    lncrna_presence <- multidrug_cluster_df[i, ]
    combo <- paste(colnames(lncrna_presence)[which(lncrna_presence == 1)], collapse=":")
    multidrug_combo <- c(multidrug_combo, combo)
  }
  multidrug_lncrna_df <- as.data.frame(cbind(multidrug_names, multidrug_clust, multidrug_num, multidrug_combo))
  return(multidrug_lncrna_df)
}


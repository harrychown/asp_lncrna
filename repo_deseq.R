### DESEQ2 SCRIPTS - HARRY CHOWN
#install.packages("DESeq2")
library("DESeq2")
subset_df <- function(input_df, input_rows, variable){
    if(length(variable) == 2){
    out_data  <- as.data.frame(input_df[input_rows,])
    rownames(out_data) <- rownames(input_df)[input_rows]
    colnames(out_data) <- variable[-1]
    } else {
    out_data <- input_df[input_rows,]
    colnames(out_data) <- variable[-1]}
    return(out_data)
}

# Clear environment variables prior to starting script
deseq_func <- function(count_directory, split_file_path, p_val, min_base_mean, min_LFC, signif_base_mean){
    # Read the split file containg information about input files and their organisation
  split_file <- read.csv(split_file_path, header = F)
  
  # First line of split file should contain condition labels
  condition_row <- split_file[1,]
  condition_row <- condition_row[!sapply(condition_row, function(x) all(x == ""))] # Removes any additional columns from the line
  
  conditions <- as.vector(unname(unlist(condition_row))) 
  number_of_conditions <- length(conditions)
  
  # Remainder of the file should be the filenames of each experiment and repetition
  input_file_names <- split_file[2:nrow(split_file),]
  input_file_names <- input_file_names[ , apply(input_file_names, 2, function(x) !any(is.na(x)))] # Removes any empty columns
  expected_reps <- ncol(input_file_names) 
  
  
  
  # Error catching to ensure there are the correct number of condtions
  if(nrow(input_file_names) != number_of_conditions){
    stop("Incorrect number of conditions")
  }
  
  sampleFiles = c()
  # Error catching to ensure there are the correct number of conditions
  for(line in 1:nrow(input_file_names)){
    file_name_row <- input_file_names[line,]
    file_name_row <- file_name_row[!sapply(file_name_row, function(x) all(x == ""))] # Removes any additional columns from the line
    file_name_row <- as.vector(unname(unlist(file_name_row)))
    if(length(file_name_row) == expected_reps){
      sampleFiles <- c(sampleFiles, file_name_row)
    } else {
      sampleFiles <- c()
      stop("Incorrect number of files/repetitions")
    }
  }
  
  # Generate condition labels for each file
  sampleCondition <- c()
  for(c in conditions){
    condition_reps <- rep(c, expected_reps)
    sampleCondition <- c(sampleCondition, condition_reps)
  } 
  
  # Data frame used to combine all deseq output data
  comb_data <- data.frame()
  
  setwd(count_directory)
  directory <- setwd(count_directory)
  
  
  
  sampleTable <- data.frame(sampleName = sampleFiles,
                            fileName = sampleFiles,
                            condition = sampleCondition)
  # Prepare data for DESeq
  
  
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                         directory = directory,
                                         design= ~ condition)
  
  
  ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = conditions[1])
  # Run DESeq and gather results
  dds <- DESeq(ddsHTSeq)
  # contrast and name is used to specify the exact contrast we want to build
  dds$condition <- relevel(dds$condition, ref = conditions[1])
  
  # List to contain all the row names
  rname_list = list()
  # List to contain all p-values for each gene
  pval_list = list()
  # List to contain all adjusted p-values for each gene
  adjpval_list = list()
  # List to contain all shrink log 2 fold changes for each gene
  shrinklog_list = list()
  # List to contain all log 2 fold changes for each gene
  log2fold_list = list()
  # List to contain all stats for each gene
  stat_list = list()
  # List to contain all lfcSE for each gene
  lfcSE_list = list()
  # Get the basemean for each transcript in each condition
  bmean_df <- sapply( levels(dds$condition), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition == lvl, drop=F] ) )
  bmean_df <- bmean_df[,-1]
  # Dataframe for each dds result
  dds_list = list()
  # Iterate through the results of pairwise comparisons and store each pval, bmean and shrinklog
  for(i in 1:length(conditions[-1])) {
    condition_name <- c("condition", conditions[-1][i], conditions[1])
    dds_results <- results(dds, condition_name)
    dds_list[[conditions[-1][i]]] <- dds_results
    pval_list[[conditions[-1][i]]] <- dds_results$padj
    log2fold_list[[conditions[-1][i]]] <- dds_results$log2FoldChange
    rname_list = list(rname = rownames(dds_results))
  }
  
  pval_df <- as.data.frame(pval_list)
  log2fold_df <- as.data.frame(log2fold_list)
  rownames(pval_df) = rownames(bmean_df) = rownames(log2fold_df) = rname_list$rname
  
  # Store gene row numbers that contain some significance
  row_numbers = c()
  # Iterate throught the dataframe of p-values and baseMeans to find rows which satisfy the
  # following conditions: p-value <0.05 or bmean > 5000, at least once, and the gene must have a logfold change greater/less than +1/-1, and have a basemean read count greater than 30
  for(i in 1:nrow(pval_df)){
    # Obtain the result for each gene
    # Identify whether a dosage obtained a significant change
    p_row = pval_df[i,]
  
    p_result = which(p_row < p_val)
     # Identify whether the gene has a high expression level
    b_row = bmean_df[i,]
  
    b_result = which(b_row > min_base_mean)
    b_high_result = which(b_row > signif_base_mean)
    
    # Identify whether the gene has a high logfold change
    l_row = log2fold_df[i,]
    l_result = which(l_row ** 2 > min_LFC ** 2)
    
  
    if((length(p_result) != 0 || length(b_high_result) != 0) && (length(l_result) != 0 ) && (length(b_result)!=0)){
      row_numbers <- c(row_numbers, i)
    }
  }
  

  signif_log <- subset_df(input_df = log2fold_df, input_rows = row_numbers, variable = conditions)
  signif_base <- subset_df(input_df = bmean_df, input_rows = row_numbers, variable = conditions)
  return(list(signif_log, signif_base, bmean_df, ddsHTSeq@assays@data@listData$counts))
  
}




  








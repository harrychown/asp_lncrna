# Perform analysis on neighbour pairs
library(dplyr)
analyse_neighbours <- function(output_prefix, out_folder, neighbour_df, signif_logchange_df, counts_df, p_val = 0.05, cutoff = 0.5, bmean_df){
  subset_dataframe <- function(input_df, subset_vector, index){
    subset_index <- c()
    for(i in 1:nrow(input_df)){
      if(input_df[i,index] %in% subset_vector){
        subset_index <- c(subset_index, i)
      }
    }
    output <- input_df[subset_index,]
    return(output)
  }
  
  proportional_regulation <- function(input_df, sub_prefix){
   prop_reg_df <-  input_df[which(input_df$r_val.cor>0),]
   positive_proportional <- prop_reg_df[which(prop_reg_df$lncRNA_r_val.cor>0),]
   negative_proportional <- prop_reg_df[which(prop_reg_df$lncRNA_r_val.cor<0),]
   
   write.csv(positive_proportional,paste(c(out_folder, output_prefix,"_", sub_prefix, "_pos_prop.csv"),sep="", collapse=""), row.names = T)
   write.csv(negative_proportional,paste(c(out_folder, output_prefix,"_", sub_prefix, "_neg_prop.csv"),sep="", collapse=""), row.names = T)
  }
  
  inverse_regulation <- function(input_df, sub_prefix){
    inverse_reg_df <-  input_df[which(input_df$r_val.cor<0),]
    positive_inverse <- inverse_reg_df[which(inverse_reg_df$lncRNA_r_val.cor>0),]
    negative_inverse <- inverse_reg_df[which(inverse_reg_df$lncRNA_r_val.cor<0),]
    
    write.csv(positive_inverse,paste(c(out_folder, output_prefix,"_", sub_prefix, "_pos_inv.csv"),sep="", collapse=""), row.names = T)
    write.csv(negative_inverse,paste(c(out_folder, output_prefix,"_", sub_prefix, "_neg_inv.csv"),sep="", collapse=""), row.names = T)
  }
  
  raw_count_correlation <- function(n_df, c_df, b_df){
    # Transpose count/mean data so that the rows = conditions/replicates and cols = transcripts
    tc_df <- t(c_df)
    tb_df <- t(b_df)
    # Iterate through the input data and perform correlation tests using the raw count/base mean data
    correlation_results_table <- list()
    for(i in 1:nrow(n_df)){
      lncRNA_id <- n_df[i,1]
      PCG_id <- n_df[i, 3]
      
      
      # Check if there is correlation between expression of transcripts
      cor_res <- cor.test(tc_df[, lncRNA_id], tc_df[, PCG_id])
      # Check if there is correlation between dosage
      y <- as.numeric(colnames(b_df))
      lncRNA_cor <- cor.test(tb_df[, lncRNA_id], y)
      PCG_cor <- cor.test(tb_df[, PCG_id], y)
  
  
      correlation_results_table[[i]] <- c("r_val" = cor_res$estimate, "sqrd_r" = cor_res$estimate ** 2, "p" = cor_res$p.value, 
                                          "lncRNA_r_val" = lncRNA_cor$estimate, "lncRNA_sqrd_r" = lncRNA_cor$estimate ** 2, "lncRNA_p" = lncRNA_cor$p.value,
                                          "PCG_r_val" = PCG_cor$estimate, "PCG_sqrd_r" = PCG_cor$estimate ** 2, "PCG_p" = PCG_cor$p.value)
    }
    # Convert results to a dataframe
    correlation_results_df <- as.data.frame(do.call(rbind, correlation_results_table))
    # Combine with the original neighbour input
    out_df <- cbind(n_df, correlation_results_df)
    filtered_df <- filter(out_df, sqrd_r.cor > cutoff  & p < p_val)
    return(filtered_df)
  }
    
  # Subset based on correlation
  # Subset the neighbour dataset based on significant DE lncRNA
  signif_lncRNA <- rownames(signif_logchange_df)[grep("MSTRG", rownames(signif_logchange_df))]

  signif_PCG <- rownames(signif_logchange_df)[grep("MSTRG", rownames(signif_logchange_df), invert = T)]
 
  subset_by_lncRNA <- subset_dataframe(input_df = neighbour_df, subset_vector = signif_lncRNA, index = 1)
  print("LNCRNA complete")
  print(nrow(subset_by_lncRNA))
  subset_by_pcg <- subset_dataframe(input_df = neighbour_df, subset_vector = signif_PCG, index = 3)
  print("PCG complete")
  print(nrow(subset_by_pcg))
  print(subset_by_pcg)
  
  subset_by_lncRNA_and_PCG <- subset_dataframe(input_df = subset_by_lncRNA, subset_vector = signif_PCG, index = 3)
  print("both complete")
  
  
  cutoff <- cutoff ** 2
  # Subsets based on DESeq results
  # signif_cor contains results where both neighbours are significantly DE
  none_neighbours <- raw_count_correlation(neighbour_df, counts_df, bmean_df)
  proportional_regulation(none_neighbours, "none")
  inverse_regulation(none_neighbours, "none")
  
  both_neighbours <- raw_count_correlation(subset_by_lncRNA_and_PCG, counts_df, bmean_df)
  proportional_regulation(both_neighbours, "both")
  inverse_regulation(both_neighbours, "both")
  
  lnc_neighbours <- raw_count_correlation(subset_by_lncRNA, counts_df, bmean_df)
  proportional_regulation(lnc_neighbours, "lnc")
  inverse_regulation(lnc_neighbours, "lnc")
  
  pcg_neighbours <- raw_count_correlation(subset_by_pcg, counts_df, bmean_df)
  proportional_regulation(pcg_neighbours, "pcg")
  inverse_regulation(pcg_neighbours, "pcg")
  
  
  
  return(list(both_neighbours, none_neighbours, pcg_neighbours, lnc_neighbours))
}
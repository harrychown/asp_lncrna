# Identifying genes neighbouring lncRNA within 5,000 bp
library("biomaRt")

# Function used to annotate PCG's in the neighbour dataframe
annotate_pcg <- function(input_df){
  # Biomart for ensembl fungi
  ensembl_fungi <- useMart(host="https://fungi.ensembl.org", 
                           biomart="fungi_mart", 
                           port = 443, dataset = "afumigatusa1163_eg_gene")
  # Obtain the attributes of the sense genes
  neigh_attr <- getBM(attributes=c('ensembl_gene_id', 'description', 'ensembl_transcript_id'),
                      filters = "ensembl_gene_id",
                      values = input_df[,2], 
                      mart = ensembl_fungi) 
  neigh_attr <- apply(neigh_attr, 2, function(x) gsub("^$|^ $", NA, x))
  
  # Order attribute results to match with the input dataframe
  ordered_desc <- c()
  for(afub in input_df[,2]){
    desc <- neigh_attr[which(neigh_attr[,1] == afub),2]
    ordered_desc <- c(ordered_desc, desc)
    
  }
  
  out_attr <- cbind(input_df, ordered_desc)
  return(out_attr)
}

get_neighbours <- function(input_gtf, lncRNA_names, PCG_names, distance_threshold, max_neighbours){
  # Generate BED files for lncRNA and PCG
  generate_bed <- function(gtf_vector){
    bed_tmp <- tempfile(pattern = "bed_tmp")
    bed_cmd <- paste(c("gtf2bed | sort -k1,1 -k2,2n > ",bed_tmp), sep = "", collapse = "")
    system(bed_cmd, input = gtf_vector, intern = FALSE)
    return(bed_tmp)
  }
  
  # Function used to identify lncRNAs with neighbours within 5,000 bp
  identify_neighbours <- function(lncRNA_bed_input, PCG_bed_input, threshold = 5000, num_neighbours = 10){
    closest_cmd <- paste(c("bedtools closest -io -k ",as.character(num_neighbours)," -t all -d -a ", lncRNA_bed_input, " -b ", PCG_bed_input), 
                         sep = "", collapse ="")
    closest_results <- system(closest_cmd, intern = TRUE)
    gene_t_list <- c()
    gene_id_list <- c()
    lncRNA_t_list <- c()
    distance_list <- c()
    for(r in closest_results){
      # Split the line
      split <- strsplit(r, "\t")
      # Check if they are classed as neighbours
      distance <-  as.integer(split[[1]][21])
      if(distance < threshold){
        lncRNA_desc <- split[[1]][10]
        gene_desc <- split[[1]][20]
        lncRNA_transcript <- str_match(lncRNA_desc, "transcript_id \"(.*?)\"")[1,2]
        gene_transcript <- str_match(gene_desc, "transcript_id \"(.*?)\"")[1,2]
        gene_id <- str_match(gene_desc, "gene_id \"(.*?)\"")[1,2]
        if(is.na(gene_id)){
          next()
        }
        gene_t_list <- c(gene_t_list, gene_transcript)
        gene_id_list <- c(gene_id_list, gene_id)
        lncRNA_t_list <- c(lncRNA_t_list, lncRNA_transcript)
        distance_list <- c(distance_list, distance)
        
      }
      
    }
    neighbour_out <- do.call(rbind, Map(data.frame, lncRNA=lncRNA_t_list, PCG_id=gene_id_list, PCG_transcript=gene_t_list, distance=distance_list))
    return(neighbour_out)
  }
  
  # Function used to annotate PCG's in the neighbour dataframe
  annotate_pcg <- function(input_df){
    # Biomart for ensembl fungi
    ensembl_fungi <- useMart(host="https://fungi.ensembl.org", 
                             biomart="fungi_mart", 
                             port = 443, dataset = "afumigatusa1163_eg_gene")
    # Obtain the attributes of the sense genes
    neigh_attr <- getBM(attributes=c('ensembl_gene_id', 'description', 'ensembl_transcript_id'),
                        filters = "ensembl_gene_id",
                        values = input_df[,2], 
                        mart = ensembl_fungi) 
    neigh_attr <- apply(neigh_attr, 2, function(x) gsub("^$|^ $", NA, x))
    
    # Order attribute results to match with the input dataframe
    ordered_desc <- c()
    for(afub in input_df[,2]){
      desc <- neigh_attr[which(neigh_attr[,1] == afub),2]
      ordered_desc <- c(ordered_desc, desc)

    }
   
    out_attr <- cbind(input_df, ordered_desc)
    return(out_attr)
  }
  
  lncRNA_bed_path <- generate_bed(lncRNA_names)
  PCG_bed_path <- generate_bed(PCG_names)
  neighbours_df <- identify_neighbours(lncRNA_bed_path, PCG_bed_path, threshold = distance_threshold, num_neighbours = max_neighbours)
  gene_attributes <- annotate_pcg(neighbours_df)
  
  return(gene_attributes)
}
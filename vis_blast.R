library(pheatmap)
raw <- read.csv("/home/harry/github/asp_lncrna/blast_res/latest_analysis/blast_binary.header.csv", row.names = 1)
blast_de <- read.table("/home/harry/github/asp_lncrna/blast_res/latest_analysis/blast_de.txt")
de_index <- which(colnames(raw) %in% blast_de[,1])
de_binary <- rep(0, (ncol(raw) -1))
de_binary[de_index] <- 1
de_binary_df <- t(as.data.frame(de_binary))
t_raw <- t(raw)
rownames(raw) <- gsub("@", " ", rownames(raw))
library(pheatmap)
pheatmap(raw,
         color = colorRampPalette(rev(c("#238210", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#FFFFFF")))(100),
         cellwidth = 4, cellheight = 20,   # changed to 3
         border_color="black",
         treeheight_row=0,
         treeheight_column=0,
         kmeans_k = NA,
         show_rownames = T, show_colnames = F,
         fontsize=15,
         scale="none",
         clustering_method = "complete",
         cluster_rows = F, cluster_cols = F,
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         legend=F
)



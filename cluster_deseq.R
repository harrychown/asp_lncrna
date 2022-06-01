# Perform hierarchical and K-means clustering on heatmap data
library(ComplexHeatmap)
cluster_fun <- function(input, x_axis_title = "Drug doasge", y_axis_title = "Genomic loci", k_val = 8){
  out_og <- Heatmap(input, 
          name = "LFC", #title of legend
          column_title = x_axis_title, row_title = y_axis_title,
          show_row_names = F
  )
  
  out_k <- Heatmap(input, 
          name = "LFC", #title of legend
          column_title = x_axis_title, row_title = y_axis_title,
          show_row_names = F,
          k = k_val
  )
  return(c(out_og, out_k))
}
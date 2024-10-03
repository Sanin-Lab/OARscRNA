#' Calculate MSE Distance
#'
#' @param x uncomplete
#' @param y uncomplete
#'
#' @return uncomplete
#' @export
#'
#' @examples 
#' \dontrun {
#' uncomplete
#' }
mse_distance <- function(x, y) {
  mean((x - y)^2)
}

#' Run MSE Distance on Seurat Object 
#'
#' @param x A Seurat Object 
#'
#' @return uncomplete
#' @export
#'
#' @examples 
#' \dontrun {
#' uncomplete
#' }
SeuratDistance <- function(x){

  # Setup matrix
  data_matrix <- GetAssayData(x, assay = "RNA", layer = "data")
  
  clustAvg <- averageColumns(
    GetAssayData(x, assay = "RNA", layer = "data"),
    by = Idents(x))
  
  distances <- future_lapply(colnames(clustAvg), function(cluster_a) {
    sapply(colnames(clustAvg), function(cluster_b) {
      mse_distance(clustAvg[ , cluster_a], clustAvg[ , cluster_b])
    })
  })
  
  # Setup idents
  cluster_ids <- unique(Idents(x))
  
  # Compute Distance
  system.time({
    distances <- future_lapply(cluster_ids, function(cluster_a) {
      sapply(cluster_ids, function(cluster_b) {
        data_a <- data_matrix[, Idents(x) == cluster_a]
        data_b <- data_matrix[, Idents(x) == cluster_b]
        # Computing MSE between cluster centroids as an example
        mse_distance(rowMeans(data_a), rowMeans(data_b))
      })
    })
  })
  
  ### Transform the list of lists into a single vector
  distance_vector <- unlist(distances)
  
  # Now, convert this vector into a matrix. The number of rows and columns should match the length of cluster_ids
  distance_matrix <- matrix(distance_vector, nrow = length(cluster_ids), ncol = length(cluster_ids), byrow = TRUE)
  
  # To make the matrix more readable, you can name the rows and columns as the cluster IDs
  rownames(distance_matrix) <- cluster_ids
  colnames(distance_matrix) <- cluster_ids
  
  # Convert the matrix to a long format
  distance_df <- as.data.frame(as.table(distance_matrix))
  names(distance_df) <- c("cluster_a", "cluster_b", "distance")
  return(distance_df)
}
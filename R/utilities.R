#' Prepare data for oar fold functions
#'
#' @param data Seurat object or gene expression matrix
#' @param tr Filtering threshold
#' @param seurat_v5 A Boolean to indicate if supplied data is a Seurat object, default is TRUE
#' @param blacklisted.genes A character vector with gene names to be excluded from the analysis. Default is empty.
#'
#' @return Data matrix with blacklisted gens removed
#' @export
#'
#' @examples 
#' \dontrun{
#' data <- oar_preprocess_data(data)
#' }
oar_preprocess_data <- function(data, tr, seurat_v5 = T, blacklisted.genes = NULL) {
  #read in Seurat object
  if(seurat_v5){
    layersList <- lapply(data@assays$RNA@layers,function(x){dim(x)}) #Identify empty layers
    data@assays$RNA@layers[names(layersList[sapply(layersList, is.null)])] <- NULL #remove empty layers
    
    if(length(layersList) >3){
      print("Data may contain multiple raw count layers")
      print("Joining data layers (time consuming!!) and extracting count tables")
    }else{
      print("Extracting count tables")}
    
    data <- SeuratObject::JoinLayers(object = data, assay = "RNA") %>%
      SeuratObject::LayerData(assay = "RNA", layer = "counts")
    
    data <- data[Matrix::rowSums(data > 0) > data@Dim[[2]]*tr,]
  }else{
    data <- data %>%
      methods::as("dgCMatrix")
    data <- data[Matrix::rowSums(data > 0) > data@Dim[[2]]*tr,]
  }
  
  #Remove blacklisted genes
  if(length(blacklisted.genes) > 0 && is.character(blacklisted.genes)){
    data <- data[!rownames(data) %chin% blacklisted.genes,]
  }
  return(data)
}

#' Identify number of folds to use for data based on desired gene-cell ratio
#'
#' @param data A matrix with low expression genes filtered
#' @param gene.ratio A numeric value indicating the ratio of genes to cells to partition data matrix. Default is 20. Values between 15 and 25 are recommended.
#' @param max.cells Total number of cells available in data
#'
#' @return A variable called folds that contains a numeric value
#' @export
#'
#' @examples
#' \dontrun{
#' folds <- oar_identify_folds(data)
#' }
oar_identify_folds <- function(data, gene.ratio, max.cells) {
  print("Identifying effective ratio:")
  if(max.cells >= data@Dim[[2]]){
    folds = 1
    max.cells = data@Dim[[2]]
  }else{
    n=1
    while(!data@Dim[[2]] %% max.cells == 0){
      max.cells = ceiling(ceiling(data@Dim[[1]]/(gene.ratio*runif(n = 1,0.85,1.15))))
      if(
        dplyr::between(
          x = data@Dim[[2]] %% max.cells,
          left = max.cells*0.80,
          right = max.cells*1.2)) break
      n = n+1
      if(n == 50){
        gene.ratio = gene.ratio+2
      }
      if(n > 100){
        stop("Maximum number of iterations reached. Choose a different starting gene.ratio\n")
      }
    }
    folds = ceiling(data@Dim[[2]]/max.cells)
  }
  
  print(paste0("Effective ratio: ",ceiling(data@Dim[[1]])/max.cells))
  
  if (ceiling(data@Dim[[1]])/max.cells < 15) {
    warning("A small gene to cell ratio may result in unreliable outputs\n")
  }
  
  if (ceiling(data@Dim[[1]])/max.cells > 30) {
    warning("A large gene to cell ratio may result in unreliable outputs\n")
  }
  
  return(folds)
}


  
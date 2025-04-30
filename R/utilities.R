##===================================================================
#Preprocess data
##===================================================================
#' Prepare data for oar fold functions
#'
#' @param data a seurat object or gene expression matrix
#' @param tr a filtering threshold. Default is 1.
#' @param seurat_v5 a boolean to indicate if supplied data is a Seurat object, default is TRUE
#' @param blacklisted.genes a character vector with gene names to be excluded from the analysis. Default is empty.
#'
#' @return Data matrix with blacklisted genes removed
#' @export
#'
#' @examples 
#' \dontrun{
#' data <- oar_preprocess_data(data)
#' }
oar_preprocess_data <- function(data, tr = 1, seurat_v5 = TRUE, blacklisted.genes = NULL) {
  #Set filtering threshold
  tr = tr/100
  
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
    data <- Matrix::Matrix(data, sparse = T)
    data <- data[Matrix::rowSums(data > 0) > data@Dim[[2]]*tr,]
  }
  
  #Remove blacklisted genes
  if(!is.null(blacklisted.genes)){
    data <- data[!rownames(data) %chin% blacklisted.genes,]
  }
  
  #save gene names
  gene_names <- rownames(data)
  
  #Convert to a dense matrix and Replace 0 with NA
  data <- data %>% as.matrix()
  colnames(data) = NULL
  rownames(data) = NULL
  data[data == 0] <- NA
  
  # data must be in .data.frame()
  if (all(complete.cases(data))) {
    stop("No missing data exists\n")
  }
  output <- list(data, gene_names)
  
  return(output)
  
}

##===================================================================
#Group missing data patterns based on tolerance with a graph
##===================================================================
#' Group missing data patterns based on tolerance with a graph
#'
#' @param dm a matrix of gene vector hamming distances.
#' @param tol a numeric value indicating the maximum fraction of mismatch in genes to group as a pattern, default is 0.05.
#'
#' @return missing data pattern vector
#' @export
#'
#' @examples 
#' \dontrun{
#' mdp <- oar_missing_data_patterns(dm, tol = 0.05)
#' }
oar_missing_data_graph <- function (dm, tol = 0.05) {
  g <- igraph::graph_from_adjacency_matrix(
    adjmatrix = dm <= tol, mode = "undirected", diag = F) # Create an graph based on similar gene patterns
  g <- igraph::decompose(g) # Split into connected nodes
  mdp <- rep(NA,nrow(dm)) # create space holder for patterns
  ps <- 1:length(g) # create pattern numbers
  n=1
  for(i in g){ # Assign pattern numbers
    mdp[as.numeric(igraph::V(i)$name)] <- ps[n]
    n=1+n
  }
  return(mdp)
}

##===================================================================
#Test gene distributions across patterns
##===================================================================
#' Kruskal-Wallis test to generate a per cell p-value based on missing data patterns
#'
#' @param x Item from list of cell gene expression vectors
#' @param mdp Matrix with gene participation per pattern
#'
#' @return list with a p-value for each cell
#' @export
#'
#' @examples
#' \dontrun{
#' pvalue <- missing_pattern_pval_kw(x, mdp)
#' }
#' 
missing_pattern_pval_kw = function(x, mdp){
  y.l <- x[!is.na(x)] # subset observed genes of the nth cell to y.l
  mdp.l <- mdp[!is.na(x)] # subset observed genes of the nth cell to y.l
  pval = kruskal.test(x = y.l, g = factor(mdp.l))$p.value
  
  return(pval)
}


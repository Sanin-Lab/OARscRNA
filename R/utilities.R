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
    stop("Data is not sparse\n")
  }
  output <- list(data, gene_names)
  
  return(output)
  
}

##===================================================================
#Group gene co-expression patterns based on tolerance with a graph
##===================================================================
#' Group gene co-expression patterns based on tolerance with a graph
#'
#' @param dm a matrix of gene vector hamming distances.
#' @param tol a numeric value indicating the maximum fraction of mismatch in genes to group as a pattern.
#'
#' @return gene co-expression pattern vector
#' @export
#'
#' @examples 
#' \dontrun{
#' cgp <- oar_gene_graph(dm, tol)
#' }
oar_gene_graph <- function (dm, tol) {
  loop = TRUE
  while(loop){
    g <- igraph::graph_from_adjacency_matrix(
      adjmatrix = dm <= tol, mode = "undirected", diag = F)
    ecc <- max(igraph::eccentricity(g, mode = "all"))
    loop = ecc > 3
    if(loop){tol = tol - 0.01}
  }
  
  g <- igraph::decompose(g)
  out <- rep(NA, nrow(dm))
  ps <- 1:length(g)
  n = 1
  for (i in g) {
    out[as.numeric(igraph::V(i)$name)] <- ps[n]
    n = 1 + n
  }
  return(out)
}

##===================================================================
#Test gene distributions across patterns
##===================================================================
#' Kruskal-Wallis test to generate a per cell p-value based on gene co-expression patterns
#'
#' @param x Item from list of cell gene expression vectors
#' @param cgp Matrix with gene participation per pattern
#'
#' @return list with a p-value for each cell
#' @export
#'
#' @examples
#' \dontrun{
#' pvalue <- pattern_pval_kw(x, cgp)
#' }
#' 
pattern_pval_kw = function(x, cgp){
  y.l <- x[!is.na(x)] # subset observed genes of the nth cell to y.l
  cgp.l <- cgp[!is.na(x)] # subset observed genes of the nth cell to y.l
  if(length(unique(cgp.l)) > 1){ # check that cell participates in more than one pattern
    pval = kruskal.test(x = y.l, g = factor(cgp.l))$p.value
  }else{
    pval = NA # return NA for all cells that do not participate in more than one pattern
  }
  
  return(pval)
}


##===================================================================
#Base Test
##===================================================================
#' Generate scores and p-values to determine heterogeneity of data by looking at whether missingness is observed-at-random (OAR)
#'
#' @param data A gene-cell expression matrix with NA values in place of 0s.
#' @param mdp A vector indicating the pattern to which each gene belongs. 
#'
#' @return Data frame with OAR-score, p-value, adjusted p-value, and percent missing data for each cell. 
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' output <- oar_base(data, mdp)
#' }
oar_base <- function (data, mdp) {
  cl = lapply(seq_len(dim(data)[2L]), function(i) data[,i]) #convert mtx to list of cell vectors
  pvalue.list.KW <- unlist(lapply(cl, FUN = missing_pattern_pval_kw, mdp = mdp)) # Calculate p value
  pvalue.KW.BH = p.adjust(pvalue.list.KW, method = "BH") # Benjamini & Hochberg correction ("BH" or its alias "fdr")
  
  #transform and scale adjusted pvalues
  score = scale(-log10(pvalue.KW.BH), center = T, scale = T)*-1
  
  #calculate missing values
  sp = unlist(lapply(cl, naniar::pct_miss))
  
  #output
  out <- data.frame(score,pvalue.list.KW,pvalue.KW.BH,sp)
  colnames(out) <- c("OARscore","KW.pvalue","KW.BH.pvalue","pct.missing")
  return(out)
}

##===================================================================
#Wrapper function to run full OAR scoring
##===================================================================
#' Single line pipeline to run complete analysis
#'
#' @param data a Seurat (v5) object or a matrix with cell barcodes as column names and genes as row names.
#' @param seurat_v5 a boolean to indicate if supplied data is a Seurat object, default is TRUE
#' @param count.filter a numeric value indicating the minimum fraction of cells expressing any given gene that will be included in the analysis, default is 1. Values between 0.5 and 2 are recommended.
#' @param blacklisted.genes a character vector with gene names to be excluded from the analysis. Default is empty.
#' @param suffix a string to append to the output variables. Default is empty
#' @param tolerance a boolean or numeric value controlling the tolerance threshold for pattern matching. If set to `TRUE`, then tolerance is automatically adjusted. Alternatively, user may supply a numeric value indicating the maximum fraction of mismatch between pairs of genes for pattern grouping. Values between 0.01 and 0.05 are recommended.
#' @param cores a numeric value indicating the number of cores to use un parallel processing. Use `parallel::detectCores()` or `parallelly::availableCores()` to identify possibilities. Default is 1.
#' @param store.hamming a boolean to control if hamming distances should be stored in the Seurat object. Default set to `TRUE`
#'
#' @return A Seurat object with OAR stats added into meta data, or a matrix with OAR stats. 
#' @export
#'
#' @examples
#' \dontrun{
#' pbmcs <- oar(pbmcs)
#' }
#' 
oar <- function (data, seurat_v5 = TRUE, count.filter = 1, 
                 blacklisted.genes = NULL, suffix = "",
                 tolerance = TRUE, cores = 1, store.hamming = TRUE) {
  
  #Check parameters were correctly supplied
  if(!is.numeric(count.filter)){stop("count.filter must be numeric\n")}
  if(!is.character(suffix)){stop("suffix must be a string\n")}
  if(!is.logical(tolerance)){if(!is.numeric(tolerance))stop("tolerance must be logical or numeric\n")}
  if(!is.logical(seurat_v5)){stop("seurat_v5 must be TRUE or FALSE\n")}
  if(!is.null(blacklisted.genes)){if(!is.character(blacklisted.genes))stop("Supplied blacklisted genes are not characters\n")}
  if(seurat_v5){if(!is.logical(store.hamming)){stop("store.hamming must be TRUE or FALSE\n")}}
  
  #Value range warnings
  if (count.filter > 2) {
    warning("Overfiltering expression matrix (count.filter > 2) may lower signal detection\n")
  }else if(count.filter == 0) {
    warning("Minimum fraction of cells expressing any given gene should be greater than 0\n")
  }
  if(is.numeric(tolerance)){
    if (tolerance == 0) {
      warning("Tolerance should be greater than 0\n")
    }else if(tolerance >= 0.1) {
      warning("Tolerances greater than 0.10 are not recommended\n")
    }
  }
  
  #Processing Warnings 
  if (parallelly::availableCores() < cores ){
    stop("Specified cores are greater than available in your machine.\nRun `parallelly::availableCores()` to identify appropriate number\n")}
  if (cores < 2 ) {
    warning("Running process in fewer than 2 cores will considerably slow down progress\n")}
  
  #store Seurat object for later
  if(seurat_v5){
    sc.data = data
    cells = colnames(sc.data)
  }else{
    cells = colnames(data)
  }
  
  #read in Seurat object & remove blacklisted genes
  print("Extracting data...")
  output <- oar_preprocess_data(data, tr = count.filter, seurat_v5, blacklisted.genes)
  data <- output[[1]]
  gene_names <- output[[2]]
  
  ot <- Sys.time()
  print("Analysis started on:")
  print(ot)
  
  #Identify missing data patterns
  print("Identifying missing data patterns...")
  if(seurat_v5){
    if("hamming" %in% names(sc.data[["RNA"]]@misc)){
      dm <- sc.data[["RNA"]]@misc[["hamming"]]
    }else{
      dm <- oar_hamming_distance(data, cores = cores)
      if(store.hamming){
        sc.data[["RNA"]]@misc[["hamming"]] <- dm
      }
    }
  }else{
    dm <- oar_hamming_distance(data, cores = cores)
  }
  mdp <- oar_missing_data_patterns(data = data, dm = dm, tolerance = tolerance)
  
  #Run missingness test
  print("Calculating scores")
  output <- oar_base(data, mdp)
  output$barcodes = cells
  
  #add names to mdp
  names(mdp) <- gene_names
  
  # Clean up output
  print("Collecting results")
  colnames(output)[-ncol(output)] <- paste0(colnames(output)[-ncol(output)],suffix)
  
  print("Analysis completed succesully!")
  print(round((Sys.time()-ot),2))
  
  #output
  if(seurat_v5){
    idx <- match(rownames(sc.data@meta.data),output$barcodes)
    
    sc.data@meta.data <- base::cbind(
      sc.data@meta.data,
      output[idx,!colnames(output) %in% "barcodes"])
    sc.data[["RNA"]]<- SeuratObject::AddMetaData(sc.data[["RNA"]], mdp, col.name = "mdp")
    
    return(sc.data)
  }else{
    return(output)
  }
}

##===================================================================
#Oar Test By Cluster
##===================================================================
#' Generate OAR score within each cluster and add them to full objects metadata
#'
#' @param data a clustered Seurat (v5) object.
#' @param count.filter a numeric value indicating the minimum fraction of cells expressing any given gene that will be included in the analysis, default is 1. Values between 0.5 and 2 are recommended.
#' @param blacklisted.genes a character vector with gene names to be excluded from the analysis. Default is empty.
#' @param suffix a string to append to the output variables. Default is empty
#' @param tolerance a boolean or numeric value controlling the tolerance threshold for pattern matching. If set to `TRUE`, then tolerance is automatically adjusted. Alternatively, user may supply a numeric value indicating the maximum fraction of mismatch between pairs of genes for pattern grouping. Values between 0.01 and 0.05 are recommended.
#' @param cores a numeric value indicating the number of cores to use un parallel processing. Use `parallel::detectCores()` or `parallelly::availableCores()` to identify possibilities.
#'
#' @return a Seurat object with OAR stats added into meta data
#' @export
#'
#' @examples
#' \dontrun{
#' pmbcs <- oar_by_cluster(pmbcs)
#' }
oar_by_cluster <- function (data, count.filter = 1,
                            blacklisted.genes = NULL, suffix = "",
                            tolerance = TRUE,
                            cores = 1) {
  
  print("Splitting data by cluster...")
  
  data_list <- Seurat::SplitObject(data) #split into list of objects by active ident
  if("hamming" %in% names(data[["RNA"]]@misc)){
    data_list <- lapply(data_list,function(x){
      x[["RNA"]]@misc <- list()
      return(x)})
    warning("Previously calculated hamming distances will not be used/stored when estimating OAR score by cluster")
  }
  
  #run oar on each object
  data_oar <- lapply(data_list, oar, count.filter = count.filter,
                     blacklisted.genes = blacklisted.genes, 
                     suffix = suffix, store.hamming = FALSE,
                     tolerance = tolerance, cores = cores)
  
  #combine objects back together
  collect <- paste0(c("OARscore","KW.pvalue","KW.BH.pvalue","pct.missing"),suffix)
  oar_combine <- lapply(data_oar, function(x){x@meta.data %>% dplyr::select(dplyr::all_of(collect))})
  oar_combine <- do.call(rbind,oar_combine)
  data <- SeuratObject::AddMetaData(data, metadata = oar_combine)
  
  mdp_combine <- lapply(data_oar, function(x){
    rn = c("mdp")
    names(rn) = paste0("mdp_cl_",levels(Seurat::Idents(x)))
    x[["RNA"]]@meta.data %>% dplyr::select(dplyr::any_of(rn))})
  mdp_combine <- do.call(cbind,mdp_combine)
  
  data[["RNA"]] <- SeuratObject::AddMetaData(data[["RNA"]], mdp_combine)
  
  return(data)
}

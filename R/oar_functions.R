##===================================================================
#Base Oar Test
##===================================================================
#' Generate scores and p-values to determine heterogeneity of data by looking at whether missingness is observed-at-random (OAR)
#'
#' @param data A gene-cell expression matrix with NA values in place of 0s.
#'
#' @return Data frame with OAR-score, p-value, adjusted p-value, and percent missing data for each cell. 
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' output <- oarbase(data)
#' }
oarbase <- function(data) {

  rownames(data) = NULL # remove gene names
  cl = lapply(seq_len(dim(data)[2L]), function(i) data[,i]) #convert mtx to list of cell vectors
  gl = lapply(seq_len(dim(data)[1L]), function(i) data[i,]) #convert mtx to list of gene vectors
  mdp <- find_unique_patterns(gl) # Identify unique missing data patterns
  pvalue.list.KW <- unlist(lapply(cl, FUN = missing_pattern_pval_kw, mdp = mdp)) # Calculate p value
  pvalue.KW.BH = p.adjust(pvalue.list.KW, method = "BH") # Benjamini & Hochberg correction ("BH" or its alias "fdr")
  
  #transform and scale adjusted pvalues
  score = scale(-log10(pvalue.KW.BH), center = T, scale = T)*-1
  
  #calculate missing values
  sp = unlist(lapply(cl, naniar::pct_miss))
  
  #Remove data to free up memory#
  rm(data)
  
  #output
  out <- data.frame(score,pvalue.list.KW,pvalue.KW.BH,sp)
  colnames(out) <- c("OARscore","KW.pvalue","KW.BH.pvalue","pct.missing")
  return(out)
}

##===================================================================
#Iterative Oar Fold Test
##===================================================================
#' Generate OAR fold stats for Seurat object and add them to meta data, run iteratively with different partitions of data and average across them. 
#'
#' @param data A Seurat (v5) object or a matrix with cell barcodes as column names and genes as row names.
#' @param seurat_v5 A Boolean to indicate if supplied data is a Seurat object, default is TRUE
#' @param count.filter A numeric value indicating the minimum fraction of cells expressing any given gene that will be included in the analysis, default is 0.01. Values between 0.005 and 0.02 are recommended.
#' @param blacklisted.genes A character vector with gene names to be excluded from the analysis. Default is empty.
#' @param gene.ratio A numeric value indicating the ratio of genes to cells to partition data matrix. Default is 20. Values between 15 and 25 are recommended.
#' @param iterations A numeric value indicating the number of times to iterate the OARscore calculation. Default and recommendation is 10. Can use 1 to not iterate the calculation. 
#' @param parallel.loop A Boolean to indicate if process should be run in parallel cores. Currently only `doParallel` and `foreach` are supported, default is FALSE
#' @param cores A numeric value indicating the number of cores to use un parallel processing. Use `detectCores()` to identify possibilities. Default is 12. Ignored if `parallel.loop` set to FALSE.
#' @param suffix A string to append to the output variables. Default is empty
#'
#' @return A Seurat object with OAR stats added into meta data
#' @export
#'
#' @examples
#' \dontrun{
#' pmbcs <- oar_fold(pmbcs)
#' }
oar_fold <- function (data, seurat_v5 = TRUE, count.filter = 1,
                      blacklisted.genes = NULL, gene.ratio = 20,
                      iterations = 10, parallel.loop = T,
                      cores = 12, suffix = "") {
  
  print("Analysis started on:")
  print(Sys.time())
  
  #Set filtering threshold
  tr = count.filter/100
  
  if(!is.numeric(count.filter)){stop("count.filter must be numeric\n")}
  if(!is.numeric(iterations)){stop("iterations must be numeric\n")}
  if(parallel.loop && !is.numeric(cores)){stop("cores must be numeric\n")}
  if(!is.numeric(gene.ratio)){stop("gene.ratio must be numeric\n")}
  if(!is.character(suffix)){stop("suffix must be a string\n")}
  if (count.filter > 2) {
    warning("Overfiltering expression matrix (count.filter > 2) may lower signal detection\n")
  }else if(count.filter == 0) {
    #stop("Minimum fraction of cells expressing any given gene must be greater than 0\n")
    warning("Minimum fraction of cells expressing any given gene should be greater than 0\n")
  }
  if (gene.ratio < 15) {
    warning("A small gene to cell ratio may result in unreliable outputs\n")
  }
  if (gene.ratio > 30) {
    warning("A large gene to cell ratio may result in unreliable outputs\n")
  }
  
  #store Seurat object for later
  if(seurat_v5){
    sc.data = data
  }
  
  #read in Seurat object & remove blacklisted genes
  data <- oarpreprocessdata(data, tr, seurat_v5, blacklisted.genes)
  
  #save Cell names
  cell.names = colnames(data)
  
  #Identify number of folds based on ratio
  max.cells = ceiling(data@Dim[[1]]/gene.ratio)
  
  folds <- oaridentifyfolds(data, gene.ratio, max.cells)
  print(paste0("Number of folds: ",folds))
  
  # Notify of time to completion
  print(paste0("Anticipated time to completion: ",
               round(folds*iterations*10/60/cores/2,0)+1," minutes, if using parallel processing"))
  
  #Generate list of cell vectors for parallel or sequential processing
  fold.data = list()
  for(j in 1:iterations){
    a = cell.names
    for(i in 1:folds){
      if(length(a) == 0) break
      s <- paste0("iteration_",j,"_","fold_",i)
      if(length(a) <= max.cells){max.cells = length(a)}
      cells = sample(a,max.cells)
      a = a[!a %chin% cells]
      fold.data[[s]] = data[,colnames(data) %chin% cells]
    }
  }
  
  print("Identifying missing patterns and scoring across folds...")
  #parallel or sequential processing
  if(parallel.loop){
    
    #register cluster
    my.cluster <- parallel::makeCluster(
      cores,
      type = "PSOCK"
    )
    
    doParallel::registerDoParallel(cl = my.cluster) #register it to be used by %dopar%
    
    #Run the loop in parallel
    output <- foreach::foreach(
      f.data = fold.data,
      folds = names(fold.data),
      .verbose = T,
      .packages = c("OAR", "dplyr")) %dopar% {
        print(f.data) #this line fixes s4 subsetting issue and I do not know why
        # Replace 0 with NA
        f.data[f.data == 0] <- NA
        
        # Convert to a dense matrix
        f.data <- as.matrix(f.data)
        
        # store column names
        cells = colnames(f.data)
        colnames(f.data) = NULL
        
        # data must be in .data.frame()
        if (all(complete.cases(f.data))) {
          stop("No missing data exists\n")
        }
        
        # Run test
        w <- oarbase(data = f.data) %>%
          dplyr::mutate(Fold = folds,
                 barcodes = cells)
        return(w)
      }
    parallel::stopCluster(cl = my.cluster)
  }else{
    output <- list()
    n = 1
    for(i in fold.data){
      # Replace 0 with NA
      i[i == 0] <- NA
      
      # Convert to a dense matrix
      i <- i %>% as.matrix()
      
      # store column names
      cells = colnames(i)
      colnames(i) = NULL
      
      # data must be in .data.frame()
      if (all(complete.cases(i))) {
        stop("No missing data exists\n")
      }
      
      # Run test
      output[[n]] <- oarbase(data = i) %>%
        dplyr::mutate(Fold = paste0("fold_",n),
               barcodes = cells)
      n=n+1
    }
  }
  
  print("Score calculation completed")
  
  #Remove data to free up memory#
  rm(fold.data)
  
  print("Consolidating results")
  
  #Consolidate results
  output = do.call(rbind,output) %>%
    dplyr::group_by(barcodes) %>%
    dplyr::summarise(
      OARscore.sd = sd(OARscore),
      OARscore = median(OARscore),
      pct.missing = median(pct.missing)
    )
  
  print("Analysis completed at:")
  print(Sys.time())
  
  # Clean up output
  colnames(output)[-1] <- paste0(colnames(output)[-1],suffix)
  
  #output
  if(seurat_v5){
    idx <- match(rownames(sc.data@meta.data),output$barcodes)
    
    sc.data@meta.data <- base::cbind(
      sc.data@meta.data,
      output[idx,!colnames(output) %in% "barcodes"])
    
    return(sc.data)
  }else{
    return(output)
  }
}

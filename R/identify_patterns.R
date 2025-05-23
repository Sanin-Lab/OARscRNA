##===================================================================
#Calculate hamming distances between genes
##===================================================================
#' Calculate hamming distances between genes
#'
#' @param data a minimal dataset ready for processing.
#' @param cores a numeric value indicating the number of cores to use un parallel processing. Use `parallel::detectCores()` or `parallelly::availableCores()` to identify possibilities.
#'
#' @return A list of matrices of hamming distances across gene bins
#' @export
#'
#' @examples 
#' \dontrun{
#' ##Automatic tolerance setting
#' dm <- oar_hamming_distance(data, cores = 2)
#'
#' }
oar_hamming_distance <- function (data, cores = 1) {
  print("Calculating Hamming distances between gene vectors using specified cores.")
  print("This operation may take several minutes")
  
  bins <- base::rowSums(!is.na(data))/ncol(data) # Define bin ranges
  bins <- base::cut(
    bins, 
    breaks = unique(as.vector(stats::quantile(bins, probs = seq(from = 0, to = 1, by = 0.1)))), # split genes evenly across 10 bins
    labels = labels[1:length(unique(as.vector(stats::quantile(bins, probs = seq(from = 0, to = 1, by = 0.1)))))],
    right = T, include.lowest = T) # Make sure all genes are assigned a bin
  
  d.list <- split(as.data.frame(data),bins) # Split data across bins
  
  dm <- lapply(d.list, function(mm){
    x <- FastHamming::hamming_distance(
      X = +(!is.na(as.matrix(mm))),
      nthreads = cores) # Calculate Hamming distance between gene vectors across data splits
    x <- x/ncol(mm)
    rownames(x) <- 1:nrow(x)
    colnames(x) <- 1:ncol(x)
    return(x)
  })
  
  return(dm)
}

##===================================================================
#Identify missing data patterns allowing for mismatch
##===================================================================
#' Identify missing data patterns allowing for mismatch
#'
#' @param data a minimal dataset ready for processing.
#' @param dm a matrix of hamming distances between gene vectors.
#'
#' @return Vector of missing data patterns
#' @export
#'
#' @examples 
#' \dontrun{
#' mdp <- oar_missing_data_patterns(data, dm)
#' }
oar_missing_data_patterns <- function (data, dm) {
  g <- c()
  for(i in dm){g = c(g,nrow(i))}
  if(!sum(g) == nrow(data)){stop("Hamming distance matrix and count matrix have different number of genes.\nWhere similar filters applied?\n")}
  
  bins <- base::rowSums(!is.na(data))/ncol(data) # Define bin ranges
  bins <- base::cut(
    bins, 
    breaks = unique(as.vector(stats::quantile(bins, probs = seq(from = 0, to = 1, by = 0.1)))), # split genes evenly across 10 bins
    labels = labels[1:length(unique(as.vector(stats::quantile(bins, probs = seq(from = 0, to = 1, by = 0.1)))))],
    right = T, include.lowest = T) # Make sure all genes are assigned a bin
  
  mdp <- lapply(names(dm), function(x){
    if(!nrow(dm[[x]]) == ncol(dm[[x]])){stop("Hamming distance matrix is not square\n")}
    if(!all(diag(dm[[x]]) == 0)){stop("Hamming distance matrix must contain 0s in its diagonal\n")}
    
    vm <- dm[[x]][base::upper.tri(dm[[x]], diag = F)]
    tol <- base::mean(vm) - 4*stats::sd(vm)
    if(tol < 0.01){tol = 0.01}
    
    out <- oar_missing_data_graph(dm[[x]], tol = tol)
    out <- paste(x, out, sep = ".P")
    return(out)
  })
  
  mdp <- unlist(mdp)
  mdp.t <- rep(NA,nrow(data))
  for(i in levels(bins)){mdp.t[bins == i] <- mdp[grepl(pattern = paste0(i,"."), x = mdp, perl = F, fixed = T)]}
  idx <- match(names(table(mdp.t)[(!table(mdp.t)>1)]),mdp.t)
  mdp.t[idx] <- "Unique"
  mdp <- mdp.t
  
  print(paste0("Identified ",length(unique(mdp))-1," non-unique missing data patterns"))
  print(paste0("A total of ",sum(table(mdp)) - sum(table(mdp)["Unique"])," genes captured in non-unique patterns"))
  
  from = unique(mdp)[!unique(mdp) %chin% "Unique"]
  to = paste0("ptn.",1:length(from))
  df <- data.frame("from" = c(from,"Unique"), "to" = c(to,"Unique"))
  idx <- match(mdp,df$from)
  mdp <- df$to[idx]
  
  return(mdp)
}

##===================================================================
# Pull Data Frame of Variable Genes Involved in Patterns
##===================================================================
#' Create list of which genes participate in each pattern.
#'
#' @param data a Seurat object that has had `oar()` or `oar_by_cluster()` run on it previously. 
#'
#' @return data.frame of genes annotated with missing data pattern they participate in (globally or by cluster). 
#' @export
#'
#' @examples 
#' \dontrun{
#' mdp <- get_missing_pattern_genes(data)
#' }
get_missing_pattern_genes <- function (data) {
  mdp <- data@assays$RNA@meta.data
  
  mdp <- mdp %>% 
    dplyr::select(dplyr::contains("mdp")) %>%
    dplyr::mutate(gene_id = rownames(data))
  
  mdp[is.na(mdp)] <- "Filtered"
  
  return(mdp)
  }

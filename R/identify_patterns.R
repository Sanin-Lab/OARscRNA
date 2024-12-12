##===================================================================
#Identify missing data patterns with exact matches
##===================================================================
#' Identify missing data patterns with exact matches and combine all unique missing data patterns
#'
#' @param data a minimal dataset ready for processing.
#'
#' @return missing data pattern vector
#' @export
#'
#' @examples 
#' \dontrun{
#' mdp <- oar_exact_missing_data_patterns(data)
#' }
oar_exact_missing_data_patterns <- function (data) {
  gl = lapply(seq_len(dim(data)[1L]), function(i) data[i,]) #convert mtx to list of gene vectors
  mdp <- unlist(lapply(gl, function(x) digest::digest(1 * is.na(x), algo = "md5"))) # Apply a hashing function to each row to get unique pattern IDs
  mdp <- as.numeric(factor(mdp))
  
  if (length(unique(mdp)) == 1) {stop("Not enough missing data patterns detected\n")}
  freq.mdp <- table(mdp)
  un.p <- names(freq.mdp)[freq.mdp < 2]
  
  if (length(un.p) < dim(freq.mdp)[1L]) {mdp[mdp %in% un.p] <- "unique"}# Combine patterns that occurred once
  if (length(unique(mdp)) == 1) {stop("All missing data patterns were unique.\nConsider allowing for mismatch")}
  
  map <- data.frame(
    "on" = names(table(mdp)),
    "rn" = c(paste0("Pattern ",1:(length(names(table(mdp)))-1)),"Unique patterns"))
  idx <- match(mdp,map$on)
  mdp <- map$rn[idx]
  
  print(paste0(
    "Identified ",
    length(unique(mdp))-1," non-unique missing data patterns, encompassing ",
    sum(table(mdp)) - sum(table(mdp)["Unique patterns"]), " genes"))
  
  return(mdp)
}

##===================================================================
#Identify missing data patterns allowing for mismatch
##===================================================================
#' Identify missing data patterns allowing for mismatch
#'
#' @param data a minimal dataset ready for processing.
#' @param tolerance a boolean or numeric value controlling the tolerance threshold for pattern matching. If set to `TRUE`, then tolerance is automatically adjusted to maximize pattern detection, while minimizing mismatch. Alternatively, user may supply a numeric value indicating the maximum fraction of mismatch between pairs of genes for pattern grouping. Values between 0.01 and 0.05 are recommended.
#' @param cores a numeric value indicating the number of cores to use un parallel processing. Use `parallel::detectCores()` or `parallelly::availableCores()` to identify possibilities.
#'
#' @return Matrix of gene vector hamming distances
#' @export
#'
#' @examples 
#' \dontrun{
#' ##Automatic tolerance setting
#' mdp <- oar_missing_data_patterns(data, tolerance = T)
#'
#' ##User defined
#' mdp <- oar_missing_data_patterns(data, tolerance = 0.05)
#' }
oar_missing_data_patterns <- function (data, tolerance = TRUE, cores = 1) {
  print("Calculating Hamming distances between gene vectors using specified cores.")
  print("This operation may take several minutes")
  dm <- parallelDist::parDist(
    x = +(!is.na(data)), method = "hamming", 
    threads = cores) %>% as.matrix() # Calculate Hamming distance between gene vectors
  
  p = 1
  tol = 0
  if(is.logical(tolerance) && tolerance){
    print("Scanning for optimal tolerance...")
    for(i in seq(0.01,0.05,by = 0.01)){
      mdp <- oar_missing_data_graph(dm, tol = i)
      if(length(unique(mdp)) > p){mdp.f = mdp; tol = i; p = length(unique(mdp))}
    }
    mdp = mdp.f
    print(paste0("Identified ",p-1," non-unique missing data patterns"))
    print(paste0("Tolerance set to ",tol))
    print(paste0("A total of ",sum(table(mdp)) - sum(table(mdp)["Unique patterns"])," genes captured in non-unique patterns"))
  } else if(is.numeric(tolerance)){
    mdp <- oar_missing_data_graph(dm, tol = tolerance)
    print(paste0("Identified ",length(unique(mdp))-1," non-unique missing data patterns"))
    print(paste0("Tolerance set to ",tolerance))
    print(paste0("A total of ",sum(table(mdp)) - sum(table(mdp)["Unique patterns"])," genes captured in non-unique patterns"))
  }else{
    stop("Tolerance must be TRUE or a mumerical value less than 1")
  }
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

##===================================================================
#Identify missing data patterns allowing for mismatch
##===================================================================
#' Identify missing data patterns allowing for mismatch
#'
#' @param data a minimal dataset ready for processing.
#' @param tolerance a boolean or numeric value controlling the tolerance threshold for pattern matching. If set to `TRUE`, then tolerance is automatically adjusted. Alternatively, user may supply a numeric value indicating the maximum fraction of mismatch between pairs of genes for pattern grouping. Values between 0.01 and 0.05 are recommended.
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
#' mdp <- oar_missing_data_patterns(data, tolerance = 0.01)
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
    i = 0.01
    test = TRUE
    mdp <- NULL
    while(test){
      mdp.out <- mdp
      mdp <- oar_missing_data_graph(dm, tol = i)
      lm.out <- oar_base(data, mdp) %>% 
        dplyr::select("x"="pct.missing","y"="KW.BH.pvalue") %>% 
        dplyr::mutate(y = -log10(y)) %>% 
        stats::lm(formula = y~x) %>% 
        summary()
      test <- lm.out$coefficients[2,4] > 0.1
      if(!test && i == 0.01){mdp.out <- mdp}
      if(!test){print(paste0("Tolerance set to ",i))}
      i = i + 0.01
      if(i == 0.05){stop("Setting tolerances avobe 0.05 percent is not recommended")}
    }
    mdp = mdp.out
    print(paste0("Identified ",length(unique(mdp))-1," non-unique missing data patterns"))
    print(paste0("A total of ",sum(table(mdp)) - sum(table(mdp)["unique"])," genes captured in non-unique patterns"))
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

#' Find unique patterns in matrix
#'
#' @param data A list of gene vectors 
#'
#' @return A matrix with counts per pattern
#' @export
#'
#' @examples 
#' \dontrun{
#' mdp <- find_unique_patterns(data)
#' }
find_unique_patterns <- function(data) {
  
  # Apply a hashing function to each row to get unique pattern IDs
  data <- unlist(lapply(data, function(x) digest::digest(1 * is.na(x), algo = "md5")))
  
  # Convert the hash values to factors and then to numeric IDs
  data <- as.numeric(factor(data))
  
  return(matrix(data, ncol = 1))
}

#' Kruskal-Wallis test to generate a per cell p-value based on missing data patterns
#'
#' @param x Item from list of cell gene expression vectors
#' @param mdp Matrix with counts per pattern
#'
#' @return list with a p-value for each cell
#' @export
#'
#' @examples
#' \dontrun{
#' pvalue.list.KW <- unlist(lapply(cl, FUN = missing_pattern_pval_kw, mdp = mdp))
#' }
#' 
missing_pattern_pval_kw = function(x, mdp){
  y.l <- x[!is.na(x)] # subset observed genes of the nth cell to y.l
  mdp.l <- mdp[!is.na(x)] # subset observed genes of the nth cell to y.l
  
  if (length(unique(mdp.l)) == 1) {stop("Not enough missing data patterns for using the tests\n")}
  
  if (sum(!duplicated(mdp.l))>1) {mdp.l[!duplicated(mdp.l)] <- sample(mdp.l[!duplicated(mdp.l)], size = 1)}# Combine patterns that occurred once
  kruskal.test(x = y.l, g = factor(mdp.l))$p.value
}

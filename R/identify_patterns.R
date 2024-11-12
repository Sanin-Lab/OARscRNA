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
  
  freq.mdp <- table(mdp.l)
  un.p <- names(freq.mdp)[freq.mdp < 2]
  if (length(un.p) == length(mdp.l)){warning("NAs generated in OARscore, may be fixed by changing gene.ratio\n")}
  if (length(un.p)<dim(freq.mdp)[1L]) {mdp.l[mdp.l %in% un.p] <- sample(mdp.l[mdp.l %in% un.p], size = 1)}# Combine patterns that occurred once
  if (length(unique(mdp.l)) == 1) {stop("Not enough missing data patterns for using the tests\n")}
  
  kruskal.test(x = y.l, g = factor(mdp.l))$p.value
}

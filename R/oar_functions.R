##========================= Base OAR test ===========================
#' Title
#'
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
oarbase <- function (data) {

  rownames(data) = NULL # remove gene names
  cl = lapply(seq_len(dim(data)[2L]), function(i) data[,i]) #convert mtx to list of cell vectors
  gl = lapply(seq_len(dim(data)[1L]), function(i) data[i,]) #convert mtx to list of gene vectors
  mdp <- find_unique_patterns(gl) # Identify unique missing data patterns
  pvalue.list.KW <- unlist(lapply(cl, FUN = MDP.KW, mdp = mdp)) # Calculate p value
  pvalue.KW.BH = p.adjust(pvalue.list.KW, method = "BH") # Benjamini & Hochberg correction ("BH" or its alias "fdr")
  
  #transform and scale adjusted pvalues
  score = scale(-log10(pvalue.KW.BH), center = T, scale = T)*-1
  
  #calculate missing values
  sp = unlist(lapply(cl, pct_miss))
  
  #Remove data to free up memory#
  rm(data)
  
  #output
  out <- data.frame(score,pvalue.list.KW,pvalue.KW.BH,sp)
  colnames(out) <- c("OARscore","KW.pvalue","KW.BH.pvalue","pct.missing")
  return(out)
}

##===================================================================
#Gene expression prediction based on OAR score
##===================================================================
#' Generate DEGs based on OAR score
#'
#' @param data a Seurat (v5) object or a data.frame with cell barcodes as column names and genes as row names. Seurat object must have a column in metadata with OARscore.
#' @param seurat_v5 a boolean to indicate if supplied data is a Seurat object, default is `TRUE`. If `FALSE`, user must supply a score vector. 
#' @param score a named numeric vector with OARscores to be used in the model. Names should be cell barcodes. Ignored if Seurat object is supplied.
#' @param count.filter a numeric value indicating the minimum fraction of cells expressing any given gene that will be included in the analysis, default is 1. Values between 0.5 and 2 are recommended.
#' @param splines a boolean to indicate if splines should be used to fit model, default is TRUE
#' @param degrees.freedom a numeric value indicating the degrees of freedom to calculate splines. Default is 5. Values between 3 and 5 are recommended.
#' @param blacklisted.genes a character vector with gene names to be excluded from the analysis. Default is empty.
#' @param auto.threshold a boolean to indicate if FDR threshold should be calculated from the data, default is TRUE
#' @param custom.tr a numeric value to use as an FDR threshold. Ignored if `auto.threshold` is set to TRUE.
#' @param score.name name of OAR score in dataset, default is "OARscore". If suffix used on previous functions, would need to include full new name here. 
#'
#' @return A data.frame with p-values for genes that significantly contribute to OAR score. 
#' @export
#'
#' @examples
#' \dontrun{
#' degs <- oar_deg(data)
#' }
oar_deg <- function (data, seurat_v5 = T, score = NULL, count.filter = 1,
                     splines = TRUE, degrees.freedom = 5,
                     blacklisted.genes = NULL, auto.threshold = TRUE, 
                     custom.tr = NULL, score.name = "OARscore")
{
  print("Analysis started on:")
  print(Sys.time())
  
  if(!is.numeric(count.filter)){stop("count.filter must be numeric\n")}
  if(!auto.threshold && !is.numeric(custom.tr)){stop("custom.tr must be numeric\n")}
  if (count.filter > 2) {
    warning("Overfiltering expression matrix (count.filter > 2) may lower signal detection\n")
  }else if(count.filter == 0) {
    stop("Minimum fraction of cells expressing any given gene must be greater than 0\n")
  }
  if(degrees.freedom > 5) {
    stop("More than 5 degrees of freedom are not recommended in spline calculation\n")
  }else if(degrees.freedom <= 0 && splines){
    stop("1 or more degrees of freedom are necessary in spline calculation\n")
  }
  if(seurat_v5 == F && is.null(score)) {stop("Must provide score vector if input is dataframe")}
  
  #Set filtering threshold
  tr = count.filter/100
  
  #Load and filter expression data####
  if(seurat_v5){
    if(!sum(colnames(data@meta.data) %chin% score.name) == 1){
      stop("OARscore column not present, or repeated, in the Seurat object. Have you run oar?\n")
    }
    
    layersList <- lapply(data@assays$RNA@layers,function(x){dim(x)}) #Identify empty layers
    data@assays$RNA@layers[names(layersList[sapply(layersList, is.null)])] <- NULL #remove empty layers
    
    if(length(layersList) >3){
      print("Data may contain multiple raw count layers")
      print("Joining data layers (time consuming!!) and extracting count tables")
    }else{
      print("Extracting count tables")}
    
    dt <- SeuratObject::JoinLayers(object = data, assay = "RNA") %>%
      SeuratObject::LayerData(assay = "RNA", layer = "counts") %>% #Extract count matrix
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column("gene_id")
    
    if(length(blacklisted.genes) > 0 && is.character(blacklisted.genes)){
      dt <- dt %>%
        filter(!gene_id %chin% blacklisted.genes)
    }
    dt <- dt %>%
      filter(!grepl("RPL",gene_id),
             !grepl("Rpl",gene_id),
             !grepl("rpl",gene_id),
             !grepl("RPS",gene_id),
             !grepl("Rps",gene_id),
             !grepl("rps",gene_id))
    
    tr = count.filter/100
    dt <- dt[Matrix::rowSums(dt[,-1]) > dim(dt)[[2]]*tr,]
    
    rownames(dt) <- dt$gene_id
    dt <- dt[,-1]
    
    df <- data@meta.data %>% #Prepare count metadata with residuals as the only information
      dplyr::mutate(ID = rownames(.)) %>%
      dplyr::select("ID", score.name)
    
  }else{
    dt <- data %>%
      rownames_to_column("gene_id")
    
    if(length(blacklisted.genes) > 0 && is.character(blacklisted.genes)){
      dt <- dt %>%
        filter(!gene_id %chin% blacklisted.genes)
    }
    
    dt <- dt %>%
      filter(!grepl("RPL",gene_id),
             !grepl("Rpl",gene_id),
             !grepl("rpl",gene_id),
             !grepl("RPS",gene_id),
             !grepl("Rps",gene_id),
             !grepl("rps",gene_id))
    
    tr = count.filter/100
    dt <- dt[Matrix::rowSums(dt[,-1]) > dim(dt)[[2]]*tr,]
    
    rownames(dt) <- dt$gene_id
    dt <- dt[,-1]
    
    df <- data.frame(
      "ID" = names(score),
      score.name = score)
  }
  
  print("FDR threshold set to:")
  if(!auto.threshold){threshold = custom.tr}else{threshold = 10^(-(ceiling(log10(dim(dt)[[2]]))*3))}
  print(threshold)
  
  #Set up models####
  if(splines){
    warning("Using splines increases calculation time by 3-5x\n")
    
    X <- splines::ns(df[, score.name], df = degrees.freedom) # Any df between 3 - 5 usually works well.
    mm <- model.matrix(~X) # Model matrix with splines
  }else{
    mm <- model.matrix(~df[, score.name]) #Model no splines
  }
  
  #Run differential analysis####
  sf <- scran::calculateSumFactors(dt %>% as.matrix())
  
  fit <- glmGamPoi::glm_gp(
    data = dt %>% as.matrix(),
    design = mm,
    col_data = df,
    size_factors = sf,
    on_disk = F)
  
  out <- glmGamPoi::test_de(
    fit, reduced_design = ~1, sort_by = "adj_pval") %>%
    filter(adj_pval < threshold)
  
  print("Analysis completed at:")
  print(Sys.time())
  
  #output
  return(out)
}


##===============================================
# Score v Percent Missing Scatter Plot
##===============================================
#' Create scatter plot of OAR score vs percent missing
#'
#' @param data a seurat v5 object that has OAR score in meta data, or a data.frame with the OAR score results.
#' @param group.by a meta data category to color data by. Default is seurat_clusters. 
#' @param seurat_v5 a boolean to indicate if supplied data is a Seurat object, default is TRUE
#' @param suffix a string that was previously appended to the output variables. Default is empty.
#' @param pt.size a numerical value for the size of points to be passed to the size argument in `geom_point`. Default is 0.5.
#' 
#' @return Scatter plot of OAR score vs percent missing, colored by grouping of choice. 
#' @export
#'
#' @examples
#' \dontrun{
#' ##Starting from a Seurat object
#' output <- scatter_score_missing(pmbcs_oar)
#' 
#' ##Starting from a oar results data.frame
#' output <- scatter_score_missing(oar, seurat_v5 = F)
#' }
scatter_score_missing <- function(
    data, group.by = 'seurat_clusters', seurat_v5 = TRUE, suffix = "",
    pt.size = 0.5) {
  if(seurat_v5){
    input_data <- data@meta.data
    
    plot <- ggplot(
      data = input_data,
      aes(x = .data[[paste0("pct.missing", suffix)]],
          y = .data[[paste0("OARscore", suffix)]],
          color = .data[[group.by]])) +
      geom_point(size = pt.size) +
      theme_classic() +
      labs(x = "% Missing Values",
           y = paste0("OAR score", suffix)) +
      theme(aspect.ratio = 1)
    
  }else{
    input_data <- data
    
    plot <- ggplot(
      data = input_data,
      aes(x = .data[[paste0("pct.missing", suffix)]],
          y = .data[[paste0("OARscore", suffix)]])) +
      geom_point(size = pt.size) +
      theme_classic() +
      labs(x = "% Missing Values",
           y = paste0("OAR score", suffix)) +
      theme(aspect.ratio = 1)
  }
  
 return(plot)
}

##===================================================================
#Plot gene patterns in data
##===================================================================
#' Plot identified missing data patterns
#'
#' @param data a gene-cell expression matrix with NA values in place of 0s and 1s everywhere else or a Seurat Object to which `oar` has been applied to.
#' @param mdp a vector indicating the pattern to which each gene belongs. Default is NULL.
#' @param seurat_v5 a boolean to indicate where or not the input data is a Seurat object. Default is TRUE. 
#'
#' @return Plot of missing data patterns. 
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' ##Starting from a Seurat Object analysed jointly
#' output <- oar_missing_data_plot(data, seurat_v5 = TRUE)
#' 
#' ##Starting from filtered and binarized expression matrix
#' output <- oar_missing_data_plot(data, mdp, seurat_v5 = FALSE)
#' }
oar_missing_data_plot <- function(data, mdp = NULL, seurat_v5 = TRUE) {
  if(seurat_v5){
    mdp <- get_missing_pattern_genes(data) %>% 
      dplyr::select(dplyr::all_of("mdp")) %>% 
      dplyr::filter(!mdp == "Filtered") %>% 
      .$mdp
    
    output <- oar_preprocess_data(data)
    data <- output[[1]]
    
  }else{
    data <- data
    mdp <- mdp
  }
  
  m <- +(!is.na(data)) %>% 
    reshape2::melt() %>% 
    dplyr::mutate(
      x = as.factor(Var2),
      group = as.factor(
        rep(mdp,
            times = length(unique(Var2))))) %>% 
    dplyr::filter(!group == "Unique") %>% 
    dplyr::group_by(group, x) %>% 
    dplyr::summarise(p = sum(value)/n(),
                     count = n()) %>% 
    dplyr::mutate(group = paste0(group,": ",count," genes")) %>% 
    reshape2::acast(formula = group ~ x, value.var = "p")
  
  print("Generating plots...")
  ht <- ComplexHeatmap::Heatmap(
    matrix = m, 
    name = NULL, heatmap_legend_param = list(
      title = NULL, at = c(0, 1), 
      labels = c("Missing", "Observed"),
      legend_height = unit(1, "cm"),
      grid_width = unit(2, "mm"),
      labels_gp = grid::gpar(fontsize = 8),
      border = "black"),
    col = circlize::colorRamp2(
      c(0,0.5,1), colors = c("white","grey","black"),
      transparency = 0),
    border = T, show_row_dend = F, show_column_dend = F,
    show_column_names = F, show_row_names = T,
    cluster_rows = T, cluster_columns = T,
    row_labels = rownames(m),
    row_names_side = "left", row_names_gp = grid::gpar(fontsize = 8),
    column_title = "Cells", column_title_side = "top",
    column_title_gp = grid::gpar(fontsize = 8),
    width = unit(ncol(m)*0.01, "mm"), height = unit(3*nrow(m), "mm"))
  
  # Output #
  return(ht)
}
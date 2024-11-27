##===============================================
# Score v Percent Missing Scatter Plot
##===============================================
#' Create scatter plot of OAR score vs percent missing
#'
#' @param data A seurat v5 object that has oar score in meta data.
#' @param group.by meta data category to color data by. Default is seurat_clusters. 
#' @param seurat_v5 A Boolean to indicate if supplied data is a Seurat object, default is TRUE
#' @param suffix A string that was previously appended to the output variables. Default is empty.
#' 
#' @return Scatter plot of OAR score vs percent missing, colored by grouping of choice. 
#' @export
#'
#' @examples
#' \dontrun{
#' pmbcs_oar <- scatter_score_missing(pmbcs_oar)
#' }
scatter_score_missing <- function(data, group.by = 'seurat_clusters', seurat_v5 = TRUE, suffix = "") {
  if(seurat_v5){
    input_data <- data@meta.data
    
    plot <- ggplot(
      data = input_data,
      aes(x = .data[[paste0("pct.missing", suffix)]],
          y = .data[[paste0("OARscore", suffix)]],
          color = .data[[group.by]])) +
      geom_point(size = 1) +
      theme_classic() +
      labs(x = "% Missing Values",
           y = paste0("OAR score", suffix)) +
      geom_hline(yintercept = 0,
                 linetype = "dashed",
                 color = "black", linewidth = 0.5) +
      theme(aspect.ratio = 1)
    
  }else{
    input_data <- data
    
    plot <- ggplot(
      data = input_data,
      aes(x = .data[[paste0("pct.missing", suffix)]],
          y = .data[[paste0("OARscore", suffix)]])) +
      geom_point(size = 1) +
      theme_classic() +
      labs(x = "% Missing Values",
           y = paste0("OAR score", suffix)) +
      geom_hline(yintercept = 0,
                 linetype = "dashed",
                 color = "black", linewidth = 0.5) +
      theme(aspect.ratio = 1)
  }
  
 return(plot)
}

##===================================================================
#Plot gene patterns in data
##===================================================================
#' Plot identified missing data patterns
#'
#' @param data A gene-cell expression matrix with NA values in place of 0s.
#' @param mdp A vector indicating the pattern to which each gene belongs. Default is Null.
#' @param seurat_v5 A Boolean to indicate where or not the input data is a Seurat object. Default is T. 
#'
#' @return Plot of missing data patterns. 
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' output <- oar_missing_data_plot(data, mdp)
#' }
oar_missing_data_plot <- function(data, mdp = NULL, seurat_v5 = TRUE) {
  if(seurat_v5){
    mdp <- data@assays$RNA@meta.data$mdp
    output <- oar_preprocess_data(data)
    data <- output[[1]]
  }else{
    data <- data
    mdp <- mdp
  }
  
  freq.mdp <- table(mdp)
  p <- names(freq.mdp)
  
  gr = c()
  for(i in p){
    gr = c(gr,
           rep(i,times=freq.mdp[i]))}
  
  dout <- +(!is.na(data)) %>% 
    reshape2::melt() %>% 
    dplyr::mutate(
      x = as.factor(Var2),
      group = as.factor(
        rep(mdp[mdp %in% p],
            times = length(unique(Var2))))) %>% 
    dplyr::group_by(group, x) %>% 
    dplyr::summarise(p = sum(value)/n(),
                     count = n())
  
  print("Generating plots...")
  pp <- ggplot2::ggplot(
    data = dout,
    aes(x = x, y = group, fill = p)) +
    ggplot2::geom_tile(height = 0.5) +
    ggplot2::scale_fill_gradient(
      low = "white", high = "black",
      limits = c(0,1),
      breaks = c(0,1),
      labels = c("Missing","Observed"),
      guide = guide_colorbar(
        title = NULL,
        frame.colour = "black",
        ticks.colour = "black")) +
    ggplot2::geom_label(
      data = dout %>% 
        dplyr::distinct(group,count) %>% 
        dplyr::mutate(count = paste0("Genes in pattern ",group,": ",count)), 
      aes(x = 1, y = group, label = count),
      inherit.aes = F, size.unit = "mm", size =3,
      nudge_y = 0.5, nudge_x = 150) +
    ggplot2::theme_minimal() +
    ggplot2::coord_cartesian(clip = "off") +
    theme(
      panel.spacing = unit(0,"mm"),
      legend.position = "right",
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      legend.key.size = unit(2,"mm")
    )
  
  # Output #
  return(pp)
}
#' Create scatter plot of OAR score vs OAR score Standard Deviation
#'
#' @param data A seurat v5 object that has oar score in meta data.
#' @param seurat_v5 A Boolean to indicate if supplied data is a Seurat object, default is TRUE
#'
#' @return Scatter plot of OAR score standard deviation vs OAR score, colored by percent missingness of each cell. 
#' @export 
#'
#' @examples
#' \dontrun{
#' pmbcs_oar <- scatter_score_sd(pmbcs_oar)
#' }
scatter_score_sd <- function(data, seurat_v5 = TRUE){
  
  if(seurat_v5){
    input_data <- data@meta.data
  }
  
  plot <- ggplot(
    data = input_data,
    aes(x = OARscore.sd,
        y = OARscore,
        color = pct.missing)) +
    ggrastr::rasterise(geom_point(size = 0.2), dpi = 300) +
    theme_classic() +
    scale_color_gradientn(
      colours = rev(c("#000004FF","#3B0F70FF","#8C2981FF",
                      "#DE4968FF", "#FE9F6DFF","#FCFDBFFF")),
      na.value = "transparent",
      name = "% Missing values",
      guide = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black")) +
    labs(x = "OAR score SD",
         y = "OAR score Median") +
    geom_hline(yintercept = 2,
               linetype = "dashed",
               color = "black", linewidth = 0.5) +
    theme(aspect.ratio = 1,
          legend.key.width = unit(3,"mm"))
  return(plot)
}

#' Create scatter plot of OAR score vs percent missing
#'
#' @param data A seurat v5 object that has oar score in meta data.
#' @param group.by meta data category to color data by. Default is seurat_clusters. 
#' @param seurat_v5 A Boolean to indicate if supplied data is a Seurat object, default is TRUE
#'
#' @return Scatter plot of OAR score vs percent missing, colored by grouping of choice. 
#' @export
#'
#' @examples
#' \dontrun{
#' pmbcs_oar <- scatter_score_missing(pmbcs_oar)
#' }
scatter_score_missing <- function(data, group.by = 'seurat_clusters', seurat_v5 = TRUE) {
  if(seurat_v5){
    input_data <- data@meta.data
  }
  
   plot <- ggplot(
    data = input_data,
    aes(x = pct.missing,
        y = OARscore,
        color = .data[[group.by]])) +
    geom_point(size = 1) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow=2, title = NULL))+
    theme_classic() +
    labs(x = "% Missing Values",
         y = "OAR score") +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color = "black", linewidth = 0.5) +
    theme(aspect.ratio = 1,
          legend.key.width = unit(3,"mm"))
   
 return(plot)
}
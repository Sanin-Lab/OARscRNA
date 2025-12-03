# Create scatter plot of OAR score vs sparsity

Create scatter plot of OAR score vs sparsity

## Usage

``` r
scatter_score(
  data,
  group.by = "seurat_clusters",
  seurat_v5 = TRUE,
  suffix = "",
  pt.size = 0.5
)
```

## Arguments

- data:

  a seurat v5 object that has OAR score in meta data, or a data.frame
  with the OAR score results.

- group.by:

  a meta data category to color data by. Default is seurat_clusters.

- seurat_v5:

  a boolean to indicate if supplied data is a Seurat object, default is
  TRUE

- suffix:

  a string that was previously appended to the output variables. Default
  is empty.

- pt.size:

  a numerical value for the size of points to be passed to the size
  argument in `geom_point`. Default is 0.5.

## Value

Scatter plot of OAR score vs sparsity, colored by grouping of choice.

## Examples

``` r
if (FALSE) { # \dontrun{
##Starting from a Seurat object
output <- scatter_score(pmbcs_oar)

##Starting from a oar results data.frame
output <- scatter_score(oar, seurat_v5 = F)
} # }
```

# Plot identified gene co-expression patterns

Plot identified gene co-expression patterns

## Usage

``` r
oar_gcp_plot(data, gcp = NULL, seurat_v5 = TRUE)
```

## Arguments

- data:

  a gene-cell expression matrix with NA values in place of 0s and 1s
  everywhere else or a Seurat Object to which `oar` has been applied to.

- gcp:

  a vector indicating the pattern to which each gene belongs. Default is
  NULL.

- seurat_v5:

  a boolean to indicate where or not the input data is a Seurat object.
  Default is TRUE.

## Value

Plot of gene co-expression patterns.

## Examples

``` r
if (FALSE) { # \dontrun{
##Starting from a Seurat Object analysed jointly
output <- oar_gcp_plot(data, seurat_v5 = TRUE)

##Starting from filtered and binarized expression matrix
output <- oar_gcp_plot(data, gcp, seurat_v5 = FALSE)
} # }
```

# Prepare data for oar fold functions

Prepare data for oar fold functions

## Usage

``` r
oar_preprocess_data(data, tr = 1, seurat_v5 = TRUE, blacklisted.genes = NULL)
```

## Arguments

- data:

  a seurat object or gene expression matrix

- tr:

  a filtering threshold. Default is 1.

- seurat_v5:

  a boolean to indicate if supplied data is a Seurat object, default is
  TRUE

- blacklisted.genes:

  a character vector with gene names to be excluded from the analysis.
  Default is empty.

## Value

Data matrix with blacklisted genes removed

## Examples

``` r
if (FALSE) { # \dontrun{
data <- oar_preprocess_data(data)
} # }
```

# Single line pipeline to run complete analysis

Single line pipeline to run complete analysis

## Usage

``` r
oar(
  data,
  seurat_v5 = TRUE,
  count.filter = 1,
  blacklisted.genes = NULL,
  suffix = "",
  cores = 1,
  store.hamming = TRUE
)
```

## Arguments

- data:

  a Seurat (v5) object or a matrix with cell barcodes as column names
  and genes as row names.

- seurat_v5:

  a boolean to indicate if supplied data is a Seurat object, default is
  TRUE

- count.filter:

  a numeric value indicating the minimum fraction of cells expressing
  any given gene that will be included in the analysis, default is 1.
  Values between 0.5 and 2 are recommended.

- blacklisted.genes:

  a character vector with gene names to be excluded from the analysis.
  Default is empty.

- suffix:

  a string to append to the output variables. Default is empty

- cores:

  a numeric value indicating the number of cores to use un parallel
  processing. Use
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  or
  [`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html)
  to identify possibilities. Default is 1.

- store.hamming:

  a boolean to control if hamming distances should be stored in the
  Seurat object. Default set to `TRUE`

## Value

A Seurat object with OAR stats added into meta data, or a matrix with
OAR stats.

## Examples

``` r
if (FALSE) { # \dontrun{
pbmcs <- oar(pbmcs)
} # }
```

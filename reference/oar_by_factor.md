# Generate OAR score within each cluster and add them to full objects metadata

Generate OAR score within each cluster and add them to full objects
metadata

## Usage

``` r
oar_by_factor(
  data,
  count.filter = 1,
  blacklisted.genes = NULL,
  suffix = "",
  factor = "ident",
  cores = 1
)
```

## Arguments

- data:

  a Seurat (v5) object.

- count.filter:

  a numeric value indicating the minimum fraction of cells expressing
  any given gene that will be included in the analysis, default is 1.
  Values between 0.5 and 2 are recommended.

- blacklisted.genes:

  a character vector with gene names to be excluded from the analysis.
  Default is empty.

- suffix:

  a string to append to the output variables. Default is empty

- factor:

  a character string mapping to a column in `meta.data`. Defaults to
  `ident`.

- cores:

  a numeric value indicating the number of cores to use un parallel
  processing. Use
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  or
  [`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html)
  to identify possibilities.

## Value

a Seurat object with OAR stats added into meta data

## Examples

``` r
if (FALSE) { # \dontrun{
pmbcs <- oar_by_cluster(pmbcs)
} # }
```

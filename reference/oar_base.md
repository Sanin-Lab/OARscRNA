# Generate scores and p-values to determine transcriptional shifts

Generate scores and p-values to determine transcriptional shifts

## Usage

``` r
oar_base(data, gcp)
```

## Arguments

- data:

  A gene-cell expression matrix with NA values in place of 0s.

- gcp:

  A vector indicating the pattern to which each gene belongs.

## Value

Data frame with OAR-score, p-value, adjusted p-value, and sparsity as a
percentage for each cell.

## Examples

``` r
if (FALSE) { # \dontrun{
output <- oar_base(data, gcp)
} # }
```

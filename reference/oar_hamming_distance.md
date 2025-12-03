# Calculate hamming distances between genes

Calculate hamming distances between genes

## Usage

``` r
oar_hamming_distance(data, cores = 1)
```

## Arguments

- data:

  a minimal dataset ready for processing.

- cores:

  a numeric value indicating the number of cores to use un parallel
  processing. Use
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  or
  [`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html)
  to identify possibilities.

## Value

A list of matrices of hamming distances across gene bins

## Examples

``` r
if (FALSE) { # \dontrun{
##Automatic tolerance setting
dm <- oar_hamming_distance(data, cores = 2)

} # }
```

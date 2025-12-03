# Group gene co-expression patterns based on tolerance with a graph

Group gene co-expression patterns based on tolerance with a graph

## Usage

``` r
oar_gene_graph(dm, tol)
```

## Arguments

- dm:

  a matrix of gene vector hamming distances.

- tol:

  a numeric value indicating the maximum fraction of mismatch in genes
  to group as a pattern.

## Value

gene co-expression pattern vector

## Examples

``` r
if (FALSE) { # \dontrun{
gcp <- oar_gene_graph(dm, tol)
} # }
```

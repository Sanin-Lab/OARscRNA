# Kruskal-Wallis test to generate a per cell p-value based on gene co-expression patterns

Kruskal-Wallis test to generate a per cell p-value based on gene
co-expression patterns

## Usage

``` r
pattern_pval_kw(x, gcp)
```

## Arguments

- x:

  Item from list of cell gene expression vectors

- gcp:

  Matrix with gene participation per pattern

## Value

list with a p-value for each cell

## Examples

``` r
if (FALSE) { # \dontrun{
pvalue <- pattern_pval_kw(x, gcp)
} # }
```

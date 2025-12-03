# Create list of which genes participate in each pattern.

Create list of which genes participate in each pattern.

## Usage

``` r
get_pattern_genes(data)
```

## Arguments

- data:

  a Seurat object that has had
  [`oar()`](https://sanin-lab.github.io/OARscRNA/reference/oar.md) or
  `oar_by_cluster()` run on it previously.

## Value

data.frame of genes annotated with gene co-expression pattern they
participate in (globally or by cluster).

## Examples

``` r
if (FALSE) { # \dontrun{
gcp <- get_pattern_genes(data)
} # }
```

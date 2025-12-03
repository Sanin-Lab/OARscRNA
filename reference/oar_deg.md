# Generate DEGs based on OAR score

Generate DEGs based on OAR score

## Usage

``` r
oar_deg(
  data,
  seurat_v5 = T,
  score = NULL,
  count.filter = 1,
  splines = TRUE,
  degrees.freedom = 5,
  blacklisted.genes = NULL,
  auto.threshold = TRUE,
  custom.tr = NULL,
  score.name = "OARscore"
)
```

## Arguments

- data:

  a Seurat (v5) object or a data.frame with cell barcodes as column
  names and genes as row names. Seurat object must have a column in
  metadata with OARscore.

- seurat_v5:

  a boolean to indicate if supplied data is a Seurat object, default is
  `TRUE`. If `FALSE`, user must supply a score vector.

- score:

  a named numeric vector with OARscores to be used in the model. Names
  should be cell barcodes. Ignored if Seurat object is supplied.

- count.filter:

  a numeric value indicating the minimum fraction of cells expressing
  any given gene that will be included in the analysis, default is 1.
  Values between 0.5 and 2 are recommended.

- splines:

  a boolean to indicate if splines should be used to fit model, default
  is TRUE

- degrees.freedom:

  a numeric value indicating the degrees of freedom to calculate
  splines. Default is 5. Values between 3 and 5 are recommended.

- blacklisted.genes:

  a character vector with gene names to be excluded from the analysis.
  Default is empty.

- auto.threshold:

  a boolean to indicate if FDR threshold should be calculated from the
  data, default is TRUE

- custom.tr:

  a numeric value to use as an FDR threshold. Ignored if
  `auto.threshold` is set to TRUE.

- score.name:

  name of OAR score in dataset, default is "OARscore". If suffix used on
  previous functions, would need to include full new name here.

## Value

A data.frame with p-values for genes that significantly contribute to
OAR score.

## Examples

``` r
if (FALSE) { # \dontrun{
degs <- oar_deg(data)
} # }
```

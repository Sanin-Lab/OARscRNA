# Package index

## Wrapper Functions

These are all the functions you need to run the full analysis on your
data.

- [`oar()`](https://sanin-lab.github.io/OARscRNA/reference/oar.md) :
  Single line pipeline to run complete analysis
- [`oar_by_factor()`](https://sanin-lab.github.io/OARscRNA/reference/oar_by_factor.md)
  : Generate OAR score within each cluster and add them to full objects
  metadata
- [`oar_deg()`](https://sanin-lab.github.io/OARscRNA/reference/oar_deg.md)
  : Generate DEGs based on OAR score
- [`get_pattern_genes()`](https://sanin-lab.github.io/OARscRNA/reference/get_pattern_genes.md)
  : Create list of which genes participate in each pattern.

### Visualize results

A few convenient functions to visualize the results of your analysis.

- [`scatter_score()`](https://sanin-lab.github.io/OARscRNA/reference/scatter_score.md)
  : Create scatter plot of OAR score vs sparsity
- [`oar_gcp_plot()`](https://sanin-lab.github.io/OARscRNA/reference/oar_gcp_plot.md)
  : Plot identified gene co-expression patterns

## Step-by-step functions

Run these in this order to achieve the same results as with the wrapper
functions.

- [`oar_preprocess_data()`](https://sanin-lab.github.io/OARscRNA/reference/oar_preprocess_data.md)
  : Prepare data for oar fold functions
- [`oar_hamming_distance()`](https://sanin-lab.github.io/OARscRNA/reference/oar_hamming_distance.md)
  : Calculate hamming distances between genes
- [`oar_patterns()`](https://sanin-lab.github.io/OARscRNA/reference/oar_patterns.md)
  : Identify gene co-expression patterns allowing for mismatch
- [`oar_base()`](https://sanin-lab.github.io/OARscRNA/reference/oar_base.md)
  : Generate scores and p-values to determine transcriptional shifts

## Utility functions

Called within other functions in the package.

- [`pattern_pval_kw()`](https://sanin-lab.github.io/OARscRNA/reference/pattern_pval_kw.md)
  : Kruskal-Wallis test to generate a per cell p-value based on gene
  co-expression patterns
- [`oar_gene_graph()`](https://sanin-lab.github.io/OARscRNA/reference/oar_gene_graph.md)
  : Group gene co-expression patterns based on tolerance with a graph

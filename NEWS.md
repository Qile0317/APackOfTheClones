# APackOfTheClones 0.1.0 
Initial semi-stable release. The main functions of the package are working with their default parameters. The package currently exports three functions:
* `integrate_tcr`
* `count_clone_sizes`
* `clonal_expansion_plot`

All of which have documentation, including on the documentation site.

Currently, Automated cluster repulsion uses a highly flawed mathematical formula loosely based on force directed graph drawing techniques, the resulting plots are not presentable and the function has been disabled temporarily

Upcoming features/plans are:
* Customizable cluster coloration, not just based on the original Seurat/ggplot palette
* User-controlled cluster shifting
* Comprehensive cluster repulsion
* better legends and optional legend border
* automated optimization of the initial parameters of `clonal_expansion_plot`, especially the `clone_scale_factor`
* BCR library integration
* User vignette
* Developer vignette
* Submission to CRAN or Bioconductor
* Rewrite of certain circle packing functions in `rust` to improve performance
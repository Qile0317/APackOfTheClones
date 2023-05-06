# APackOfTheClones 0.1.1
## Additions
* A preliminary user vignette can be viewed on the documentation website, but the code chunks are missing on the site. It is unclear currently if it can be viewed in R either. Future patches will  improve this.

## Changes
* An initial version of automated cluster repulsion is correctly implemented and the related arguments in `clonal_expansion_plot` has been "unlocked".
* Size legends are now correctly implemented, and the related arguments have also been "unlocked" in `clonal_expansion_plot`.
* Function documentation and the README in general has been slightly improved. (automatically updated on the pkgdown website)

# APackOfTheClones 0.1.0 
Initial semi-stable release. The main functions of the package are working with their default parameters. The package currently exports three functions:

* `integrate_tcr`
* `count_clone_sizes`
* `clonal_expansion_plot`

All of which have documentation, including on the documentation site.

Currently, Automated cluster repulsion uses a highly flawed mathematical formula loosely based on force directed graph drawing techniques, the resulting plots are not presentable and the function has been disabled temporarily

##  Future features/plans:
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
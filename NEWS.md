# APackOfTheClones (development version)

# APackOfTheClones 0.1.2
## Additions
* Two new datasets previously in the test folder have now moved to `\data` and are available to users and can be called via `data("mini_seurat_obj")` and `data("mini_clonotype_data")`.

## Changes
* Released package to CRAN .
* Minor documentation changes to adhere to CRAN standards.

# APackOfTheClones 0.1.1
## Additions
* A preliminary user vignette can be viewed on the documentation website, but not locally.

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
[ ] Customizable cluster coloration, not just based on the original Seurat/ggplot palette
[ ] User-controlled cluster shifting
[x] Comprehensive cluster repulsion
[x] better legends and optional legend border
[ ] automated optimization of the initial parameters of `clonal_expansion_plot`, especially the `clone_scale_factor`
[ ] BCR library integration
[x] User vignette
[ ] Developer vignette
[ ] Submission to CRAN or Bioconductor
[ ] Rewrite of certain circle packing functions in `rust` or `c++` to improve performance
[ ] Increase interoperability with `Seurat` and `scRepertoire`

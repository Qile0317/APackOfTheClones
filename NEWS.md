# APackOfTheClones (development version)
## Additions
* There is now a new argument in `clonal_expansion_plot` called `scramble` which allows the user to optionally make the clones within each clsuter randomly distributed.

## Changes
* When working with size legend positions within the plot with the `legend_position` argument in `clonal_expansion_plot`, it can now simply just be a numeric vector of length 2 indicating the x and y coordinate of the "top center" of the legend.

# APackOfTheClones 0.1.4
## Additions
* The clonal expansion plotting function now has a new optional argument `reduction` which now allows the user to choose which reduction the circle clonal clusters should be based on. Defaults to `'umap'` but can now be changed to `'tsne'` or `'pca'` given that they have been ran already on the seurat object. The vignette has been updated accordingly as well
* New CRAN badge showing date of latest release on `README.md`

## Changes
* Some testcases have been altered/fixed to pass CRAN's R-CMD check
* Package startup message has been improved

# APackOfTheClones 0.1.3
## Additions
* `clonal_expansion_plot` now displays time elapsed if `verbose = TRUE`
* There are now package startup messages

## Changes
* Many functions have been rewritten in C++, improving the circle packing runtime of `clonal_expansion_plot` 70 ~ 140 fold. For reference, for a dataset with around 2500 viable T cells, the runtime averages 0.1s

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

##  Future features/plans: (~ means in progress)
[~] Customizable cluster coloration, not just based on the original Seurat/ggplot palette

[~] User-controlled cluster shifting

[ ] automated optimization of the initial parameters of `clonal_expansion_plot`, especially the `clone_scale_factor`

[ ] BCR library integration

[ ] Developer vignette / manuscript

[~] Rewrite of certain circle packing functions in `rust` or `c++` to improve performance

[ ] Increase interoperability with `Seurat` and `scRepertoire`

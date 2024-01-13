# APackOfTheClones 1.0.0

## Additions

- Add `vizAPOTC` to replace `clonal_expansion_plot` with some similar arguments but works completely differently with an additional set of argument requirements, to ultimately produce a similar but better visualization.
- Add `RunAPOTC` to store multiple APackOfTheClones runs on the seurat object.
- Add `containsApotcRun` to check for the existence of some run based on an ID.
- Add `getApotcDataIds` to get all APackOfTheClones run IDs.
- Add `getLastApotcDataId` to get the latest APackOfTheClones run ID.
- Add `deleteApotcData` to delete all data associate with some APackOfTheClones run.
- Add `APOTCPlot` to plot APackOfTheClones run results into a ggplot.
- Add `AdjustAPOTC` to modify certain attributes of prior APackOfTheClones runs.
- Add `countCloneSizes` as a convenience function to count clonotype sizes based on the clonecall and allows subsetting of the seurat metadata
- Add `getReductionCentroids` as a convenience function to get the physical centroids of each cell in some reduction in a seurat object.
- Added two new vignettes are on the `pkgdown` website and indexed, and should be viewed in order of `vignette("APackOfTheClones")` and `vignette("APackOfTheClones-runs")`.

## Changes

- ***ALL v0.1.x functions are deprecated/defunct***, and the clonal expansion workflow has changed.
- The latest version two of `scRepertoire` is now a "dependency", and the package is now structured as an extension of `scRepertoire` when it comes to integrating VDJ data, as `scRepertoire` already possesses a mature and common method to do so.
- Package startup message has been modified to notify users about the deprecation and breaking changes.
- All clonal expansion plots will no longer have a title for the default theme

# APackOfTheClones 0.1.4 (github release)

## Additions

- The clonal expansion plotting function now has a new optional argument `reduction` which now allows the user to choose which reduction the circle clonal clusters should be based on. Defaults to `'umap'` but can now be changed to `'tsne'` or `'pca'` given that they have been ran already on the seurat object. The vignette has been updated accordingly as well
- New CRAN badge showing date of latest release on `README.md`

## Changes

- Some testcases have been altered/fixed to pass CRAN's R-CMD check
- Package startup message has been improved

# APackOfTheClones 0.1.3

## Additions

- `clonal_expansion_plot` now displays time elapsed if `verbose = TRUE`
- There are now package startup messages

## Changes

- Many functions have been rewritten in C++, improving the circle packing runtime of `clonal_expansion_plot` 70 ~ 140 fold. For reference, for a dataset with around 2500 viable T cells, the runtime averages 0.1s

# APackOfTheClones 0.1.2

## Additions

- Two new datasets previously in the test folder have now moved to `\data` and are available to users and can be called via `data("mini_seurat_obj")` and `data("mini_clonotype_data")`.

## Changes

- Released package to CRAN.
- Minor documentation changes to adhere to CRAN standards.

# APackOfTheClones 0.1.1 (github release)

## Additions

- A preliminary user vignette can be viewed on the documentation website, but not locally.

## Changes

- An initial version of automated cluster repulsion is correctly implemented and the related arguments in `clonal_expansion_plot` has been "unlocked".
- Size legends are now correctly implemented, and the related arguments have also been "unlocked" in `clonal_expansion_plot`.
- Function documentation and the README in general has been slightly improved. (automatically updated on the pkgdown website)

# APackOfTheClones 0.1.0 (github release)

Initial semi-stable release. The main functions of the package are working with their default parameters. The package currently exports three functions:

- `integrate_tcr`
- `count_clone_sizes`
- `clonal_expansion_plot`

All of which have documentation, including on the documentation site.

Currently, Automated cluster repulsion uses a highly flawed mathematical formula loosely based on force directed graph drawing techniques, the resulting plots are not presentable and the function has been disabled temporarily

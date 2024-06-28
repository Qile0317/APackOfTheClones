# APackOfTheClones <img src="man/figures/logo.png" width="20%" align="right" />

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/APackOfTheClones)](https://CRAN.R-project.org/package=APackOfTheClones)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/APackOfTheClones?color=brightgreen)](https://www.r-pkg.org/pkg/APackOfTheClones)
[![R-CMD-check](https://github.com/Qile0317/APackOfTheClones/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Qile0317/APackOfTheClones/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/Qile0317/APackOfTheClones/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Qile0317/APackOfTheClones?branch=main)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://qile0317.github.io/APackOfTheClones/)
[![JOSS status](https://joss.theoj.org/papers/eb9d2b2fc3ffb1d440bd7120b66f1f71/status.svg)](https://joss.theoj.org/papers/eb9d2b2fc3ffb1d440bd7120b66f1f71)
[![Developmental Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://qile0317.github.io/APackOfTheClones/dev/)
[![R Universe APackOfTheClones status badge](https://qile0317.r-universe.dev/badges/APackOfTheClones)](https://qile0317.r-universe.dev/APackOfTheClones)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/Qile0317/APackOfTheClones/blob/main/LICENSE.md)
<!-- badges: end -->

```APackOfTheClones``` is an R package that extends the bioconductor ```scRepertoire``` package to produce easily customizable "ball-packing" visualizations of the clonal expansion of T-cells/B-cells in a `Seurat` object, based on its receptor library and single cell RNA sequencing data (for example outputs from 10X genomics' single-cell immune profiling).

The baseline concept was first implemented in a study Ma et al.[[1]](#1) by Murray Christian and Ben Murrell (@Murrellb) for nasal polyp $T_{H}$ cells. `APackOfTheClones` counts the clonotype frequencies for each seurat/umap cluster and produces a circle packing of the frequencies to intuitively represent clonal expansion. The packing for each cluster are then plotted with roughly the same coordinates as the original dimensional reduction and color. Below is an example of `APackOfTheClones` being used on an [example seurat object for scRepertoire](https://drive.google.com/file/d/1_YuRraDyg8UgF3oasjF0-jgPnwox-B24/view?usp=share_link) with its corresponding VDJ library:

<img src="man/figures/example.png" width="100%" align="center" alt="An example APackOfTheClones plot" />

## Installation

`APackOfTheClones` is registered on CRAN. To install the latest stable release, simply run the following

```R
install.packages("APackOfTheClones")
```

For more details on installation methods and information on alternative/development versions, see [```vignette("APackOfTheClones-install")```](https://qile0317.github.io/APackOfTheClones/articles/APackOfTheClones-install.html)

## Usage

The package extends the functionality of `scRepertoire` ***v2*** by working with a seurat object's corresponding T/B cell receptor library. To do this, read the [scRepertoire vignette](borch.dev/uploads/screpertoire). Briefly, an scTCR-seq/scBCR-seq experiment (e.g. from a 10X genomics single cell immune profiling run) should be processed with ```scRepertoire::combineTCR``` / ```scRepertoire::combineBCR``` first. Then, it should be integrated into the corresponding seurat object either with ```scRepertoire::combineExpression```.

To quickly produce the visualization, the ```vizAPOTC(your_combined_seurat_object)``` should give a reasonable visualization. There is an example seurat object included in the package which can be used with ```data("combined_pbmc")```. The following codechunk is an example of how it can be done:

```R
library(Seurat)
library(scRepertoire) # ensure v2 is installed: devtools::install_github("ncborcherding/scRepertoire")
library(APackOfTheClones)

# integrate the contigs with scRepertoire example data - this is identical to "combined_pbmc"
pbmc <- get(data("mini_contig_list", package = "scRepertoire")) %>%
    combineTCR(
        samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L")
    ) %>%
    combineExpression(get(data("scRep_example", package = "scRepertoire")))

# produce the ball-packing plot with the default parameters
vizAPOTC(pbmc)

# there are many parameters to adjust, but most importantly, one can filter for
# subsets of the seurat object easily with keyword arguments corresponding to
# columns in the meta.data dataframe in the seurat object. Here, we filter for
# only the clones of the sample "P17" in a column named "orig.ident" that are only
# in seurat clusters 1, 3, and 4.
vizAPOTC(pbmc, orig.ident = c("P17B", "P17L"), seurat_clusters = c(1, 3, 4))
```

For finer control of the visualization / other parameters, the data for a run can be stored and adjusted with ```RunAPOTC()``` and ```AdjustAPOTC```, then visualized with ```APOTCPlot()```, like so:

```R
# using the same seurat object as before, run and store associated data,
# identifying it with a run id. If not, one will be generated.
pbmc <- RunAPOTC(pbmc, run_id = "default_run_1")

# plot the corresponding data with the run_id. If not provided, will use the latest run.
APOTCPlot(pbmc)

# adjust the sizes of the clones and move the 2nd cluster to (0, 0). Run id this time was
# left empty for no particular reason, but it will use "default_run_1" since its the latest
pbmc <- AdjustAPOTC(
    pbmc,
    new_clone_scale_factor = 0.25,
    relocate_cluster = 2,
    relocation_coord = c(0, 0)
)

# one can also do another run with different parameters and a different id.
pbmc <- RunAPOTC(pbmc, run_id = "run2", try_place = TRUE, repulse = FALSE)

# plot both runs side by side
library(cowplot)
cowplot::plot_grid(
    APOTCPlot(pbmc, run_id = "default_run_1"), APOTCPlot(pbmc, run_id = "run2")
)

```

### Package conventions

Most exported functions are named with `camelCase` with the exception of those that modify Seurat objects having `PascalCase` to mimic their conventions. All function arguments follow `snake_case`.

## Documentation

Comprehensive documentation, vignettes, and a changelog is deployed at <https://qile0317.github.io/APackOfTheClones/>

There are also the following vignettes that should be read in the following order:

- [```vignette("APackOfTheClones")```](https://qile0317.github.io/APackOfTheClones/articles/APackOfTheClones.html)
- [```vignette("APackOfTheClones-runs")```](https://qile0317.github.io/APackOfTheClones/articles/APackOfTheClones-runs.html)
- [```vignette("APackOfTheClones-shared")```](https://qile0317.github.io/APackOfTheClones/articles/APackOfTheClones-shared.html)
- [```vignette("APackOfTheClones-utils")```](https://qile0317.github.io/APackOfTheClones/articles/APackOfTheClones-utils.html)

All exported functions have function level documentation.

## Contributing

Github pull requests from forked branches are more than welcome as it is mostly a solo-project at the moment. For major changes, please open an issue first to discuss what you would like to change. Please also make sure to update tests as appropriate.

An explanation of the algorithms will eventually be made publicly available in a either a pre-print or developer vignette.

## Contact

Qile Yang - qile.yang \[at\] berkeley.edu

## References

<a id="1">[1]</a>
Ma, J., Tibbitt, C. A., Georén, S. K., Christian, M., Murrell, B., Cardell, L. O., Bachert, C., & Coquet, J. M. (2021). Single-cell analysis pinpoints distinct populations of cytotoxic CD4+ T cells and an IL-10+CD109+ TH2 cell population in nasal polyps. Science immunology, 6(62), eabg6356. https://doi.org/10.1126/sciimmunol.abg6356

## Acknowledgements

Thanks for Ben Murrell (@murrelb) at the Karolinska Institute for introducing the idea, implementing julia code, debug support, and giving suggestions. Thanks to Nick Borcherding (@ncborcherding) for providing more insights,  suggestions, and promoting the package.

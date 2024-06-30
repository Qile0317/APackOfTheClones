---
title: 'APackOfTheClones: Visualization of clonal expansion with circle packing'
tags:
  - R
  - dimensional reduction
  - clonal expansion
  - immunology
  - immunoinformatics
  - bioinformatics
  - computational biology
  - seurat
  - single cell
authors:
  - name: Qile Yang
    orcid: 0009-0005-0148-2499
    affiliation: 1
affiliations:
  - name: University of California, Berkeley, Berkeley, CA 94720, United States of America
    index: 1
date: 21 April 2024
bibliography: paper.bib
---

# Summary

T and B lymphocytes exhibit diverse surface receptor repertoires that interact with antigens to communicate and protect the body from foreign invaders. [@den2014activation] Modern methods allows the experimental characterization of the collective identity of the immune cells in tissue samples by sequencing RNA transcripts and cell receptors with scRNA-seq and scTCR-seq / scBCR-seq. Such methods generate high-dimensional data which are analyzed with the aid of dimensional reduction techniques to allow annotation of cell subsets. [@huang2022role]

The R [@r2023r] package APackOfTheClones implements and *extends* a novel method to visualize clonal expansion at a single-cell resolution using circle packing, along with many clonal analysis utilities. Clonotype frequencies for each cell subset is counted and the values are used as radii and packed into one circular cluster with the largest circles near the center, and overlaid onto a the centroids of each cell subset on the corresponding 2D dimensional reduction plot. The package easily integrates into existing analysis pipelines using the *Seurat* [@hao2023dictionary] and *scRepertoire* [@borcherding2020screpertoire] packages. The original implementation was created in the julia language [@benzanson2017julia] by Christian Murray and Ben Murrell and successfully utilized in two immunology papers. [@ma2021single; @stark2022recombinant]

<!-- should those last two sentences be here or in statement of need? -->

# Statement of need

Cells in single-cell experiments are conventionally visualized on a reduction scatterplot where each point represents a cell after all its features have been projected into two dimensions. The points can be colored by different factors to display useful information, including its assigned identity. [@andrews2021tutorial]

Understanding the role of specific cells in various contexts requires understanding of relationships between cell populations and their behaviors. One attribute is the clonal expansion dynamics, which are inferred from downstream analyses. [@adams2020clonal] Overlaying clonal expansion information on a by-clonotype basis onto the reduction scatterplot of the cells of an experiment aids in adding an additional layer of insight, allowing for a swift, intuitive understanding of how clonal dynamics relate to the identified cellular subsets. For example, it can help gauge the presence of hyper expanded clones for each cell type; compare potential changes in frequencies after certain therapeutic treatments, etc.

There is no standardized convention to visualize this single cell level clonotype data on the dimensional reduction. Some of the current approaches include 1. Using a color gradient corresponding to each frequency to highlight each individual point by the clonal expansion, implemented in *scRepertoire* and *scirpy* [sturm2020scirpy] 2. Overlaying a 2D contour where points representing clones with higher frequencies have elevated levels [@andreatta2023tcell] 3. Increasing sizes of points based on clonal frequency, used in figure 2c of @wang2021single.

These approaches capture approximate, global trends but do not provide an exact representation. They lack precision in depicting the true diversity and abundance of clonal populations for every cell subset. From a visual standpoint, This limitation hinders subtle variations in clonal expansion patterns. APackOfTheClones solves this issue by representing exact sizes of each clonotype, in a manner that corresponds exactly to the relevant cell subset. This level of granularity is helpful uncovering hidden patterns, identifying rare clonal populations, and precisely quantifying the impact of therapeutic interventions on immune responses.

APackOfTheClones also offers a suite of methods for visualizing and analyzing single-cell clonal data. Novelty features include functions for highlighting certain clones, and the filtering and visualization of clonotypes shared between subsets by linking circles on the APackOfTheClones clonal expansion plot.

# Results

The main clonal expansion visualization that APackOfTheClones implements is shown in the following figure \autoref{fig:example}, using scRNA-seq + scTCR-seq data from @borcherding2021mapping.

![A single-cell experiment's dimensional reduction projection and its corresponding APackOfTheClones clonal expansion plot. The plot on the right is the projection of all cells of its *Seurat* object based on Uniform Manifold Approximation and Projection (UMAP), colored by unknown identities 1 through 17. On the left is the clonal expansion plot generated by APackOfTheClones on the same data. Each cell identity on the seurat object corresponds to a cluster of circles with a similar geometric placement, and the size of each individual circle is the clonotype frequency of some clonotype within that cell subset. Note that the largest clones are packed near the origins of each cluster to accentuate their difference with the rest of all clonotypes. \label{fig:example}](figures/main_example.png)

<!-- ```R
# the following code was used to generate the plot above - should it be included?

# assume the example `contig_list` and `pbmc` (seurat object) are loaded

contig_list %>%
  combineTCR(
    samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"),
    removeNA = FALSE,
    removeMulti = FALSE,
    filterMulti = FALSE
  ) %>%
  combineExpression(pbmc) %>%
  RunAPOTC() %>%
  AdjustAPOTC(
    new_rad_scale_factor = 0.9,
    nudge_cluster = c(8, 6, 2, 16, 12, 13, 10, 15, 14, 17, 11),
    nudge_vector = list(c(-3,0), c(-1,2), c(-1.25,-1), c(0,-1), c(2, 0), c(1, 0), c(1.3, 0), c(1,0.5), c(2.6, 0), c(2, 0.2), c(3, 0.2))
  ) %>%
  APOTCPlot(
    legend_position = "bottom left",
    legend_sizes = c(1, 50, 200),
    add_legend_background = FALSE,
    use_default_theme = FALSE,
    retain_axis_scales = TRUE
  )
``` -->

The visualization gives the immediate insight that certain cell subsets such as those in cluster four contains many more expanded clones both by quantity and frequency.

The package extends objects and the functionality of the *Seurat* and *scRepertoire* package, and given a correctly processed seurat object of scRNA-seq data that was combined with paired TCR/BCRs, only a few functions need to be used to as little or as much customization of function arguments as needed to produce a `ggplot` object [@wickham2016ggplot2] that fits into the conventional plotting ecosystem of R. Functions are accelerated with a `c++` layer via the *Rcpp* package [@eddelbuettel2011rcpp] to deliver all plots and R objects quickly in time complexity roughly linearly proportional to the number of cells, with the main time bottleneck being the plot display time.

<!--
- should I write about other novelty features - customizing, highlighting, clone links?
- maybe an "implementation details" section about how the user can store "runs" of this plot with different parameters and manually customize them?
-->

# Conclusion

APackOfTheClones offers a fast, and simple interface to produce an intuitive, easily extendible, and *publication ready* visualization of clonal expansion on a cell-by-cell basis, and slots seamlessly into existing analysis pipelines. It can be a useful sub-figure in any immunological/therapeutic study involving single cell omics and immune repertoire to provide an additional degree of understanding for readers and researchers like. However, it should be noted that precise statistical/biological statements about clonal dynamics still obviously require rigorous analysis.

# Acknowledgements

Thanks to Ben Murrell for introducing the idea, implementing the original julia code along with Christian Murray, as well as giving debug support and suggestions. Thanks to Nick Borcherding for providing more insights, suggestions, and promotion.

# References
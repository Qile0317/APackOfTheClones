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
  - name: University of California, Berkeley, Berkeley, CA 94720, USA
    index: 1
date: 10 April 2024
bibliography: paper.bib
---

# Summary

T and B lymphocytes are critical immune cells which express surface receptors that interact with foreign antigens to protect the body. [@den2014activation] There exist a variety of methods to experimentally characterize these receptors at a single-cell resolution, conventionally known as "scTCR-seq / scBCR-seq". Combined with single-cell RNA sequencing ("scRNA-seq"), most cell types and clonotypes (cells proliferated from a common ancestor) can be computationally analyzed to form an immune profile. scRNA-seq generates high-dimensional data where the expression level of each RNA transcript is an additional dimension - which involves using dimensional reduction techniques to (basically section about why and how its commonly used in visualizations for each cell) (TODO refer to paper)

[@huang2022role] The clonal expansion (observed frequency) of each clonotype can be additionally inferred from downstream analysis of any profile and is useful in many studies. [@adams2020clonal] The R [@r2023r] package APackOfTheClones implements and *extends* a novel method to visualize clonal expansion at a single-cell resolution using circle packing, that is easily integratable into existing analysis pipelines using the R packages Seurat [@hao2023dictionary] and scRepertoire [@borcherding2020screpertoire]. Its genesis was done in the julia language [@Julia-2017] by Christian Murray and Ben Murrell and utilized in two papers. [@ma2021single; @stark2022recombinant]

# Statement of need

scRNA-seq naturally generate high-dimensional data for each cell where the expression level of each RNA transcript

One common necessity of many single-cell experiments is to assign identities manually and/or via unsupervised to each cell based on significant transcript features, which frequently involves using some dimensionality reduction technique. This is then conventionally visualized as a reduction scatterplot where each point represents a cell after all its features have been projected into two dimensions. The scatterplot is usually then colored by certain factors and/or annotated on to present more useful information, such as the cell cluster. [@andrews2021tutorial] One such information that is ... (maybe put into summary)

- need to prove why this is truly useful - e.g. can be used to identiy cell type - hyperexpanded cells are more likely to be certian types - can get quick intuition on how certain 

There is no standard method to produce this visualization. Some of the current approaches include 1. Using a color gradient corresponding to each frequency to highlight each individual point by the clonal expansion 2. Overlaying a 2D countour where points representing clones with higher frequencies have elevated levels, created in @andreatta2023tcell 3. Increasing sizes of points based on clonal frequency, used in figure 2c of @wang2021single. TODO state improvement

- https://scirpy.scverse.org/en/latest/generated/scirpy.pl.clonal_expansion.html
- https://www.sciencedirect.com/science/article/pii/S0923753419353852 (evolution of clonotypes too)

APackOfTheClones also makes it simple for anyone to produce the visualization - many single cell pipelines are reliant on the Seurat package, and users only need to add one or two functions in the pipeline that act on the seurat object to produce the plot, with plentiful optional customizations.

# Results

- integrates in seurat pipeline

![A single-cell experiment's UMAP projection and its corresponding APackOfTheClones clonal expansion plot. \label{fig:example}](figures/example.svg)

<!---
Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }
-->

# Acknowledgements

Thanks to Ben Murrell for introducing the idea, implementing the original julia code along with Christian Murray, as well as giving debug support and suggestions. Thanks to Nick Borcherding for providing more insights, suggestions, and promotion.

# References
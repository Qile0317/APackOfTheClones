---
title: 'APackOfTheClones: Novel visualizations of clonal expansion with circle packing'
tags:
  - R
  - dimension reduction
  - clonal expansion
  - immunology
  - immunoinformatics
  - bioinformatics
  - computational biology
authors:
  - name: Qile Yang
    orcid: 0009-0005-0148-2499
    equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: University of California, Berkeley, Berkeley, CA 94720, USA
   index: 1
date: 15 August 2023
bibliography: paper.bib
---

# Summary
APackOfTheClones is an R implementation of a novel method (and utilities) for visualizing clonal expansion of single cell RNA and clonotype data from 10X genomics' Single Cell Immune Profiling, in R. The package can enumerate clone sizes for each cell cluster, and visualize heterogenity by representing each clone as a circle with an area representing relative clone size, and packing all clones into circular clusters with similar coordinates and color to the corresponding cluster in the dimension reduction plot.

APackOfTheClones' approach was originally implemented in the julia language [@Julia-2017] by Christian Murray and Ben Murrell and utillized in two papers. [@ma2021single; stark2022recombinant]

# Statement of need
- state issue to be solved
- show improvement & ease of use over julia implementation

# Method and features
- summary of the method
- integrate_tcr
- getting clonotype frequency
- scRepertoire and Seurat inter operability

# Usage example
- simple pipeline and side-by-side comparison

<!---
Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }
-->

# Performance
- runtime benchmark table & graph

# Discussion
- intuitive
- customizable
- fast
- seamless integration into any analysis pipeline
- simple, no domain knowledge required
- publication ready visualizations
- however, requires just a little extra user work 

# Acknowledgements
Thanks to Ben Murrell at the Karolinska Institute for introducing the idea, implementing julia code along with Christian Murray, debug support, and giving suggestions. Thanks to Nick Borcherding at Washington University in St. Louis for providing more insights, suggestions, and promoting the package.

# References
<!---
For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"
-->
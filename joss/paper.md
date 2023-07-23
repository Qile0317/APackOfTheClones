---
title: 'APackOfTheClones: Visualization of clonal expansion with circle packing'
tags:
  - R
  - C++
  - dimension reduction
  - clonal expansion
  - immunology
  - immunoinformatics
  - bioinformatics
authors:
  - name: Qile Yang
    orcid: 0009-0005-0148-2499
    equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: University of California, Berkeley
   index: 1
date: 15 August 2023
bibliography: paper.bib
---

# Summary

# Statement of need

# method

# result

<!---
Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }
-->

# performance

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Acknowledgements

Thanks for Ben Murrell at the Karolinska Institute for introducing the idea, implementing julia code, debug support, and giving suggestions. Thanks to Nick Borcherding for providing more insights,  suggestions, and promoting the package.

# References
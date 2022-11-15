# work in progress
translation and potential improvement of Ben Murrell's Julia ball packing visualization algorithm in R, and package publication is underway. It was originally implemented in this paper (figure 2d):

```
Ma, J., Tibbitt, C. A., Georén, S. K., Christian, M., Murrell, B., Cardell, L. O., ... & Coquet, J. M. (2021). Single-cell analysis pinpoints distinct populations of cytotoxic CD4+ T cells and an IL-10+ CD109+ TH2 cell population in nasal polyps. Science Immunology, 6(62), eabg6356.
```
<img src="https://github.com/Qile0317/scRNAseq-CircularPacking/blob/main/GoalPlot.png" />

The upcoming package aims to take in TCR/BCR library data and sc-RNAseq data to produce the above "ball-packing" plot that represents clonal expansion. 

The main packing and repulsion algorithm is completed.

Final idea is to allow seamless integration into common scRNAseq pipelines to allow instant visualization with a single function. This will depend on the seurat package for now. There is a script ```scRNAseq_Pipeline``` in main where scRNAseq of the 4 biological samples in the original paper was analyzed. 

# Example of 2 randomly generated clusters
<img src="https://github.com/Qile0317/scRNAseq-CircularPacking/blob/main/Example.png" />

# docs - unfinished
```plot_API(sizes)``` is the main function at the moment. It takes in a list of vectors that are just arrays of numbers representing circle radii. The G argument specifies how hard clusters repell eachother.

### sc-RNAseq integration with the seurat package
As mentioned before, it is underway to get the TCR library data to be able to make the "sizes of clone within a cluster" a reality. The main issue I am having is integrating the scrnaseq and tcr/bcr data. The lack of avaliable packages/knowledge makes it very hard to work with.

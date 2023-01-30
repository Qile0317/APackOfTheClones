#VDJ integration https://www.biostars.org/p/384640/

#load in pre-analyzed seurat object
library(Seurat)
library(tidyverse)
suppressPackageStartupMessages(library(data.table))

#load in the readily-analyzed 10X dataset(This should be done by the user, with UMAP)
# pbmc <- readRDS("C:/Users/lu_41/Desktop/BM_internship/data/10x.rds")
# umap_plt <- DimPlot(object = pbmc, reduction = 'umap')
# umap_plt
# coords <- group_clusters(pbmc) #function from cluster_API.R. In this case i dont think its ac needed.

#load in TCR
tcr <- read.csv("data/fca.csv")
#tcr$barcode <- gsub("-1", "", tcr$barcode)

#takes in a seurat object and tcr dataframe and integrates them
integrate_tcr <- function(seurat_obj, tcr_file) {
  tcr <- as.data.table(tcr_file)

  #cell barcodes are duplicated in rows because
  #the sequencing of TRA and TRB genes creates
  #multiple data points for the same cell.

  # Prepare a progress bar to monitor progress (helpful for large aggregations)
  grpn <- uniqueN(tcr$barcode)
  pb <- txtProgressBar(min = 0, max = grpn, style = 3)

  # Generate a function that will concatenate unique data entries and collapse duplicate rows
  # To do this, I first factorize the data and then get factor levels as unique data points
  # Then I paste these  data points together separated with "__" to access later on if needed

  data_concater <- function(x){
    x <- levels(factor(x))
    paste(x, collapse = "__")
  }

  # This code applies data_concater function per  barcodes to create a
  # concatenated string with  the information we want to keep
  tcr_collapsed <- tcr[, {setTxtProgressBar(pb,.GRP); lapply(.SD, data_concater)} , by=barcode]

  #assign rownames for integration
  rownames(tcr_collapsed) <- tcr_collapsed$barcode
  pbmc <- AddMetaData(pbmc, metadata = tcr_collapsed)
  return(pbmc)
}

# pbmc <- integrate_tcr(pbmc, tcr)

#thats it. Lets see if theres alot of NAs
# head(pbmc@assays[["RNA"]]@data@Dimnames[[2]])
# head(pbmc@meta.data[["barcode"]])
#simple looping
percentNA <- function(d) {
  len <- 0
  na <- 0
  for (i in d) {
    len <- len + 1
    if (is.na(i)) {na <- na + 1}
  }
  return(na/len)
}
# percentNA(pbmc@meta.data[["barcode"]]) #69% ?!?!?!? ALthough idk how many T cells are in here.

# To calculate the clone sizes, can use the TCR information 
# to identify the unique TCR sequences for each cell 
# and then count the number of cells with each TCR sequence.

# takes in seurat obj, calculates clone sizes (radii) and make the list of radii for the plot_API
# in the future can add color and label easily.
# also need to get UMAP cluster centroids
get_clone_sizes <- function(integrated_seurat_obj) { #needs to have UMAP ran already
  centroid_list <- find_centroids(group_clusters(integrated_seurat_obj))
  print("unfinished")
}


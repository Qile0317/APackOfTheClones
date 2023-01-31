#VDJ integration https://www.biostars.org/p/384640/

#load in pre-analyzed seurat object
library(Seurat)
suppressPackageStartupMessages(library(tidyverse))
library(hash) # conflicts w "copy" from data.table
suppressPackageStartupMessages(library(data.table))

#shortcut function
get_umap <- function(seurat_obj) {
  return(DimPlot(object = seurat_obj, reduction = 'umap'))
}

#takes in a seurat object and tcr dataframe and integrates them into a NEW object
integrate_tcr <- function(seurat_obj, tcr_file) {
  new_seurat_obj <- seurat_obj
  tcr <- as.data.table(tcr_file)

  #cell barcodes are duplicated in rows because
  #the sequencing of TRA and TRB genes creates
  #multiple data points for the same cell.

  # Prepare a progress bar to monitor progress (helpful for large aggregations)
  grpn <- uniqueN(tcr$barcode)
  pb <- txtProgressBar(min = 0, max = grpn, style = 3)

  # Generate a function that will concatenate unique data entries and collapse duplicate rows
  # To do this, I first factorize the data and then get factor levels as unique data points
  # Then data points are pasted together separated with "__" to access later on if needed

  data_concater <- function(x){
    x <- levels(factor(x))
    paste(x, collapse = "__")
  }

  # This code applies data_concater function per  barcodes to create a
  # concatenated string with  the information we want to keep
  tcr_collapsed <- tcr[, {setTxtProgressBar(pb,.GRP); lapply(.SD, data_concater)} , by=barcode]

  #assign rownames for integration
  rownames(tcr_collapsed) <- tcr_collapsed$barcode
  new_seurat_obj <- AddMetaData(new_seurat_obj, metadata = tcr_collapsed)
  return(new_seurat_obj)
}

# test on real data and integrate
pbmc <- readRDS("C:/Users/lu_41/Desktop/BM_internship/scballpack/data/10x.rds")
umap_plt <- get_umap(pbmc)
tcr <- read.csv("data/fca.csv")
pbmc <- integrate_tcr(pbmc, tcr)

# view a few parts
# head(pbmc@assays[["RNA"]]@data@Dimnames[[2]])
# head(pbmc@meta.data[["barcode"]])

percent_na <- function(seurat_obj) {
  d <- seurat_obj@meta.data[["barcode"]]
  len <- 0
  na <- 0
  for (i in d) {
    len <- len + 1
    if (is.na(i)) {
      na <- na + 1
    }
  }
  return(na / len)
}

# percent_na(pbmc) #69% ?!?!?!? ALthough idk how many T cells are in here.

# To calculate the clone sizes, can use the TCR information
# to identify the unique TCR sequences for each cell
# and then count the number of cells with each TCR sequence.

# takes in seurat obj, calculates clone sizes (radii) and make the list of radii for the plot_API
# in the future can add color and label easily.

#split seurat obj into UMAP coords and clusters into df
group_clusters <- function(SeuratObj){
  return(data.frame(
    SeuratObj@reductions[["umap"]]@cell.embeddings,
    clusters = SeuratObj$seurat_clusters))
}

# generate a hashmaps of barcode:cluster pairs for clone sizes (is putting in the title equivalent to hashing?)
#first a simple custom function
dict <- function(keys, values) {
  output <- hash()
  for (i in 1:length(keys)) {
    output[[keys[i]]] <- values[i]
  }
  return(output)
}

#generate the hashmap for get_clone sizes
gen_cluster_dict <- function(seurat_obj) {
  barcodes <- seurat_obj@meta.data[["barcode"]]
  clusters <- as.numeric(seurat_obj$seurat_clusters)
  return(dict(barcodes,clusters))
}

#processes the intergers in a string separated by "__" and sum them
sum_read_str <- function(numstr) {
  return(sum(as.numeric(unlist(strsplit(numstr, "__")))))
}

# wrapper to get the number of identified clusters in the UMAP:
umap_cluster_count <- function(seurat_obj) {
  return(length(levels(seurat_obj@meta.data[["seurat_clusters"]])))
}

#generate the data for the plot_API(), however it just counts the reads for each barcode... its not the clone size. how do I get that??
# alot of speed improvements can be made from more hashing

get_clone_sizes <- function(integrated_seurat_obj, normalize = TRUE) {

  #initialize outputs
  radii <- replicate(umap_cluster_count(integrated_seurat_obj), c())
  min_read_count <- Inf

  #label_radii_index_dict <- dict() #if in the future i wanna assign label to cluster i can do this

  #work in progress:
  label_list <- list()
  color_list <- list()

  #variables for iteration
  cluster_dict <- gen_cluster_dict(integrated_seurat_obj)
  meta_data_barcodes <- integrated_seurat_obj@meta.data[["barcode"]]
  meta_data_readcount <- integrated_seurat_obj@meta.data[["reads"]]

  #iteration
  for (i in 1:length(meta_data_barcodes)) {
    curr_barcode <- meta_data_barcodes[i]
    if (!is.na(curr_barcode)) {
      current_cluster <- cluster_dict[[curr_barcode]] # assumes they all exist in the dict
      current_reads <- sum_read_str(meta_data_readcount[i])

      if (normalize) {
        min_read_count <- min(min_read_count, current_reads)
      }

      radii[[current_cluster]] <- c(radii[[current_cluster]],current_reads)
    }
  }

  if(!normalize){return(radii)}

  # normalization to 1
  norm_factor = 1/min_read_count
  for (i in 1:length(radii)) {
    for (j in 1:length(radii[[i]])) {
      if (!is.null(radii[[i]][j])) {
        radii[[i]][j] <- radii[[i]][j] * norm_factor # doesn't work
      }
    }
  }
  return(radii)
}

#get centroids (wrapper)
get_cluster_centroids <- function(seurat_obj) {
  return(find_centroids(group_clusters(seurat_obj))) # indexing is a bit off but no big deal
}

#these can be passed to plot_APi
pbmc_size_list <- get_clone_sizes(pbmc,normalize=FALSE)
pbmc_norm_size_list <- get_clone_sizes(pbmc)
pbmc_centroid_list <- get_cluster_centroids(pbmc)


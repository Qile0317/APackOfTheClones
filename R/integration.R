# creating simpler (tho prob unessecary) dictionary syntax
dict <- function(keys, values) {
  output <- hash()
  for (i in 1:length(keys)) {
    output[[keys[i]]] <- values[i]
  }
  return(output)
}

#takes in a seurat object and tcr dataframe and integrates them into a NEW object
integrate_tcr <- function(seurat_obj, tcr_file) {
  new_seurat_obj <- seurat_obj
  tcr <- as.data.table(tcr_file)

  #cell barcodes are duplicated in rows because
  #the sequencing of TRA and TRB genes creates
  #multiple data points for the same cell.

  # Prepare a progress bar to monitor progress (helpful for large aggregations)
  message("integrating TCR library into seurat object")
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
  tcr_collapsed <- tcr[, {setTxtProgressBar(pb,.GRP); lapply(.SD, data_concater)} , by = barcode]

  #assign rownames for integration
  rownames(tcr_collapsed) <- tcr_collapsed$barcode
  new_seurat_obj <- AddMetaData(new_seurat_obj, metadata = tcr_collapsed)

  return(new_seurat_obj)
}

#optional function to check how many didn't have matches
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
  return(100 * na / len)
}

# wrapper to get the number of identified clusters in the UMAP:
count_umap_clusters <- function(seurat_obj) {
  return(length(levels(seurat_obj@meta.data[["seurat_clusters"]])))
}

#take in a list of dicts and return a list of vectors of only values
dict_list_vals <- function(dict_list, scale_factor = 0.001) {

  output <- list()

  # multiply by scale factor for sapply
  multiply_by_factor <- function(x) {return(x * scale_factor)}

  for (i in 1:length(dict_list)) {
    curr_hashmap <- dict_list[[i]]
    if (is.null(curr_hashmap)) {
      output[[i]] <- NULL
    }else{
      output[[i]] <- sapply(unname(values(curr_hashmap)), multiply_by_factor)
    }
  }
  return(output)
}

# counting clonotypes to make the sizes of clones within cluster data. Needs testing
get_clone_sizes <- function(integrated_seurat_obj, scale_factor = 0.001) {

  clonotype_dict_list <- replicate(count_umap_clusters(integrated_seurat_obj), hash())
  min_read_count <- Inf

  #work in progress
  #label_radii_index_dict <- dict()
  #label_list <- list()
  #color_list <- list()

  #variables for iteration
  barcodes <- integrated_seurat_obj@meta.data[["barcode"]]
  clusters <- integrated_seurat_obj@meta.data[["seurat_clusters"]]
  clonotype_ids <- integrated_seurat_obj@meta.data[["raw_clonotype_id"]]

  #iteration, and construct a lsit of clonotype count dicts
  for (i in 1:length(clusters)) {
    if (!is.na(barcodes[i])) {

      curr_clonotype_id <- clonotype_ids[i] # redundant?
      curr_cluster <- clusters[i] # assumes its numeric
      curr_dict_val <- clonotype_dict_list[[curr_cluster]][[curr_clonotype_id]]

      #add to dictionary
      if (!is.null(curr_dict_val)){
        clonotype_dict_list[[curr_cluster]][[curr_clonotype_id]] <- curr_dict_val + 1
      }
      clonotype_dict_list[[curr_cluster]][[curr_clonotype_id]] <- 1
    }
  }
  #convert the list of dicts into list of radii
  return(dict_list_vals(clonotype_dict_list, scale_factor))
}

#get centroids (wrapper) # indexing is a bit off but no big deal for the function
get_cluster_centroids <- function(seurat_obj) {
  return(find_centroids(data.frame(
    seurat_obj@reductions[["umap"]]@cell.embeddings,
    clusters = seurat_obj$seurat_clusters)
    ))
}

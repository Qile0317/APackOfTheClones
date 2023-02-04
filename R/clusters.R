#script for functions to deal with centroids and cluster coords

#centroid finder for a whole dataframe. returns dataframe.
find_centroids <- function(df, return_mode = "list") { # or "df"
  cll <- split(df, factor(df[, 3])) #the last cluster column becomes redundant
  l <- length(cll)

  nameset <- rep(c(""), times = l)
  xset <- rep(c(0), times = l)
  yset <- xset

  for (i in 1:l){
    nameset[i] <- cll[[i]][,3][1] #i dont think this works lol
    xset[i] <- sum(cll[[i]][, 1])/length(cll[[i]][, 1])
    yset[i] <- sum(cll[[i]][, 2])/length(cll[[i]][, 2])
  }
  if (return_mode != "list") { # return df, although its not rlly needed...
    return(data.frame(cluster = nameset, x = xset, y = yset))
  }
  #return list (theres definetely a smarter way to do this)
  list_output <- list()
  for (i in 1:l) {
    list_output[[nameset[i]]] <- c(xset[i],yset[i])
  }
  return(list_output)
}

#get centroids from seurat obj # indexing is a bit off but no big deal for the function
get_cluster_centroids <- function(seurat_obj) {
  return(find_centroids(data.frame(
    seurat_obj@reductions[["umap"]]@cell.embeddings,
    clusters = seurat_obj$seurat_clusters)
  ))
}

#td <- data.frame(pbmc@reductions[["umap"]]@cell.embeddings,clusters = pbmc$seurat_clusters)

#transform coordinates of a cluster from c(0, 0) to its own new centroid, or MOVE to new coord from current
trans_coord <- function(cluster, new_coord = NULL) {
  len <- length(cluster[[1]]) # x,y,r have the same length

  if (!is.null(new_coord)) {
    dx <- new_coord[1] - cluster[[4]][1]
    dy <- new_coord[2] - cluster[[4]][2]
    cluster[[4]] <- new_coord

  }else {
    dx <- cluster[[4]][1]
    dy <- cluster[[4]][2]
  }

  ansc <- list("x" = rep(0, len),
               "y" = rep(0, len),
               "rad" = cluster[[3]],
               "centroid" = cluster[[4]],
               "clRad" = cluster[[5]])

  for(i in 1:len) {
    ansc[[1]][i] <- cluster[[1]][i] + dx
    ansc[[2]][i] <- cluster[[2]][i] + dy
  }

  return(ansc)
}

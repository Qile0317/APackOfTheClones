#script for function to manipulate and manage all clusters. Once again user doesn't need this script

#centroid finder for a whole dataframe. returns dataframe.

# 1st row must be UMAP_1 (x)
# 2nd row must be UMAP_2 (y)
# 3rd row of the dataframe must be clusters
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

#transform coordinates of a cluster by its own centroid
trans_coord <- function(cluster) {
  ansc <- list(x = c(), y = c(),
               rad = cluster[[3]],
               centroid = cluster[[4]],
               clRad = cluster[[5]])

  for(i in 1:length(cluster[[1]])) {
    ansc[[1]][i] <- cluster[[1]][i] + cluster[[4]][1] # x of centroid
    ansc[[2]][i] <- cluster[[2]][i] + cluster[[4]][2] # y of centroid
  }
  return(ansc)
}
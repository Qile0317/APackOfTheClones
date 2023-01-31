#script for function to manipulate and manage all clusters. Once again user doesn't need this script

#centroid finder for a whole dataframe. returns dataframe.
#it also needs to be able to incorporate color. the labelling is also off somehow.
#Run group_clusters before imputting df.

#the FIRST row must be UMAP_1 (x)
#the SECOND row must be UMAP_2 (y)
#the THIRD row of the dataframe must be clusters
find_centroids <- function(df, return_mode = "list"){ # or "df"
  cll <- split(df, factor(df[,3])) #the last cluster column becomes redundant
  l <- length(cll)
  nameset <- rep(c(""),times=l)
  xset <- rep(c(0),times=l)
  yset <- xset
  for (i in 1:l){
    nameset[i] <- cll[[i]][,3][1] #i dont think this works lol
    xset[i] <- sum(cll[[i]][,1])/length(cll[[i]][,1])
    yset[i] <- sum(cll[[i]][,2])/length(cll[[i]][,2])
  }
  if (return_mode != "list") { # return df, although its not rlly needed...
    return(data.frame(cluster = nameset, x = xset, y = yset))
  }
  #return list (theres definetely a smarter way to do this)
  list_output <- list()
  for (i in 1:l){
    list_output[[nameset[i]]] <- c(xset[i],yset[i])
  }
  return(list_output)
}

#t <- find_centroids(clsdf)
#ggplot(t,aes(x=x,y=y,color=cluster)) + geom_point() + umapPlot
#looks right!

#transform coordinates of a cluster by the transvec, c(x,y)
trans_coord <- function(cluster, transvec){
  ansc <- list(x=c(),y=c(),rad=cluster[[3]],centroid=cluster[[4]]+transvec,clRad=cluster[[5]])
  for(i in 1:length(cluster[[1]])){
    ansc[[1]][i] <- cluster[[1]][i] + transvec[1]
    ansc[[2]][i] <- cluster[[2]][i] + transvec[2]
    }
  return(ansc)
  }


trans_cluster <- function(cluster_list){
  ansc <- list()
  for(i in 1:length(cluster_list)){
    ansc[[i]]<-trans_coord(cluster_list[[i]],cluster_list[[i]][[4]])
  }
  return(ansc)
}

#script for function to manipulate and manage all clusters. Once again user doesn't need this script

#split seurat obj into UMAP coords and clusters
group_clusters <- function(SeuratObj){
  return(data.frame(
    SeuratObj@reductions[["umap"]]@cell.embeddings,
    clusters = SeuratObj$seurat_clusters))
}

#centroid finder for a whole dataframe. returns dataframe.
#it also needs to be able to incorporate color. the labelling is also off somehow.
#Run group_clusters before imputting df.

#the FIRST row must be UMAP_1 (x)
#the SECOND row must be UMAP_2 (y)
#the THIRD row of the dataframe must be clusters
find_centroids <- function(df){
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
  return(data.frame(cluster=nameset,x=xset,y=yset))
}

#t <- find_centroids(clsdf)
#ggplot(t,aes(x=x,y=y,color=cluster)) + geom_point() + umapPlot
#looks right!

#functions for moving cluster to NEWLY APPOINTED centroid
#this is for 1 that creates a new list which is inefficient but R is very dumb with mutation. might try in the future with the Lepl() function
trans_coord <- function(cluster, new_centroid = c("none","none")){
  ansc <- list(x=c(),y=c(),rad=cluster[[3]],centroid=cluster[[4]],clRad=cluster[[5]])
  #cluster <- deparse(substitute(cluster))
  if(!identical(c("none","none"),new_centroid)){ansc[[4]] <- new_centroid} #mutate_list(deparse(substitute(cluster)),4,new_centroid)
  for(i in 1:length(cluster[[1]])){
    #mutate_list_el(deparse(substitute(cluster)),1,i,cluster[[1]][i] + cluster$centroid[1])
    ansc[[1]][i] <- cluster[[1]][i] + cluster[[4]][1]
    #mutate_list_el(deparse(substitute(cluster)),2,i,cluster[[2]][i] + cluster$centroid[2])
    ansc[[2]][i] <- cluster[[2]][i] + cluster[[4]][2]
    }
  return(ansc)
  }

trans_cluster <- function(cluster_list){
  ansc <- list()
  for(i in 1:length(cluster_list)){
    ansc[[i]]<-trans_coord(cluster_list[[i]])
  }
  return(ansc)
}

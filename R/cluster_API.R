#script for function to manipulate and manage all clusters

#The packing algo is basically done. Now to manage centroid locations of both clusters is easy, coords can be transformed. 
#split seurat obj into UMAP coords and 
group_clusters <- function(SeuratObj){
  return(data.frame(
    SeuratObj@reductions[["umap"]]@cell.embeddings,
    clusters = SeuratObj$seurat_clusters))
}

#clsdf <- group_clusters(pbmc)
#ggplot(data=clsdf,aes(x=UMAP_1,y=UMAP_2,color=clusters)) + geom_point()

#centroid finder for a whole dataframe. returns dataframe. I think it might be bugged.
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

# result plotting function. clusters is a list of lists transformed into a dataframe, which are clusters. 
# A cluster list includes $x, $y, $rad, $centroid.
#the clusters imput is a dataframe. 
plot_clusters <- function(clusters, n=360, linetype="blank", #linewidth=1, #linewidth doesnt work lol.
                          title = "Sizes of clones within each cluster",
                          haslegend=TRUE, void=TRUE,
                          origin=FALSE){
  if(!origin){
    p1 <- ggplot() + geom_circle(data = clusters, mapping = aes(
      x0 = x, y0 = y, r=r, fill=label),  n=n, linetype=linetype) + #higher n is basically higher resolution
      labs(title = title) + coord_fixed()  
    #geom_text(data = clusters, aes(x,y, label = label)) #this should only be near a cluster. can make simple function to put it on bottom right.
    if(void){p1 <- p1 + theme_void()}
    if(haslegend){p1 <- p1 + theme(legend.position="none")}
    p1
  }else{
    p1 <- ggplot(clusters,mapping=aes(x,y)) + geom_point() + labs(title = title) + coord_fixed()
    if(void){p1 <- p1 + theme_void()}
    p1
  }
}

#now: be able to group by color and have an actual good title
#repulsion is coming

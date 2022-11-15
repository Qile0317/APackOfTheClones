#main plotting function

# result plotting function. clusters is a list of lists transformed into a dataframe, which are clusters.
# A cluster list includes $x, $y, $rad, $centroid.
#the clusters imput is a dataframe. #move into seperate script

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

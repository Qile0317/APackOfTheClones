#main plotting function

# result plotting function. clusters is a list of lists transformed into a dataframe, which are clusters.
# A cluster list includes $x, $y, $rad, $centroid.
#the clusters imput is a dataframe. #move into seperate script

plot_clusters <- function(clusters, n=360, linetype="blank", #linewidth=1, #linewidth doesnt work lol.
                          title = "Sizes of clones within each cluster",
                          haslegend=TRUE, void=TRUE,
                          origin=FALSE){ #, label=TRUE (lavels each individual circle, yikes)
  if(!origin){
    p1 <- ggplot() + geom_circle(data = clusters, mapping = aes(
      x0 = x, y0 = y, r=r, fill= .data[["label"]]), # not sure if .data shoulf be here
      n=n, linetype=linetype) +
      labs(title = title) + coord_fixed()
    #if(label){p1 <- p1 + geom_text(data = clusters, aes(x,y, label = .data[["label"]]))} #this should only be near a cluster. can make simple function to put it on bottom right.
    if(void){p1 <- p1 + theme_void()}
    if(haslegend){p1 <- p1 + theme(legend.position="none")}
  }else{
    p1 <- ggplot(clusters,mapping=aes(x,y)) + geom_point() + labs(title = title) + coord_fixed()
    if(void){p1 <- p1 + theme_void()}
  }
  return(p1)
}

#now: be able to group by color and have an actual good title
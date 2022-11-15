#The API functions for generating the ggplot, with repulsion

#API for plotting size vectors . Integrating the data is work in progress.
plot_API <- function(sizes, # list of size vectors,[[1]] c(a,b,..) and centroids [[2]] c(x,y)
                     ORDER = TRUE, try_place = TRUE, progbar = TRUE, # packing
                     repulse=TRUE, thr=1, G=0.05, max_iter=100, #repulsion
                     n=360, linetype="blank",
                     ggtitle = "Sizes of clones within each cluster",
                     haslegend=TRUE, void=TRUE, origin=FALSE){ #prev 3 lines were for plotting function
  ans <- list()

  #circle layout
  for(i in 1:length(sizes)){
    if(progbar){print(paste("packing cluster", as.character(i))," \n")}
    ans[[i]] <- circle_layout(sizes[[i]],
                              ORDER = ORDER,
                              try_place = try_place,
                              progbar = progbar)
  }
  #repulsion
  if(repulse){ans <- repulse_cluster(ans, thr = thr, G = G, max_iter = max_iter)}
  print("repulse done")
  #joining list into df
  ans <- df_full_join(ans)

  #plotting
  plot_clusters(ans, n = n, linetype = linetype, title = title,
                haslegend = haslegend, void = void, origin = origin)
}

#the following is unfinished.
#' @export
scballpack <- function(seuratobj, repulse=TRUE, thr=3, n=360){
  #functions to get the data from the seuratobj, return lists of centroid and areas. return in a list called mainli
  for(i in 1:length(mainli)){
    mainli[[i]] <- circle_layout(mainli[[i]])
  }
  #transform
  #repulse
  maindf <- data.frame()
  plt <-
  return(plot_clusters(maindf))
}

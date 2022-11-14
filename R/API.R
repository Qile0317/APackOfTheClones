#The main API function for generating the ggplot, with repulsion

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
  return(plot_clusters(maindf),n=n)
}

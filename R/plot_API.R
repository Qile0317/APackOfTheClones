#The API functions for generating the ggplot

library(dplyr)
library(data.table)
library(utils)
library(ggplot2)
library(ggforce)

#full join a list of lists of (x,y,r) vectors into a dataframe with generated labels.
df_full_join <- function(clstr_list) {
  df <- data.frame(label = character(0),
                   x = numeric(0),
                   y = numeric(0),
                   r = numeric(0))
  
  seurat_cluster_index <- 0 # impportant to skip empty cluster indicies for insert_colors, need testing
  for(i in 1:length(clstr_list)){
    if (!any(is.na(clstr_list[[i]]))) {
      df <- dplyr::full_join(
        df, data.frame("label" = rep(paste("cluster", seurat_cluster_index),
                                     length(clstr_list[[i]][["x"]])),
                       "x" = clstr_list[[i]][["x"]],
                       "y" = clstr_list[[i]][["y"]],
                       "r" = clstr_list[[i]][["rad"]]),
        by = dplyr::join_by("label", "x", "y", "r"))
    }
    seurat_cluster_index <- seurat_cluster_index + 1
  }
  df
}

#main plotting function

# result plotting function. clusters is a list of clusterlists TRANSFORM into a dataframe, which are clusters.
# A cluster list includes $x, $y, $rad, $centroid.
#the clusters imput is a dataframe. #move into seperate script

plot_clusters <- function(clusters, n = 360, linetype ="blank", #linewidth=1, #linewidth doesnt work lol.
                          title = "Sizes of clones within each cluster",
                          haslegend = FALSE, void = TRUE,
                          origin = FALSE){ #, label=TRUE (labels each individual circle, yikes)
  if (!origin) {
    p1 <- ggplot2::ggplot(data = clusters) +
      ggforce::geom_circle(aes_string(
        x0 = "x", y0 = "y", r = "r", fill = "color"),
        n = n, linetype = linetype) +
      ggplot2::scale_fill_identity() +
      
      ggplot2::labs(title = title) +
      ggplot2::coord_fixed()
  }else {
    p1 <- ggplot2::ggplot(
      clusters, 
      mapping = aes_string(x = "x", y = "y")
      ) + 
      ggplot2::geom_point() +
      ggplot2::labs(title = title) +
      ggplot2::coord_fixed()
  }
  
  if (void) {
    p1 <- p1 + ggplot2::theme_void()
  }
  
  return(p1 + ggplot2::theme(legend.position = "none"))
}

#now: be able to group by color...

#API for plotting size vectors .
plot_API <- function(sizes, # list of size vectors,[[1]] c(a,b,..)
                     centroids, # centroids of size vectors [[1]] c(x,y)
                     num_clusters,
                     rad_decrease = 1,
                     ORDER = TRUE,
                     try_place = TRUE,
                     progbar = TRUE, 
                     repulse = FALSE,
                     thr = 1, G = 0.05, 
                     max_repulsion_iter = 100, #repulsion parameters
                     n = 360, linetype = "blank",
                     plot_title = "Sizes of clones within each cluster",
                     haslegend = TRUE,
                     void = TRUE,
                     origin = FALSE){
  ans <- list()
  
  #circle layout
  for(i in 1:length(sizes)){
    if (is.null(sizes[[i]])) {
      ans[[i]] <- NA # important!
    }else{
      if(progbar){
        message(paste("\npacking cluster", as.character(i-1)))
      }
      ans[[i]] <- circle_layout(
        sizes[[i]],
        centroid = centroids[[i]],
        rad_decrease = rad_decrease,
        ORDER = ORDER,
        try_place = try_place,
        progbar = progbar
      )
    }
  }
  
  if (repulse) {
    if(progbar){message(paste("\nrepulsing clusters | max iterations =", max_repulsion_iter))}
    ans <- repulse_cluster(ans, thr = thr, G = G, max_iter = max_repulsion_iter, verbose = progbar)
  }
  
  #joining list into df and color for plotting. in future make customizable
  ans <- df_full_join(ans)
  ans <- insert_colors(ans, num_clusters)
  
  #plotting
  plt <- plot_clusters(ans, n = n, linetype = linetype, title = plot_title,
                       haslegend = haslegend, void = void, origin = origin)
  plt
}

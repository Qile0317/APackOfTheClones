#The API functions for generating the ggplot

library(dplyr)
library(hash)
library(data.table)
library(utils)
library(ggplot2)
library(ggforce)

# shortcut to get the umap plot
get_umap <- function(seurat_obj) {
  return(Seurat::DimPlot(object = seurat_obj,
                         reduction = 'umap'))
}

#full join a list of lists of (x,y,r) vectors into a dataframe with generated labels.
df_full_join <- function(clstr_list) {
  df <- data.frame(label = character(0),
                   x = numeric(0),
                   y = numeric(0),
                   r = numeric(0))

  for(i in 1:length(clstr_list)){
    if (!is.null(clstr_list[[i]])) {

      df <- dplyr::full_join(
        df, data.frame("label" = paste("cluster", as.character(i-1)),
                       "x" = clstr_list[[i]]$x,
                       "y" = clstr_list[[i]]$y,
                       "r"=clstr_list[[i]]$rad),
        by = dplyr::join_by("label", "x", "y", "r"))
    }
  }
  return(df)
}

#main plotting function

# result plotting function. clusters is a list of clusterlists TRANSFORM into a dataframe, which are clusters.
# A cluster list includes $x, $y, $rad, $centroid.
#the clusters imput is a dataframe. #move into seperate script

plot_clusters <- function(clusters, n = 360, linetype="blank", #linewidth=1, #linewidth doesnt work lol.
                          title = "Sizes of clones within each cluster",
                          haslegend = TRUE, void = TRUE,
                          origin = FALSE){ #, label=TRUE (lavels each individual circle, yikes)
  if(!origin){
    p1 <- ggplot2::ggplot() + ggforce::geom_circle(
      data = clusters, mapping = ggplot2::aes(
        x0 = ggplot2::.data[["x"]],
        y0 = ggplot2::.data[["y"]],
        r = ggplot2::.data[["r"]],
        fill = ggplot2::.data[["label"]]),

      n = n, linetype = linetype) +
      ggplot2::labs(title = title) +
      ggplot2::coord_fixed()

    if(void){p1 <- p1 + ggplot2::theme_void()}

    if(!haslegend){
      p1 <- p1 + ggplot2::theme(legend.position="none")
    }else{
      #p1 <- p1 + guide_legend(title = NULL, )
    }
  }else{
    p1 <- ggplot2::ggplot(
      clusters,
      mapping = ggplot2::aes(
        ggplot2::.data[["x"]],
        ggplot2::.data[["y"]])
      ) +
      ggplot2::geom_point() + ggplot2::labs(title = title) +
      ggplot2::coord_fixed()

    if(void){p1 <- p1 + ggplot2::theme_void()}
  }
  return(p1)
}

#now: be able to group by color...

#API for plotting size vectors .
#I could also export this
plot_API <- function(sizes, # list of size vectors,[[1]] c(a,b,..)
                     centroids, # centroids of size vectors [[1]] c(x,y)
                     rad_decrease = 0,
                     ORDER = TRUE,
                     try_place = TRUE,
                     progbar = TRUE, # packing
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

    if (length(sizes[[i]]) == 0) {
      ans[[i]] <- NULL

    }else{
      if(progbar){
        message("")
        message(paste("packing cluster", as.character(i)))
      }

      ans[[i]] <- circle_layout(sizes[[i]],
                                centroid = centroids[[i]],
                                rad_decrease = rad_decrease,
                                ORDER = ORDER,
                                try_place = try_place,
                                progbar = progbar)
    }
  }

  #repulsion, not sure how it handles nulls
  if (repulse) {
    if(progbar){message("repulsing clusters")}
    ans <- repulse_cluster(ans, thr = thr, G = G, max_iter = max_repulsion_iter)
  }

  #joining list into df for plotting
  ans <- df_full_join(ans)

  #return(ans) # debug. but everything works untill this point

  #plotting
  ans <- plot_clusters(ans, n = n, linetype = linetype, title = plot_title,
                       haslegend = haslegend, void = void, origin = origin)
  return(ans)
}

# change the axis scales to fit the original plot approximately.
retain_scale <- function(seurat_obj, ball_pack_plt) {

  umap_plt <- get_umap(seurat_obj)

  umap_xr <- ggplot2::ggplot_build(umap_plt)$layout$panel_scales_x[[1]]$range$range
  umap_yr <- ggplot2::ggplot_build(umap_plt)$layout$panel_scales_y[[1]]$range$range

  ball_pack_xr <- ggplot2::ggplot_build(ball_pack_plt)$layout$panel_scales_x[[1]]$range$range
  ball_pack_yr <- ggplot2::ggplot_build(ball_pack_plt)$layout$panel_scales_y[[1]]$range$range

  # new ranges
  min_xr <- max(ball_pack_xr[1], umap_xr[1])
  max_xr <- max(ball_pack_xr[2], umap_xr[2])

  min_yr <- max(ball_pack_yr[1], umap_yr[1])
  max_yr <- max(ball_pack_yr[2], umap_yr[2])

  return(ball_pack_plt +
           ggplot2::xlim(min_xr, max_xr) +
           ggplot2::ylim(min_yr, max_yr)
         )
}

# A more advanced version could multiply axses by a small amount to retain ratios exactly
# also this isnt perfect, in my own testcase 1 row was removed

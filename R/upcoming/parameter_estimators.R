# script to estimate some okay-looking parameters for the plot

#samples <- c(mini_seurat_obj, pbmc)
#cell_count <- c(80, 2500)
#desirable_factor <- c(1,0.2)
#plot(cell_count, desirable_factor)
#with these parameters the line is y = -0.0003305785x + 1
#a more complicated model can use clonotype distribution as third param

estimate_clone_scale_factor <- function(seurat_obj) {
  # get valid cell count
  # return result made by regression line
}

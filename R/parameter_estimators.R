#cell_count <- c(80, 2500)
#desirable_factor <- c(1,0.2)
#plot(cell_count, desirable_factor)
#with these parameters the line is y = -0.0003305785x + 1
#a more complicated model can use clonotype distribution as third param

estimate_clone_scale_factor <- function(seurat_obj, verbose = FALSE) {
	num <- (-0.003305785 * count_valid_barcodes(seurat_obj)) + 1
	if (verbose) {
		message(paste("setting clone_scale_factor to", num))
	}
	num
}

convert_rad_decrease <- function(seurat_obj) {

}

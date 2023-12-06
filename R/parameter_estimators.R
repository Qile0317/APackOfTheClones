#cell_count <- c(80, 2500)
#desirable_factor <- c(1,0.2)
#plot(cell_count, desirable_factor)
#with these parameters the line is y = -0.0003305785x + 1
#a more complicated model can use clonotype distribution as third param

estimate_clone_scale_factor <- function(seurat_obj, verbose = FALSE) {
	num <- max((-0.0003305785 * count_valid_barcodes(seurat_obj)) + 1, 0.1)
	if (verbose) {
		message(paste("\nsetting clone_scale_factor to", num))
	}
	num
}

# for readability - assumes rad_scale_factor between 0 & 1
convert_to_rad_decrease <- function(rad_scale_factor, clone_scale_factor) {
	clone_scale_factor * (1 - rad_scale_factor)
}

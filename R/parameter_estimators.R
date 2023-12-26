# # a super simple regression of manually determined clone scale based on cell
# # count a much more complicated model can use other facts about the seurat
# # obj to improve how visually plesant is it. Some overlap factor could 
# # probably be eestimated with the raw clone counts
# cell_count <- c(80, 365, 2500)
# desirable_factor <- c(1, 0.3, 0.2)
# plot(cell_count, desirable_factor)

estimate_clone_scale_factor <- function(seurat_obj, clonecall) {
	num_clones <- count_clones(seurat_obj, clonecall)

	if (num_clones <= 365) {
		approx_clone_scale_factor <- (-0.002456 * num_clones) + 1.196491
	} else {
		approx_clone_scale_factor <- (-4.684e-05 * num_clones) + 3.171e-01
	}
	
	bound_num(approx_clone_scale_factor, lowerbound = 0.05, upperbound = 1)
}

# test ggobject for label addition

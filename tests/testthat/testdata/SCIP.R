data("mini_clonotype_data","mini_seurat_obj")

test_pbmc <- mini_seurat_obj
test_tcr <- mini_clonotype_data
test_integrated_pbmc <- readRDS("testdata/mini_integrated_seurat_obj.rds")

rm("mini_clonotype_data","mini_seurat_obj")

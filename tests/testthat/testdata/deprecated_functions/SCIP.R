data("mini_clonotype_data", "mini_seurat_obj")

test_pbmc <- mini_seurat_obj
test_tcr <- mini_clonotype_data
test_integrated_pbmc <- getdata(
    "deprecated_functions", "mini_integrated_seurat_obj"
)

rm("mini_clonotype_data", "mini_seurat_obj")

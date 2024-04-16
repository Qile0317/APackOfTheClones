test_pbmc <- getdata("v0", "mini_seurat_obj")
test_tcr <- getdata("v0", "mini_clonotype_data")

test_that("gg_color_hue() works", {
  expect_identical(gg_color_hue(1), "#F8766D")
  expect_identical(gg_color_hue(2), c("#F8766D", "#00BFC4"))
  expect_identical(gg_color_hue(13)[13], "#FF65AC")
})

test_that("gen_cluster_color_hashmap() works", {
  
  test_hashmap <- hash::hash()
  test_hashmap["cluster 0"] <- "#F8766D"
  test_hashmap["cluster 1"] <- "#00BA38"
  test_hashmap["cluster 2"] <- "#619CFF"
  
  expect_identical(gen_cluster_color_hashmap(3), test_hashmap)
})

test_that("insert_colors() works", {
  test_cluster_df <- df_full_join(test_cluster_lists)
  new_df <- insert_colors(test_cluster_df, 5)
  
  expect_identical(
    new_df %>% dplyr::select(color),
    test_cluster_df %>%
      dplyr::mutate(
        color = c(
          "#F8766D", "#F8766D", "#F8766D", "#F8766D",
          "#F8766D", "#F8766D", "#F8766D", "#F8766D",
          "#F8766D", "#A3A500", "#A3A500", "#A3A500", 
          "#A3A500", "#A3A500", "#A3A500", "#A3A500",
          "#A3A500", "#A3A500", "#00BF7D", "#00BF7D", 
          "#00BF7D", "#00BF7D", "#00BF7D", "#00BF7D", 
          "#00BF7D", "#00BF7D", "#00BF7D", "#00B0F6",
          "#00B0F6", "#00B0F6", "#00B0F6", "#00B0F6", 
          "#00B0F6", "#E76BF3", "#E76BF3", "#E76BF3", 
          "#E76BF3", "#E76BF3"
        )
      ) %>%
      dplyr::select(color)
  )
})

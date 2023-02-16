# test_that() is the most basic testing func
# describe() should verify that the right things were implemented

# test data, three cluster lists. c1 and c2 overlap while c3 is far away
c1 <- list(
  "x" = c(2.7, 5.3, 4, 4, 4, 5.977272, 2.022728, 3.271739, 5.249012),
  "y" = c(4.367544, 4.367544, 6.264911, 2.470178, 6.264911, 6.565567,
          6.565567, 8.127607, 8.428263),
  "rad" = c(1.3, 1.3, 1, 1, 1, 1, 1, 1, 1),
  "centroid" = c(4, 5),
  "clRad" = 2.977272
  )

c1_shifted_by_4_5 <- list(
  "x" = c(6.7, 9.3, 8, 8, 8, 9.977272, 6.022728,
          7.271739, 9.249012),
  "y" = c(9.367544, 9.367544, 11.264911, 7.470178,
          11.264911, 11.565567, 11.565567,
          13.127607, 13.428263),
  "rad" = c(1.3, 1.3, 1, 1, 1, 1, 1, 1, 1),
  "centroid" = c(4, 5),
  "clRad" = 2.977272
  )

c1_shifted_to_9_0 <- list(
  "x" = c(7.7, 10.3, 9, 9, 9, 10.977272, 7.022728,
          8.271739, 10.249012),
  "y" = c(-0.632456, -0.632456, 1.264911, -2.529822,
          1.264911, 1.565567, 1.565567, 3.127607,
          3.428263),
  "rad" = c(1.3, 1.3, 1, 1, 1, 1, 1, 1, 1),
  "centroid" = c(9,0),
  "clRad" = 2.977272
)

c2 <- list("x" = c(1.0196581, 4.9196581, 3.0606838, 3.0581197, 3.0581197,
                   0.1059401),
           "y" = c(3.0114986, 3.0114986, 5.9770027, 0.1663845, 5.8566128,
                   6.3901255),
           "rad" = c(2, 1.9, 1.6, 1.5, 1.5, 1.5),
           "centroid" = c(3, 4),
           "clRad" = 3.819658)

c3 <- list("x" = c(19.2381, 20.6381, 20.12381, 20.10952, 20.10952),
           "y" = c(19.6386, 19.6386, 20.72281, 18.67391, 20.60328),
           "rad" = c(0.8, 0.6, 0.6, 0.5, 0.5),
           "centroid" = c(20, 20),
           "clRad" = 1.238095)

# aggregate them
test_cluster_lists <- list(c1, c1_shifted_by_4_5, c1_shifted_to_9_0,c2,c3)

test_centroids <- list(
  "c1" = c(4, 5),
  "c1_shifted_by_4_5" = c(4, 5),
  "c1_shifted_to_9_0" = c(9, 0),
  "c2" = c(3, 4),
  "c3" = c(20, 20)
)

test_radii <- list(
  "c1" = c(), "c1_shifted_by_4_5" = c(),
  "c1_shifted_to_9_0" = c(), "c2" = c(3, 4),
  "c3" = c()
)

for (i in 1:length(test_cluster_lists)) {
  test_radii[[i]] <- test_cluster_lists[[i]][[3]]
}

# plot_clusters(df_full_join(list(c1,c1_shifted_by_4_5,
#                                 c1_shifted_to_9_0, c2,c3)))

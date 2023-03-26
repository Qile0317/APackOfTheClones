# it creates the OPPOSITE vector so there is no need for G to be negative
distV <- function(c1,c2){
  return(c(c1$centroid[1]-c2$centroid[1],
           c1$centroid[2]-c2$centroid[2]))
} # works

#polar form conversion from component form. its with respect to x axis. - works
polV <- function(vec) {
  return(c("magnitude" = sqrt(sum(vec^2)),
           "direction" = atan2(vec[2], vec[1])))
  }

#polar distance vector - works
pdV <- function(c1 ,c2) {
  return(polV(distV(c1, c2)))
}

#converts polar form to component form of vector - works
comV <- function(Pvec) {
  return(unname(c(Pvec[1] * cos(Pvec[2]), Pvec[1] * sin(Pvec[2]))))
}

#sum vectors in a list, very specific for this function
sumL <- function(vec_list) {
  ans <- c(0, 0)
  for (element in vec_list) {
    ans <- ans + element
  }
  return(ans)
}

#function to check if 2 cluster lists overlap, with a threshold.
do_cl_intersect <- function(Cn, Cm, thr = 1) {
  if (any(is.na(Cn)) || any(is.na(Cm))) {return(FALSE)}

  #calculate euclidean distance of centroids
  centroid_xdif <- (Cn$centroid[1] - Cm$centroid[1])
  centroid_ydif <- (Cn$centroid[2] - Cm$centroid[2])
  centroid_euc_dist <- sqrt((centroid_xdif^2) + (centroid_ydif^2))

  return((centroid_euc_dist + thr) < (Cn$clRad + Cm$clRad))
}

# helper function for iterative repulsion to check if the current input is valid
do_proceed <- function(inp, i, j, thr) {
  if (i != j) {
    if (!any(is.na(inp[[i]]))) {
      if (!any(is.na(inp[[j]]))) {
        if (do_cl_intersect(inp[[i]], inp[[j]], thr)) {
          return(TRUE)
        }
      }
    }
  }
  return(FALSE)
}

# iterative repulsion. inp is a list of clusterlists
repulse_cluster <- function(
    inp,
    thr = 0,
    G = 6e-11,
    max_iter = 100,
    dist_adjust = 2) {

  #init variables
  inp_len <- length(inp)
  iter_count <- 0

  transvec <- list()
  for (i in 1:inp_len) {transvec[[i]]<-c(0, 0)}
  curr_repulsion_vec <- transvec
  blank_vec <- transvec

  proceed <- FALSE
  change <- TRUE

  while(iter_count <= max_iter){
    if (!change) {return(inp)}
    iter_count <- iter_count + 1

    for(i in 1:inp_len) {
      currChange <- FALSE

      for(j in 1:inp_len) {

          proceed <- do_proceed(inp, i, j, thr)

          if (proceed) {
            polar_dist_vec <- pdV(inp[[i]], inp[[j]]) # find polar distance vector between centroids

            #Find polar repulsion vec
            polar_repulsion_vec <- c(
              (G * inp[[i]][[5]] * inp[[j]][[5]]) /
                (unname(polar_dist_vec["magnitude"]) + dist_adjust)^2,
              polar_dist_vec["direction"]
            )

            #store as component form
            curr_repulsion_vec[[j]] <- comV(polar_repulsion_vec)
            currChange <- TRUE

          }

          if (!currChange) {change <- FALSE}else {change <- TRUE}
      }

      # inner j looping complete. summing transforming vectors
      if (proceed) {
        transvec[[i]] <- sumL(curr_repulsion_vec)

        for(i in 1:inp_len){
          inp[[i]] <- trans_coord(inp[[i]], transvec[[i]])
        }

        curr_repulsion_vec <- blank_vec #reset repulsion
      }
    }
  }
  return(inp)
}

# TODO: could make it so that if a cluster already isn't touching anything, that repulsion wont apply to it (no vectors added)

#also if i ever wanted to optimize G, i can make a cost function for how close new centroids are to old centroids and account for how many iterations there were

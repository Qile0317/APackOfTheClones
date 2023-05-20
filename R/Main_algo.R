# This script contain functions to pack a single cluster from a vector of 
# numbers that represent radii of circles. The radii vector would be generated 
# from the integrated sc-RNAseq and TCR library.

# Node constructor for the circularly doubly linked list
node <- R6::R6Class(
  "Node",
  public = list(
    val = NULL,
    nxt = NULL,
    prv = NULL,
    
    initialize = function(
      val = NULL,
      nxt = NULL,
      prv = NULL
    ) {
      self$val <- val
      self$nxt <- nxt
      self$prv <- prv
    },
    
    print = function() {
      cat("x:", self$val$x, "y:", self$val$y, "rad:", self$val$rad, "\n")
    }
  )
)

# link the Nodes circularly and to each other in a given list of nodes "a"
init_boundary <- function(a) {
  for (i in 1:(length(a) - 1)) {
    a[[i + 1]]$prv <- a[[i]]
    a[[i]]$nxt <- a[[i + 1]]
  }

  a[[length(a)]]$nxt <- a[[1]]
  a[[1]]$prv <- a[[length(a)]]
}

#function to find "distance" between 2 elements of linked list
fwd_dist <- function(c, d) {
  count <- 0
  circ <- c
  while (!identical(circ, d)) {
    count <- count + 1
    circ <- circ$nxt
  }
  count
}

#initializing/inserting circle
insert_circle <- function(c, d, e) {
  if ((!identical(c$nxt, d)) || (!identical(d$prv, c))) {
    stop("Two circles not adjacent")
  }else {
    c$nxt <- e
    e$prv <- c
    d$prv <- e
    e$nxt <- d
  }
}

#removes segment between c d as one moves forwards
fwd_remove <- function(c, d) {

  if (identical(c, d)) {stop("Circles are the same.")}
  if (identical(c$nxt, d)) {stop("Circles are consecutive.")}

  circ <- c$nxt
  while (!identical(circ, d)) {
    circ$prv$nxt <- circ$nxt
    circ$nxt$prv <- circ$prv
    circ <- circ$nxt
  }
}

# Functions related to the geometry of circles ###############################

#distance from the centre of a circle to the origin.
centre_dist <- function(c) {
  return(sqrt((c$val[[2]])^2 + (c$val[[3]])^2))
}

#fit tangent circle function
fit_tang_circle <- function(C1, C2, C3) {

  x1 <- C1$val$x
  x2 <- C2$val$x
  
  y1 <- C1$val$y
  y2 <- C2$val$y
  
  r1 <- C1$val$rad
  r2 <- C2$val$rad
  r <- C3$val$rad

  distance <- sqrt((x1 - x2)^2 + (y1 - y2)^2)

  if (distance > (r1 + r2 + (2 * r))) {
    stop("Gap too large.")
  }

  cos_sig <- (x2 - x1)/distance
  sin_sig <- (y2 - y1)/distance
  cos_gam <- (distance^2 + (r + r1)^2 - (r + r2)^2)/(2*distance*(r + r1))
  sin_gam <- sqrt(1 - (cos_gam)^2)

  C3$val$x <- x1 + (r + r1) * (cos_sig * cos_gam - sin_sig * sin_gam)
  C3$val$y <- y1 + (r + r1) * (cos_sig * sin_gam + sin_sig * cos_gam)
  return(C3)
}

#Note: fit_tang_circle! fits C3 such that C1,C2,C3 are arranged counterclockwise
# place three circles in the center
place_starting_three <- function(C1, C2, C3) {
  C1$val$x <- -1 * (C1$val$rad)
  C2$val$x <- C2$val$rad

  fit_tang_circle(C2, C1, C3) #BM's original note: it seems like it might be necessary to initialise with opposite orientation

  #calculate the centroid of their centers, and translate each circle by it
  centroid_x <- (C1$val[[2]] + C2$val[[2]] + C3$val[[2]])/3
  centroid_y <- (C1$val[[3]] + C2$val[[3]] + C3$val[[3]])/3

  C1$val[[2]] <- C1$val[[2]] - centroid_x
  C2$val[[2]] <- C2$val[[2]] - centroid_x
  C3$val[[2]] <- C3$val[[2]] - centroid_x

  C1$val[[3]] <- C1$val[[3]] - centroid_y
  C2$val[[3]] <- C2$val[[3]] - centroid_y
  C3$val[[3]] <- C3$val[[3]] - centroid_y
}

################## circle packing algos ##############################

#finds the closest circle to the origin in the linked list containing c
closest <- function(c){
  closest <- c
  circ <- c$nxt
  while (!identical(circ, c)) {
    if (centre_dist(closest) > centre_dist(circ)) {
      closest <- circ
    }
    circ <- circ$nxt
  }
  closest
}

# alias function for closest_place that DOESNT mutate
tang_circle_dist <- function(c,d,e){
  centre_dist(fit_tang_circle(c$clone(),d$clone(),e$clone()))
}

#' place circle for try_place = TRUE - bugged maybe
#'
#' Locating the pair of successive circles, c, c$nxt, with the following
#' property: amongst all pairs of successive circles on the boundary, this pair
#' minimizes distance from the center of d to the origin, when d is fitted
#' tangent to this pair.
#' 
#' @noRd
#' 
closest_place <- function(c, d){
  closest <- c
  circ <- c$nxt
  while (!identical(circ,c)) {
    if (tang_circle_dist(closest, closest$nxt, d) > tang_circle_dist(circ, circ$nxt, d)) {
      closest <- circ
    }
    circ <- circ$nxt
  }
  return(closest)
}

#Checking for overlaps between a fitted circle and others on the boundary.
do_intersect <- function(c, d) {
  xdif <- (c$val$x - d$val$x)^2
  ydif <- (c$val$y - d$val$y)^2
  euc_dist <- sqrt(xdif + ydif)
  return(euc_dist < ((c$val$rad) + (d$val$rad)))
} #could reformulate as discrepancy > threshold

#convenience function to tidy up overlap_check, called quite alot 
geod_dist <- function(Cm, Cn, C) {
  min(fwd_dist(Cn, C), fwd_dist(C, Cm))
}

#overlap check of 3 circles. From profiling, this seems to be the function called the most and longest
overlap_check <- function(Cm, Cn, C) {
  C_em <- Cm
  C_en <- Cn
  obstruct <- list()

  #collect circles that C intersects(if any)by adding intersectors to 'obstruct'
  circ <- Cn$nxt
  while (!identical(circ, Cm)) {
    if (do_intersect(circ, C)) {
      obstruct <- c(obstruct, circ)
    }
    circ <- circ$nxt
  }

  LenObs <- length(obstruct)  #if there are any intersections

  if (LenObs > 0) { #find the one closest to {Cm, Cn}, (distance is in number of steps)
    nearest <- obstruct[[1]]
    for(i in 1:LenObs){
      if(geod_dist(Cm,Cn,obstruct[[i]]) < geod_dist(Cm,Cn,nearest)){
        nearest <- obstruct[[i]]
      }
    }
    if (fwd_dist(Cn, nearest) <= fwd_dist(nearest, Cm)) {  #if the distance is realised fwd, change C_en
      C_en <- nearest
    }else { #if distance is realised bkwd and not fwd, change C_em
      C_em <- nearest
    }
  }

  if ((identical(C_em, Cm)) && (identical(C_en, Cn))) {
    return("clear")
  }
  return(c(C_em, C_en))
}

# handling of degenerate cases of 1 or 2 circles for circle layout
is_degenerate_case <- function(lenCirc) {
  lenCirc == 1 || lenCirc == 2
} 

handle_degenerate_cases <- function(
  lenCirc, circles, centroid, rad_scale, verbose
) {
  if (verbose) {progress_bar(1,1)}
  if (lenCirc == 1) {
    return(list(
      "x" = circles[[1]]$val[[2]] + centroid[1],
      "y" = circles[[1]]$val[[3]] + centroid[2],
      "rad" = circles[[1]]$val[[6]] * rad_scale,
      "centroid" = centroid,
      "clRad" = circles[[1]]$val[[6]] * rad_scale
    ))
  }
  
  # if lenCirc == 2
  # transform the x coordinates to place left and right of center
  circles[[1]]$val[[2]] <- -1 * (circles[[1]]$val[[6]])
  circles[[2]]$val[[2]] <- circles[[2]]$val[[6]]
  
  return(list(
    "x" = c(circles[[1]]$val[[2]], circles[[2]]$val[[2]]) + centroid[1],
    "y" = rep(centroid[2], 2),
    "rad" = c(circles[[1]]$val[[6]], circles[[2]]$val[[6]]) * rad_scale,
    "centroid" = centroid,
    "clRad" = 0.5 * rad_scale * (
      circles[[1]]$val[[6]] + circles[[2]]$val[[6]]
    )
  ))
}

# removed original circle_layout

# vectorized circle_layout - outputs list
pack_into_clusterlists <- function(
  sizes, centroids, num_clusters, rad_scale = 1,
  ORDER = TRUE, try_place = FALSE, verbose = TRUE
){
  output_list <- list()
  for(i in 1:num_clusters){
    if (is.null(sizes[[i]])) {
      output_list[[i]] <- list()
    }else{
      if(verbose){
        message(paste("\npacking cluster", as.character(i-1)))
      }
      output_list[[i]] <- circle_layout(
        sizes[[i]],
        centroid = centroids[[i]],
        rad_scale_factor = rad_scale,
        ORDER = ORDER,
        try_place = try_place,
        progbar = verbose
      )
    }
  }
  output_list
}

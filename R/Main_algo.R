#This script is functions to mathematically form a SINGLE cluster from an array
#of numbers that represent radii of circles, Of course, this radii array is
#generated from the integrated sc-RNAseq and TCR library.

#function to initialize list into the circular doubly linked list.
init_boundary <- function(a){

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
  while (!identical(circ,d)) {
    count <- count + 1
    circ <- circ$nxt
  }
  count
}

#initializing/inserting 3 circles
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
  sqrt((c$val[[2]])^2 + (c$val[[3]])^2)
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

  if (distance > (r1 + r2 + 2 * r)) {
    stop("Gap too large.")
  }else {distance} # ?

  cos_sig <- (x2 - x1)/distance
  sin_sig <- (y2 - y1)/distance
  cos_gam <- (distance^2 + (r + r1)^2 - (r + r2)^2)/(2*distance*(r + r1))
  sin_gam <- sqrt(1 - (cos_gam)^2)

  C3$val[[2]] <- x1 + (r + r1) * (cos_sig * cos_gam - sin_sig * sin_gam)
  C3$val[[3]] <- y1 + (r + r1) * (cos_sig * sin_gam + sin_sig * cos_gam)
  C3
}

# place three circles in the center
place_starting_three <- function(C1, C2, C3) {
  C1$val[[2]] <- -1 * (C1$val[[6]])
  C2$val[[2]] <- C2$val[[6]]

  fit_tang_circle(C1, C2, C3) #BM's original note: it seems like it might be necessary to initialise with opposite orientation

  centroid_x <- (C1$val[[2]] + C2$val[[2]] + C3$val[[2]])/3
  centroid_y <- (C1$val[[3]] + C2$val[[3]] + C3$val[[3]])/3

  C1$val[[2]] <- C1$val[[2]] - centroid_x
  C2$val[[2]] <- C2$val[[2]] - centroid_x
  C3$val[[2]] <- C3$val[[2]] - centroid_x

  C1$val[[3]] <- C1$val[[3]] - centroid_y
  C2$val[[3]] <- C2$val[[3]] - centroid_y
  C3$val[[3]] <- C3$val[[3]] - centroid_y
}

################## circle packing algos #############################################

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
  return(closest)
}

#Locating the pair of successive circles, c, c.s, with the following property: amongst all pairs of
#successive circles on the boundary, this pair minimizes distance from the center of d to the origin, when d is
#fitted tangent to this pair.

closest_place <- function(c, d){
  closest <- c
  circ <- c$nxt
  while (!identical(circ,c)) {
    if (centre_dist(fit_tang_circle(closest, closest$nxt, d)) > centre_dist(fit_tang_circle(circ, circ$nxt, d))) {
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

#convenience function to tidy up overlap_check
geod_dist <- function(Cm, Cn, C) {
  min(fwd_dist(Cn, C), fwd_dist(C, Cm))
}

#overlap check of 3 circles.
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

  LenObs <- length(obstruct)  #if there are any intersectiosn

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
  }else {
    return(c(C_em, C_en))
  }
}

#Cluster radius estimator that assumes imput is a list of x,y,r. Runs in linear time by 1 scan as list isnt sorted
#the x or y MUST BE the first in the list, and assumes centeroid at 0,0
est_rad <- function(coords){
  cc <- 0
  max_x <- 0
  for(i in 1:length(coords$x)){ # x should be the first in the coord list
    current_coord <- coords$x[i]
    if(current_coord > max_x){
      max_x <- current_coord
      cc <- i
    }
  }
  max_radius <- coords$rad[cc] #radii should be third in the list
  centroid_x <- coords$centroid[1] #centroid should be fourth in the list
  return(max_x + max_radius - centroid_x)
}

#The circle layout function.###################################
# could export this just as a pure function
#It takes an input vector of radii, and returns a vector of centre coordinates of the corresponding circles in the layout.

#Optional arguments are:
# rad decrease: after packing, if the circles should be slightly smaller to have a small gap of the input value between them
#"order": default = true
#if true it sorts the input vector in descending order before packing
#"try_place":  default is true
#if true the algorithm tries to place each new circle to be added to the packing as close to the origin as possilble,
#if false the algorithm tries to place each new circle to be added to the packing tangent to the closest circle on the boundary.

#   IMPORTANT:
#   this function does not incorporate colors! functionality will be added later. its very simple to do in ggplot

circle_layout <- function(input_rad_vec, centroid = c(0, 0),
                          rad_decrease = 1,
                          ORDER = TRUE, try_place = TRUE,
                          progbar = TRUE, print_BL = FALSE) {

  if (identical(input_rad_vec, list())) {return(NULL)}

  if (ORDER) {input_rad_vec <- sort(input_rad_vec, decreasing = TRUE)}

  # Initialise the circles with radii (not areas) as specified in input_rad_vec, and no boundary relations.
  circles <- list() #not sure if list of vector is better/faster here
  for (i in 1:length(input_rad_vec)) {
    currN <- paste("Circle", as.character(i), sep = "_")
    currCirc <- node$new(val = list(area = NULL, x = 0, y = 0,
                         color = NULL, label = currN,
                         rad = input_rad_vec[i]))
    circles <- append(circles, currCirc)
  }

  lenCirc <- length(circles)

  #Taking care of "degenerate" cases when there are only one or two circles
  if (lenCirc == 1) {
    return(list("x" = circles[[1]]$val[[2]] + centroid[1],
                "y" = circles[[1]]$val[[3]] + centroid[2],
                "rad" = circles[[1]]$val[[6]] * rad_decrease,
                "centroid" = centroid,
                "clRad" = circles[[1]]$val[[6]]))
  }

  if (lenCirc == 2) {
    # transform the x coordinates to place left and right of center
    circles[[1]]$val[[2]] <- -1 * (circles[[1]]$val[[6]])
    circles[[2]]$val[[2]] <- circles[[2]]$val[[6]]

    return(list("x" = c(circles[[1]]$val[[2]] + centroid[1],
                        circles[[2]]$val[[2]] + centroid[1]),
                "y" = c(centroid[2], centroid[2]),
                "rad" = c(circles[[1]]$val[[6]],
                          circles[[2]]$val[[6]]) * rad_decrease,
                "centroid" = centroid,
                "clRad" = 0.5 * (circles[[1]]$val[[6]] +
                                  circles[[2]]$val[[6]])))
  }

  # Place the first three circles to be mutually tangent, with centroid at the origin.
  place_starting_three(circles[[1]], circles[[2]], circles[[3]])

  # Initialise the boundary
  init_boundary(list(circles[[1]],circles[[2]],circles[[3]]))

  #Loop through the remaining circles,fitting them
  j <- 4
  while (j <= lenCirc) {
    if (try_place) {
      cl <- closest_place(circles[[j-1]], circles[[j]])
    }else {
      cl <- closest(circles[[j-1]])
    }

    fit_tang_circle(cl, cl$nxt, circles[[j]])

    # Check for overlaps and update, refit and recheck until "clear"
    check <- overlap_check(cl, cl$nxt, circles[[j]])

    if(identical(check, "clear")) {
      insert_circle(cl, cl$nxt, circles[[j]])
      j <- j + 1
      if (progbar) {progress_bar(j,lenCirc)}

    }else{
      while(!identical(check,"clear")){
        Cm <- check[[1]]
        Cn <- check[[2]]

        fwd_remove(Cm,Cn)

        fit_tang_circle(Cm, Cn, circles[[j]])

        check <- overlap_check(Cm, Cn, circles[[j]])

        if (identical(check, "clear")){
          insert_circle(Cm, Cn, circles[[j]])
          j <- j + 1
        }
      }
    }
    if (print_BL) {print(clength(circles[[j-1]]))} #for debugging
  }

  ans <- list() #in the future i can put the colors in the prior functions.
  cc <- 1
  Rvec <- c()
  Xvec <- c()
  Yvec <- c()

  for (c in circles) { #better practisce would be to initialize vectors of zeros first
      Rvec <- c(Rvec, c$val[[6]])
      Xvec <- c(Xvec, c$val[[2]])
      Yvec <- c(Yvec, c$val[[3]])
  }

  # construct output cluster list
  ans <- list("x" = Xvec,
              "y" = Yvec,
              "rad" = Rvec,
              "centroid" = centroid,
              "clRad" = 0)

  # transform the x and y coordinates to the centroid if its not 0,0.
  if(!identical(centroid, c(0, 0))){ans <- trans_coord(ans)}

  # estimate radius of cluster for repulsion
  ans[[5]] <- est_rad(ans)

  # scale radius
  if (rad_decrease != 1) {ans[[3]] <- Rvec * rad_decrease}

  # return
  message("")
  return(ans)
}

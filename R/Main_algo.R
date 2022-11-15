#This script is functions to mathematically form a SINGLE cluster from an array
#of numbers that represent radii of circles, Of course, this radii array is
#generated from the integrated sc-RNAseq and TCR library.

#ill get rid of these soon.
library(cli)
library(tidyverse)
library(ggdag)
library(R6)
library(pryr)
library(ggplot2)
library(ggraph)
library(ggforce)

#function to convert list to the circular doubly linked list.
init_boundary <- function(a){
  for(i in 1:(length(a)-1)){
    a[[i+1]]$prv <- a[[i]]
    a[[i]]$nxt <- a[[i+1]]}
  a[[length(a)]]$nxt <- a[[1]]
  a[[1]]$prv = a[[length(a)]]
}

#function to find "distance" between 2 elements of linked list
fwd_dist <- function(c,d){
  count <- 0
  circ <- c
  while(identical(circ,d)==FALSE){
    count <- count + 1
    circ <- circ$nxt
  }
  count
}

#initializing/inserting 3 circles
insert_circle <- function(c,d,e){
  if((identical(c$nxt, d) == FALSE) | (identical(d$prv, c) == FALSE)){
    stop("Two circles not adjacent")
  }else{
    c$nxt <- e
    e$prv <- c
    d$prv <- e
    e$nxt <- d
  }
}

#removes segment between c d as one moves forwards
fwd_remove <- function(c,d){
  if(identical(c,d)==TRUE){
    stop("Circles are the same.")
  }else if(identical(c$nxt,d)){
    stop("Circles are consecutive.")
  }else{
    circ <- c$nxt
    #removed <- c()
    while(identical(circ,d)==FALSE){
      circ$prv$nxt <- circ$nxt
      circ$nxt$prv <- circ$prv
      circ <- circ$nxt
    }
  }
}

# Functions related to the geometry of circles ###############################

#distance from the centre of a circle to the origin.
centre_dist <- function(c){
  sqrt((c$val[[2]])^2 + (c$val[[3]])^2)
}

#fit tangent circle function
fit_tang_circle <- function(C1,C2,C3){
  x1 = C1$val$x
  x2 = C2$val$x
  y1 = C1$val$y
  y2 = C2$val$y
  r1 = C1$val$rad
  r2 = C2$val$rad
  r = C3$val$rad

  distance = sqrt((x1 - x2)^2 + (y1 - y2)^2)

  if (distance > (r1 + r2 + 2*r)){
    stop("Gap too large.")}
  else{distance}

  cos_sig = (x2 - x1)/distance
  sin_sig = (y2 - y1)/distance
  cos_gam = (distance^2 + (r + r1)^2 - (r + r2)^2)/(2*distance*(r + r1))
  sin_gam = sqrt(1 - (cos_gam)^2)

  C3$val[[2]] = x1 + (r + r1)*(cos_sig*cos_gam - sin_sig*sin_gam)
  C3$val[[3]] = y1 + (r + r1)*(cos_sig*sin_gam + sin_sig*cos_gam)
  C3
}

place_starting_three <- function(C1,C2,C3){
  C1$val[[2]] = -1*(C1$val[[6]])
  C2$val[[2]] = C2$val[[6]]

  fit_tang_circle(C1,C2,C3) #BM's original note: it seems like it might be necessary to initialise with opposite orientation

  centroid_x <- (C1$val[[2]]+C2$val[[2]]+C3$val[[2]])/3
  centroid_y <- (C1$val[[3]]+C2$val[[3]]+C3$val[[3]])/3

  C1$val[[2]] <- C1$val[[2]] - centroid_x
  C2$val[[2]] <- C2$val[[2]] - centroid_x
  C3$val[[2]] <- C3$val[[2]] - centroid_x

  C1$val[[3]] <- C1$val[[3]] - centroid_y
  C2$val[[3]] <- C2$val[[3]] - centroid_y
  C3$val[[3]] <- C3$val[[3]] - centroid_y
}

################## circle packing algos #############################################

#finds the closest circle to the origin in the linked list containing c
closest <- function(c, showProcess=FALSE){
  closest <- c
  circ <- c$nxt
  if(showProcess){print("--------------")}
  while(identical(circ,c)==FALSE){

    if(showProcess){
      print(closest$val$label)
      print(centre_dist(closest))
      print(circ$val$label)
      print(centre_dist(circ))
      print("--------------")}

    if(centre_dist(closest) > centre_dist(circ)){
      closest <- circ}
    circ <- circ$nxt
  }
  return(closest)
}

#Locating the pair of successive circles, c, c.s, with the following property: amongst all pairs of
#successive circles on the boundary, this pair minimizes distance from the center of d to the origin, when d is
#fitted tangent to this pair.


closest_place <- function(c,d){
  closest = c
  circ = c$nxt
  while(identical(circ,c)==FALSE){
    if(centre_dist(fit_tang_circle(closest, closest$nxt, d)) > centre_dist(fit_tang_circle(circ, circ$nxt, d))){
      closest <- circ
    }
    circ <- circ$nxt
  }
  return(closest)
}

#Checking for overlaps between a fitted circle and others on the boundary.
do_intersect <- function(c,d){
  xdif <- (c$val$x-d$val$x)^2
  ydif <- (c$val$y-d$val$y)^2
  dista <- sqrt(xdif+ydif)
  return(dista < ((c$val$rad)+(d$val$rad))) #could reformulate as discrepancy > threshold
}

#convenience function to tidy up overlap_check
geod_dist <- function(Cm,Cn,C){
  min(fwd_dist(Cn,C), fwd_dist(C,Cm))
}

#overlap check of 3 circles.
overlap_check <- function(Cm,Cn,C){
  C_em <- Cm
  C_en <- Cn
  obstruct <- list()

  #collect circles that C intersects, if any by adding intersectors to an 'obstruct' list.
  circ <- Cn$nxt
  while(identical(circ,Cm)==FALSE){
    if(do_intersect(circ,C)){
      obstruct <- c(obstruct,circ) #root of my problems
    }
    circ <- circ$nxt
  }

  LenObs <- length(obstruct)  #if there are any intersectiosn

  if(LenObs > 0){             #find the one closest to {Cm, Cn}, where distance is in number of steps
    nearest <- obstruct[[1]]
    for(i in 1:LenObs){
      if(geod_dist(Cm,Cn,obstruct[[i]]) < geod_dist(Cm,Cn,nearest)){
        nearest <- obstruct[[i]]
      }
    }
    if(fwd_dist(Cn,nearest) <= fwd_dist(nearest,Cm)){     #if the distance is realised fwd, change C_en
      C_en <- nearest
    }else{                                                #if distance is realised bkwd and not fwd, change C_em
      C_em <- nearest
    }
  }

  if((identical(C_em,Cm))&(identical(C_en,Cn))){
    return("clear")
  }else{
    return(c(C_em, C_en))
  }
}

#The circle layout function.###################################
#It takes an input vector of radii, and returns a vector of centre coordinates of the corresponding circles in the layout.
#Optional arguments are:
#"order": default = true
#if true it sorts the input vector in descending order before packing
#"try_place":  default is true
#if true the algorithm tries to place each new circle to be added to the packing as close to the origin as possilble,
#if false the algorithm tries to place each new circle to be added to the packing tangent to the closest circle on the boundary.

#   IMPORTANT:
#   this function does not incorporate colors! functionality will be added later. its very simple to do in ggplot

circle_layout <- function(input_rad_vec, centroid = c(0,0), ORDER=TRUE, try_place=TRUE, progbar=TRUE, print_BL=FALSE){
  if(ORDER){input_rad_vec <- rev(sort(input_rad_vec))}

  # Initialise the circles with radii (not areas) as specified in input_rad_vec, and no boundary relations.
  circles = list() #not sure if list of vector is better/faster here
  for(i in 1:length(input_rad_vec)){
    currN <- paste("Circle",as.character(i),sep="_")
    currCirc <- node$new(val=list(area=NULL,x=0,y=0,color=NULL,label=currN,rad=input_rad_vec[i]))
    circles <- append(circles,currCirc)
  }
  lenCirc = length(circles)

  #Taking care of "degenerate" cases when there are one or two circles (btw im pretty sure it doesnt work lol)
  if(lenCirc==1){
    return(c(circles[[1]]$val[[2]],circles[[1]]$val[[3]],circles[[1]]$val[[6]]))
  }else if(lenCirc==2){
    circles[[1]]$val[[2]] = -1*(circles[[1]]$val[[6]])
    circles[[2]]$val[[2]] = circles[[2]]$val[[6]]
    return(c(circles[[1]]$val,circles[[2]]$val))
  }

  # Place the first three circles to be mutually tangent, with centroid the origin.
  place_starting_three(circles[[1]], circles[[2]], circles[[3]])

  # Initialise the boundary
  init_boundary(list(circles[[1]],circles[[2]],circles[[3]]))

  #if(progbar){
  #for(i in 1:length(circles)){print(paste(circles[[i]]$val$label,":",circles[[i]]$val$rad,sep=""))}
  #}

  #Loop through the remaining circles,fitting them
  j <- 4
  while(j<=lenCirc){
    if(try_place==TRUE){
      cl <- closest_place(circles[[j-1]],circles[[j]])
    }else {
      cl <- closest(circles[[j-1]])
    }

    fit_tang_circle(cl,cl$nxt,circles[[j]])

    # Check for overlaps and update, refit and recheck until "clear"
    check <- overlap_check(cl, cl$nxt, circles[[j]])
    if(identical(check,"clear")){
      insert_circle(cl,cl$nxt,circles[[j]])
      j <- j + 1
      if(progbar){progress_bar(j,lenCirc)}
    }else{
      while(!identical(check,"clear")){
        Cm <- check[[1]]
        Cn <- check[[2]]

        fwd_remove(Cm,Cn)

        fit_tang_circle(Cm,Cn,circles[[j]]) #not sure abt what happens if a thing is returned? havent checked anything yet.

        check <- overlap_check(Cm,Cn,circles[[j]])

        if(identical(check,"clear")){
          insert_circle(Cm,Cn,circles[[j]])
          j <- j + 1
        }
      }
    }
    if(print_BL==TRUE){print(clength(circles[[j-1]]))} #for debugging
  }

  ans <- list() #in the future i can put the colors in the prior functions.
  cc <- 1
  Rvec <- c()
  Xvec <- c()
  Yvec <- c()

  for(c in circles){
      Rvec <- c(Rvec,c$val[[6]])
      Xvec <- c(Xvec,c$val[[2]])
      Yvec <- c(Yvec,c$val[[3]])
  }
  ans <- list(x=Xvec,y=Yvec,rad=Rvec,centroid=c(0,0),clRad=0)
  if(!identical(centroid,c(0,0))){ans <- trans_coord(ans,centroid)} #didnt test this lol
  ans[["clRad"]] <- est_rad(ans) #estimated radius of cluster, function found in utils.r
  return(ans)
}

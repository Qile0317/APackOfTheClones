# Repulsion = -G(m_1 + m_2)/d^2
#with iteration, this algo is O(mn^2) where m is the iteration threshold

#prereq function to find the repulsion force vector of 2 circles,
#the imputs are lists for 1 cluster with $x, $y, $rad, $centroid, $clstr_rad
#$centroid must be c(x,y)
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
sumL <- function(li) {
  ans <- c(0, 0)
  for (el in li) {
    ans[1] <- ans[1] + el[1]
    ans[2] <- ans[2] + el[2]
  }
  return(ans)
}

#function to check if 2 cluster lists overlap, with threshold.
do_cl_intersect <- function(Cn, Cm, thr = 1) {
  if (is.null(Cn) || is.null(Cm) || length(Cn) <= 0 || length(Cm) <= 0) {return(FALSE)}

  #calculate euclidean distance of centroids
  centroid_xdif <- (Cn$centroid[1] - Cm$centroid[1])
  centroid_ydif <- (Cn$centroid[2] - Cm$centroid[2])
  centroid_euc_dist <- sqrt((centroid_xdif^2) + (centroid_ydif^2))

  return((centroid_euc_dist - thr) < (Cn$clRad + Cm$clRad)) #something is wrong with Cn and Cm
}

# Error in Cn[[5]] + Cm[[5]] : non-numeric argument to binary operator
# I dont think it the fault of the function its something wrong with the cluster list creation

#However, there arent that many clusters so it should be fine.
#ASSUMES NO 2 CENTROIDS ARE THE SAME!!!!!!!!
#thr is how much overlap is acceptable.
#G is the strength of the repulsion. should be a small number
#max_iter is there to prevent it from running forever as its not super optimized.

#inp is a list of cluster lists. thr is the allowed border that isnt implemented.

repulse_cluster <- function(inp, thr = 1, G = 0.05, max_iter = 100){ # so the inp is a list of ?
  #keep a variable of if ANY vector is added. if the prev and curr are unchanged, then return.
  li <- length(inp)
  transvec <- list()
  for(i in 1:li){transvec[[i]]<-c(0,0)}
  ctvec <- transvec
  otvec <- transvec
  change <- TRUE
  i <- 0

  while(i <= max_iter){
    if(!change){
      return(inp)

      }else{ #this else is ac quite unessecary nesting
      i <- i + 1

      for(i in 1:li){
        currChange <- FALSE

        for(j in 1:li){

          #check edge case
          if (!is.null(inp[[i]]) && !is.null(inp[[j]])) {
            proceed <- TRUE
          }else {proceed <- FALSE}

          #iterative repulsion
          if (proceed) {
            if (do_cl_intersect(inp[[i]],inp[[j]],thr)) {
              polDV <- pdV(inp[[i]],inp[[j]]) # find distance vector between centroids

              if(polDV[1] != 0){
                #Find polar repulsion vec (with newtons formula) + converting to and storing component form
                ctvec[[j]] <- comV(c(G*((pi*(inp[[i]][[5]]^2 + inp[[j]][[5]]^2))/(polDV[1])^2), polDV[2])) #tbh pi and squaring are completely unessecary and one could just adjust G
                currChange <- TRUE
              }
            }
          }
          if(!currChange){change = FALSE}else{change = TRUE} #unessecary?
        }

        if (proceed) {
          transvec[[i]] <- sumL(ctvec)
          ctvec <- otvec #reset

          for(i in 1:li){
            inp[[i]] <- trans_coord(inp[[i]],transvec[[i]]) #transforms cluster to new centroid
          }
        }
        }
      }
  }
  return(inp)
}

# TODO: could make it so that if a cluster already isn't touching anything, that repulsion wont apply to it (no vectors added)

#also if i ever wanted to optimize G, i can make a cost function for how close new centroids are to old centroids and account for how many iterations there were

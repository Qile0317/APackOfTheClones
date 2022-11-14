# Repulsion = -G(m_1 + m_2)/d^2
#the coding train did this and assumed all masses are equivalent, which saves iteration. HOwever its not the case here.
#with iteration, this algo is O(mn^2) where m is the number of iterationsz needed to be made.

#prereq function to find the repulsion force vector of 2 circles,
#the imputs are lists for 1 cluster with $x, $y, $rad, $centroid, $clstr_rad
#$centroid must be c(x,y)
distV <- function(c1,c2){
  return(c(c1$centroid[1]-c2$centroid[1],
           c1$centroid[2]-c2$centroid[2]))
} # works

#polar form conversion from component form. its with respect to x axis. - works
polV <- function(vec){
  return(c("magnitude"=sqrt(sum(vec^2)),
           "direction"=atan2(vec[2],vec[1])))
  }

#polar distance vector - works
pdV <- function(c1,c2){return(polV(distV(c1,c2)))}

#converts polar form to component form of vector - works
comV <- function(Pvec){
  return(unname(c(Pvec[1]*cos(Pvec[2]),Pvec[1]*sin(Pvec[2]))))
}

#sum vectors in a list, very specific for this function
sumL <- function(li){
  ans <- c(0,0)
  for(el in li){
    ans[1] <- ans[1] + el[1]
    ans[2] <- ans[2] + el[2]
  }
  return(ans)
}

#function to check if 2 cluster lists overlap, with threshold. (cluster overlap check (coc))
do_cl_intersect <- function(Cn,Cm,thr){
  xdif <- (Cn[[4]][1]+Cm[[4]][1])^2
  ydif <- (Cn[[4]][2]+Cm[[4]][2])^2
  distance <- sqrt(xdif+ydif)
  return(distance+thr < Cn[[5]] + Cm[[5]])
}

#at the end, its possible to find optimal G by seeing how many iterations each G took. but unessecary

#However, there arent that many clusters so it should be fine.
#ASSUMES NO 2 CENTROIDS ARE THE SAME!!!!!!!!
#thr is how much overlap is acceptable.
#G is the strength of the repulsion. should be a small number
#max_iter is there to prevent it from running forever as its not super optimized.

#inp is a list of cluster lists. thr is the allowed border that isnt implemented.
repulse_cluster <- function(inp, thr = 1, G = 0.05, max_iter = 100){ #not sure what G should be, need to test.
  #keep a variable of if ANY vector is added. if the prev and curr are unchanged, then return.
  li <- length(inp)
  transvec <- list()
  for(i in 1:li){transvec[[i]]<-c(0,0)}
  ctvec <- transvec
  otvec <- transvec
  change <- TRUE
  i <- 0
  while(i <= max_iter){
    if(change == FALSE){return(inp)
      }else{
      i <- i + 1
      for(i in 1:li){
        currChange <- FALSE
        for(j in 1:li){
          #find distance vector between centroids
          if(do_cl_intersect(inp[[i]],inp[[j]],thr)==TRUE){
            polDV <- pdV(inp[[i]],inp[[j]])
            if(polDV[1] != 0){
              #Find repulsion vec (newtons formula) + storing component
              ctvec[[j]] <- comV(c(G*((pi*( #getting rid of *-1 helped. Lol looks like the original direction was already pointing opposite.
                inp[[i]][[5]]^2 +
                  inp[[j]][[5]]^2))/(polDV[1])^2), #"magnitude" (5 is cluster radius clRad)
                polDV[2])) #"direction"
              currChange <- TRUE #update that there is change within this iteration
            }
          }
          if(currChange == FALSE){change = FALSE}else{change = TRUE} #unessecary?
        }
        transvec[[i]] <- sumL(ctvec)
        ctvec <- otvec #reset
        for(i in 1:li){
          inp[[i]] <- trans_coord(inp[[i]],transvec[[i]]) #transforms cluster to new centroid
        }
        }
      }
  return(inp)
  }
}
#IT WORKS!! :D

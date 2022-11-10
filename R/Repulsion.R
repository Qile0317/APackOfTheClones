# Repulsion = -G(m_1 + m_2)/d^2
#the coding train did this and assumed all masses are equivalent, which saves iteration. HOwever its not the case here. 
#with iteration, this algo is O(mn^2) where m is the number of iterationsz needed to be made.

#prereq function to find the repulsion force vector of 2 circles,
#the imputs are lists for 1 cluster with $x, $y, $rad, $centroid, $clstr_rad
#$centroid must be c(x,y) 
distV <- function(c1,c2){ #does the order matter? have to investigate
  return(c(c1$centroid[1]-c2$centroid[1],
           c1$centroid[2]-c2$centroid[2]))
} #Error: $ operator is invalid for atomic vectors

#polar form conversion from component form. its with respect to x axis.
polV <- function(vec){
  return(c("magnitude"=sqrt(sum(vec^2)), 
           "direction"=atan2(vec[2],vec[1])))
  }

#polar distance vector
pdV <- function(c1,c2){return(polV(distV(c1,c2)))}

#converts polar form to component form of vector
comV <- function(Pvec){
  return(unname(c(Pvec[1]*cos(Pvec[2]),Pvec[1]*sin(Pvec[2]))))
}

#at the end, its possible to find optimal G by seeing how many iterations each G took. but unessecary

#However, there arent that many clusters so it should be fine.
#ASSUMES NO 2 CENTROIDS ARE THE SAME!!!!!!!!
#this can even be visualized
#UNIFNISHED! need to deal with prevC, currC, thr, and test G.

#inp is a list of cluster lists. thr is the allowed border that isnt implemented. 
repulse_cluster <- function(inp, thr = 1, G = 0.5, max_iter = 10){ #not sure what G should be, need to test.
  if(max_iter > 100){max_iter <- 100}
  #keep a variable of if ANY vector is added. if the prev and curr are unchanged, then return. 
  li <- length(inp)
  transvec <- rep(c(0,0), li)
  ctvec <- transvec
  otvec <- transvec
  prevC <- FALSE
  currC <- TRUE
  i <- 0 
  while(i <= max_iter){
    i <- i + 1 
    if(prevC == FALSE && currC == FALSE){ #if no changes, then its done
      return(inp) #idk if return is nessecary
    }else{
      prevC <- currC #not sure abt this
      for(i in 1:li){
        for(j in 1:li){
          #find distance vector between centroids
          polDV <- pdV(inp[[i]],inp[[j]])
          if(polDV[1] != 0){ #if its not the same, 
            #Find repulsion vec (newtons formula) + storing component
            ctvec[j] <- comV(c(-1*G*((pi*(
              inp[[i]]$clstr_rad^2 +
                inp[[j]]$clstr_rad^2))/(polDV[1]+thr)^2), #"magnitude"
              polDV[2])) #"direction"
            #update that theres change 
            if(ctvec[j] != c(0,0)){
              currC <- TRUE
            }else{currC <- FALSE}
          }
        }
        transvec[i] <- sum(ctvec[i])
        ctvec <- otvec #reset
        for(i in 1:li){
          inp[[i]]$centroid <- inp[[i]]$centroid+transvec[i]
          inp[[i]] <- trans_coord(inp[[i]]) #unfortunately doesnt mutate
        } 
        prevC <- currC
      }
    }
  }
  return(inp) #idk if return is nessecary
  }

repulse_cluster(clusterlist)  
  
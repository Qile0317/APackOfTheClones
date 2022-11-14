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

#at the end, its possible to find optimal G by seeing how many iterations each G took. but unessecary

#However, there arent that many clusters so it should be fine.
#ASSUMES NO 2 CENTROIDS ARE THE SAME!!!!!!!!
#this can even be visualized
#UNIFNISHED! need to deal with prevC, currC, thr, and test G.

#inp is a list of cluster lists. thr is the allowed border that isnt implemented.
repulse_cluster <- function(inp, thr = 1, G = 0.5, max_iter = 10){ #not sure what G should be, need to test.
  if(max_iter > 30){
    print("Max iterations too high, reducing to 30...")
    max_iter <- 30}
  #keep a variable of if ANY vector is added. if the prev and curr are unchanged, then return.
  li <- length(inp)
  transvec <- rep(c(0,0), li)
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
          if(!identical(inp[[i]][[4]],inp[[j]][[4]])){
            polDV <- pdV(inp[[i]],inp[[j]])
            if(polDV[1] != 0){ #if its not the same,
              #Find repulsion vec (newtons formula) + storing component
              ctvec[j] <- comV(c(-1*G*((pi*(
                inp[[i]]$clstr_rad^2 +
                  inp[[j]]$clstr_rad^2))/(polDV[1]+thr)^2), #"magnitude"
                polDV[2])) #"direction"
              currChange <- TRUE #update that there is change within this iteration
            }
          }
          if(currChange == FALSE){change = FALSE}
        }
        transvec[i] <- sum(ctvec[i])
        ctvec <- otvec #reset
        for(i in 1:li){
          inp[[i]]$centroid <- inp[[i]]$centroid+transvec[i]
          inp[[i]] <- trans_coord(inp[[i]]) #unfortunately doesnt mutate
        }
        }
      }
  return(inp)
  }
}

repulse_cluster(clusterlist)
main <- df_full_join(clusterlist)
plot_clusters(main) #lmao didnt work, also had 12 warnings

# script of miscalleneous/utility functions

#simple progress bar.
#In the future it should be possible to include ETA but it usually doesnt take too long anyway
progress_bar <- function (x, max = 100) {
  percent <- x / max * 100
  cat(sprintf('\r[%-50s] %d%%',
              paste(rep('=', percent / 2), collapse = ''),
              floor(percent)))
  if (x == max)
    cat('\n')
}

#area from radii vector
areaFromRad <- function(c){
  ans <- c()
  for(i in c){
    area <- pi*(i^2)
    ans <- c(ans,area)
  }
  return(ans)
}

#crucial in circle_layout!
#Cluster radius estimator that assumes imput is a list of x,y,r. Runs in linear time by 1 scan as list isnt sorted
#the x or y MUST BE the first in the list, and assumes centeroid at 0,0
est_rad <- function(coords){
  cc <- 0
  max <- 0
  for(i in 1:length(coords[[1]])){ #$x
    if(coords[[1]][i] > max){ #$if term in x is bigger than max
      max <- coords[[1]][i]
      cc <- i
    }
  }
  maxr <- coords[[3]][cc] #the third should be radius
  return(max+maxr) #+coords[[4]][1] fourth should be centroid
}

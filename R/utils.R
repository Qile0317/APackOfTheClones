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

#full join a list of lists into a dataframe with generated labels. memory inefficient
df_full_join <- function(clstr_list){ #i need option to put custom labels
  df <- data.frame(label="cluster 1", x=clstr_list[[1]]$x,
                   y=clstr_list[[1]]$y,r=clstr_list[[1]]$rad)
  for(i in 2:length(clstr_list)){
    df <- full_join(df,data.frame(label=paste("cluster", as.character(i)),
                            x=clstr_list[[i]]$x, y=clstr_list[[i]]$y,
                            r=clstr_list[[i]]$rad))
  }
  return(df)
}


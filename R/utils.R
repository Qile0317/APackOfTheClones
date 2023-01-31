# script of miscalleneous/utility functions

#simple progress bar
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
    area <- pi * (i^2)
    ans <- c(ans, area)
  }
  return(ans)
}

#crucial in circle_layout!
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
  return(max_x + centroid_x + max_radius)
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

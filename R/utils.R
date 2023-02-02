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

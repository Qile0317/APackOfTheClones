#script of the circular, doubly linked list implementation and also some utility functions for debugging

library(R6)

# Node constructor for the circularly linked list
node <- R6Class("Node",

                list(
                  val = NULL,
                  nxt = NULL,
                  prv = NULL,
                  initialize = function(val = NULL,
                                        nxt = NULL,
                                        prv = NULL){
                    self$val <- val
                    self$nxt <- nxt
                    self$prv <- prv
                  },

                  print = function(...){ # for debug
                    cat("Circle: \n")
                    cat("  label: ", self$val[[5]], "\n", sep = "")
                    cat("  x-coordinate: ", self$val[[2]], "\n", sep = "")
                    cat("  y-coordinate:  ", self$val[[3]], "\n", sep = "")
                    cat("  radius:  ", self$val[[6]], "\n", sep = "")
                    cat("  color: ", self$val$color, "\n", sep = "")
                    cat("  predessor: ", self$prv$val[[5]], "\n", sep = "")
                    cat("  successor: ", self$nxt$val[[5]], "\n", sep = "")
                    invisible(self)
                  }
                )
)

# miscalleneous/utility functions

#simple progress bar
progress_bar <- function (x, max = 100) {
  percent <- x / max * 100
  cat(sprintf('\r[%-50s] %d%%',
              paste(rep('=', percent / 2), collapse = ''),
              floor(percent)))
  if (x == max)
    cat('\n')
}

#area from radii vector (not even used)
areaFromRad <- function(c){
  ans <- c()
  for(i in c){
    area <- pi * (i^2)
    ans <- c(ans, area)
  }
  return(ans)
}

################ debugging functions ##############################################
#traverse and print every element of the list recursively.
traverse <- function(node){
  OV <- node$val[[5]]
  node$val[[5]] <- "ORIGIN"

  trav <- function(nd){
    if(nd$nxt$val[[5]]=="ORIGIN"){
      nd$nxt$val[[5]] <- OV
      return(nd)
    }
    print(nd)
    trav(nd[["nxt"]])
  }
  trav(node)
}

#getting length of circular linked list
clength <- function(node, showRT = FALSE){
  OV <- node$val$label
  ccc <- 1
  node$val$label <- "PKNBVDFGYHJNBVCXSDFGHJB" #if this is one of the names of the labels then RIP
  current <- node

  while(!identical(current$nxt$val$label,"PKNBVDFGYHJNBVCXSDFGHJB")){
    ccc <- ccc + 1
    current <- current[['nxt']]
  }
  node$val$label <- OV
  if (showRT) {
    print(starting-ending)
  }
  return(ccc)
}

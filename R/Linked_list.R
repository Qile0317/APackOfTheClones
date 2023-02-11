#script of the circular, doubly linked list implementation and also some utility functions for debugging


suppressPackageStartupMessages(library(dplyr))
library(R6)

# Node constructor for the circularly linked list
node <- R6::R6Class("Node",
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
                      }
                    )
)

# check if R6 so that R CMD doesnt bug me ab not using R6
invisible(R6::is.R6Class(1))

# miscalleneous/utility/debugging functions

#simple progress bar
progress_bar <- function (x, max = 100) {
  percent <- x / max * 100
  cat(sprintf('\r[%-50s] %d%%',
              paste(rep('=', percent * 0.5), collapse = ''),
              floor(percent)))
  if (x == max)
    cat('\n')
}

#traverse and print every element of the list recursively.
traverse <- function(node) {
  OV <- node$val[[5]]
  node$val[[5]] <- "ORIGIN"

  trav <- function(nd){
    if (nd$nxt$val[[5]] == "ORIGIN") {
      nd$nxt$val[[5]] <- OV
      return(nd)
    }
    print(nd)
    trav(nd[["nxt"]])
  }
  trav(node)
}

#getting length of circular linked list
clength <- function(node) {
  OV <- node$val$label
  ccc <- 1
  node$val$label <- "PKNBVDFGYHJNBVCXSDFGHJB" #if this is one of the names of the labels then RIP
  current <- node
  while (!identical(current$nxt$val$label,"PKNBVDFGYHJNBVCXSDFGHJB")){
    ccc <- ccc + 1
    current <- current[['nxt']]
  }
  node$val$label <- OV
  return(ccc)
}

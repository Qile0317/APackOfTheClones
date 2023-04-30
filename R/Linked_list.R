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

#simple progress bar
progress_bar <- function (x, max = 100) {
  percent <- x / max * 100
  cat(sprintf('\r[%-50s] %d%%',
              paste(rep('=', percent * 0.5), collapse = ''),
              floor(percent)))
  if (x == max)
    cat('\n')
}
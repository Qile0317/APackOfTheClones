#script of the circular, doubly linked list implementation and also some utility functions for debugging

#simple progress bar. In the future it should be possible to include ETA
progress_bar <- function (x, max = 100) {
  percent <- x / max * 100
  cat(sprintf('\r[%-50s] %d%%',
              paste(rep('=', percent / 2), collapse = ''),
              floor(percent)))
  if (x == max)
    cat('\n')
}

mut_list_el <- function(x) {
  nx <- as.character(deparse(substitute(x)))
  print(nx)
  globalenv()[[nx]][[1]][1] = 1
}
test()
a 

#function to mutate a list in place which doesnt work bc r 
mutate_list <- function(x, el, value, env = globalenv()) { 
  nx <- deparse(substitute(x))
  if (!exists(nx, env, mode = "list")) {
    stop("There is no list in ", sQuote("env"), " named ", sQuote(nx), ".")
  }
  if (!is.null(x[[el]])) {
    nv <- deparse(substitute(value))
    #cat("replacing value of", sQuote(el), "in", sQuote(nx), "with", sQuote(nv), "\n")
    env[[nx]][[el]] <- nv
  }
}

#mutate an element in the list.
mutate_list_el <- function(x, el, ind, value, env = globalenv()) { 
  nx <- deparse(substitute(x))
  if (!exists(nx, env, mode = "list")) {
    stop("There is no list in ", sQuote("env"), " named ", sQuote(nx), ".")
  }
  if (!is.null(x[[el]])) {
    nv <- deparse(substitute(value))
    nind <- deparse(substitute(ind))
    env[[nx]][[el]][nind] <- value
  }
}

#############################################################################################
# Node constructor
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
                  print = function(...){
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

################ unnecessary functions ##############################################
#traverse and print every element of the list.
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

#getting length of circular linked list iteratively
#note! does not work if its not circularly linked!! will run forever

clength <- function(node, showRT = FALSE){
  starting <- Sys.time()
  OV <- node$val$label
  ccc <- 1
  node$val$label <- "PKNBVDFGYHJNBVCXSDFGHJB" #if this is one of the names of the labels then RIP
  current <- node
  
  while(!identical(current$nxt$val$label,"PKNBVDFGYHJNBVCXSDFGHJB")){
    ccc <- ccc + 1
    current <- current[['nxt']]
  }
  node$val$label <- OV
  ending <- Sys.time()
  if(showRT==TRUE){
    print(starting-ending)
  }
  return(ccc)
}


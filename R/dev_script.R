# plotting for devs to test a single clusteer

plot_single_cluster_layout <- function(sizes) {
  a <- circle_layout(sizes)
  a <- df_full_join(list(a))
  plot_clusters(a)
}

plot_single_cluster <- function(cluster_list) {
  plot_clusters(df_full_join(list(cluster_list)))
}

# linked list debugging

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

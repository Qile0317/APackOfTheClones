getdata <- function(dir, name) {
	readRDS(paste("testdata/", dir, "/", name, ".rds", sep = ""))
}

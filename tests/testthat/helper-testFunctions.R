# helpers for data

getdata <- function(dir, name) {
	readRDS(paste("testdata/", dir, "/", name, ".rds", sep = ""))
}

sourcedata <- function(dir, name) {
	source(paste("testdata/", dir, "/", name, ".R", sep = ""))
}

# hacky helper for quietly testing
quietly <- function(e) suppressMessages(capture.output(e))

# helpers for data

getdata <- function(dir, name) {
	readRDS(paste("testdata/", dir, "/", name, ".rds", sep = ""))
}

sourcedata <- function(dir, name) {
	source(paste("testdata/", dir, "/", name, ".R", sep = ""))
}

# hacky trick for testing with no messages in the terminal

quietly <- function(e) suppressMessages(capture.output(e))

quietly_test_that <- function(desc, code) {
	test_that(desc, {quietly(code)})
}

# skippers

skip_if_r_version_leq <- function(version) {
	skip_if(is_curr_r_version_leq(version))
}

is_curr_r_version_leq <- function(version) {
	package_version(R.version) <= package_version(version)
}

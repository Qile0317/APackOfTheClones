# helpers for data

getdata <- function(dir, name) {
    readRDS(testthat::test_path("testdata", dir, paste0(name, ".rds")))
}

# this is honestly quite hacky, will rm in the future
sourcedata <- function(dir, name) {
    source(testthat::test_path("testdata", dir, paste0(name, ".R")))
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

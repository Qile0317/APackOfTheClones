# This dummy function definition is included with the package to ensure that
# 'tools::package_native_routine_registration_skeleton()' generates the required
# registration info for the 'run_testthat_tests' symbol
(function() {
  .Call("run_testthat_tests", FALSE, PACKAGE = "APackOfTheClones")
})

progress_bar <- function (x = 0, max = 100) {
  percent <- 100 * (x / max)
  cat(sprintf(
    '\r[%-50s] %d%%',
    paste(rep('=', percent * 0.5), collapse = ''),
    floor(percent))
  )
}

start_progress_bar <- function(verbose = TRUE) {
  if (verbose) {
    progress_bar(0, 1)
  }
}

end_progress_bar <- function(verbose = TRUE) {
  if (verbose) {
    progress_bar(1, 1)
  }
}

isnt_empty <- function(inp) {
  !identical(inp, list())
}

isnt_na <- function(inp) {
  !any(is.na(inp))
}

isnt_empty_nor_na <- function(inp) {
  isnt_empty(inp) && isnt_na(inp)
}

is_int <- function(num) {
  return(num == as.integer(num))
}

print_completion_time <- function(start_time, digits = 2) {
  end_time <- Sys.time()
  message(paste(
    "\nCompleted successfully, time elapsed:",
    round(as.numeric(end_time - start_time), digits),
    "seconds\n"
  ))
}

# helper function to correct reduction string
attempt_correction <- function(s) {
  s <- tolower(s)
  if (identical(s, "t-sne")) {
    s <- "tsne"
  }
  s
}
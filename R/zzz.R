.onAttach <- function(libname, pkgname) {
    packageStartupMessage("\nThank you for using APackOfTheClones! \n")
    packageStartupMessage("Please read the vignette found at https://qile0317.github.io/APackOfTheClones/articles/web_only/Clonal_expansion_plotting.html")
    packageStartupMessage("if you use this package in your paper, please cite this package with the citation() function\n")
    #packageStartupMessage(paste(
    #  "R package Version", 
    #  packageVersion("APackOfTheClones"),
    #  "https://CRAN.R-project.org/package=APackOfTheClones\n"
    #))
}

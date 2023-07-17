.onAttach <- function(libname, pkgname) {
    packageStartupMessage(paste(
        "\nThank you for using APackOfTheClones version",
        utils::packageVersion("APackOfTheClones"), "\n"
    ))
    packageStartupMessage(
        "Please read the vignette found at https://qile0317.github.io/APackOfTheClones/articles/web_only/Clonal_expansion_plotting.html"
    )
    packageStartupMessage(
        '\nif you use this package in your paper/study, please cite this package. Run `citation("APackOfTheClones")` for citation information\n'
    )
}

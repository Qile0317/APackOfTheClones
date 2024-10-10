.onAttach <- function(libname, pkgname) {
    packageStartupMessage(paste(
        "\nThank you for using APackOfTheClones v",
        utils::packageVersion("APackOfTheClones"), "\n\n",

        "if you use this package in your paper/study, please cite with\n",
        "`citation('APackOfTheClones')`\n",
        sep = ""
    ))
}

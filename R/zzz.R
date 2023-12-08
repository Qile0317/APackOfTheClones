.onAttach <- function(libname, pkgname) {
    packageStartupMessage(paste(
        "\nThank you for using APackOfTheClones v",
        utils::packageVersion("APackOfTheClones"), "\n\n",

        'if you use this package in your paper/study, please cite with\n',
        '`citation("APackOfTheClones")`\n',
        sep = ""
    ))

    # new deprecation note for the release
    packageStartupMessage(paste(
        '*** DEPRECATION NOTICE TO RETURNING USERS OF VERSION 0.1.x *** \n',
        deprecation_docstring(), "\n",
        sep = ""
    ))
}

.onAttach <- function(libname, pkgname) {
    packageStartupMessage(paste(
        "\nThank you for using APackOfTheClones version",
        utils::packageVersion("APackOfTheClones"), "\n\n",

        'if you use this package in your paper/study, please cite with',
        ' `citation("APackOfTheClones")`\n'
    ))

    # new deprecation note for the release
    packageStartupMessage(paste(
        'RETURNING USERS OF VERSION 0.1.x:\n',
        'the `clonal_expansion_plot()` function is deprecated. Please read',
        'the new vignette which will demonstrate a better way to achieve the',
        'same visualization.\n'
    ))
}

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

    # new deprecation note for the release
    packageStartupMessage(
    'NOTE FOR RETURNING USERS OF VERSION 0.1.x: the `clonal_expansion_plot()` function is being deprecated. Please read the new vignette which will demonstrate a better way to achieve the same visualization. Namely, the workflow now involves using the functions `RunAPOTC` and `APOTCPlot`. Please read their function level documentation as well.'
    )
}

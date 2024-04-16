#' @title
#' Highlight specific clones on an APackOfTheClones ggplot
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This is an analogue for `scRepertoire::highlightClones` that can highlight
#' certain clonotypes on an APackOfTheClones clonal expansion plot. For most
#' combinations of the arguments, there will be a ggplot fill legend on the
#' right side that correspond to each (existing) clonotype.
#'
#' @param apotc_ggplot A ggplot object that is the output of [APOTCPlot] or
#' [vizAPOTC] of an APackOfTheClones plot to be highlighted on.
#' @param clonotype character vector of the sequence(s) to highlight. Note
#' that it must be of the clonecall of the code that created the plot. A
#' warning will be shown if any of the sequences are not present.
#' @param color_each Either a logical of length 1, or a character(s). It is
#' `TRUE` by default, which assigns a unique default ggplot color to each
#' highlighted clone. If `FALSE`, each highlighted clone will retain its
#' current color and no legend based on color is shown. A possible application
#' here is to simply gauge the distribution of any shared clone.
#' It can also indicate the uniform color of each highlighted clone: if it is
#' a character of length 1, all highlighted clones will be of that color.
#' Else it must be a character vector of the same length as
#' `clonotype`, with each color corresponding to the clone. Here is a suitable
#' place to use any palette function from the many other CRAN palette packages
#' such as `"viridis"` or `"RColorBrewer"`. Note that currently, the user must
#' ensure `clonotype` contains only unique characters.
#' @param default_color A character of length 1 or `NULL` indicating the color
#' of non-highlighted clones. If `NULL`, all un-highlighted sequences will
#' retain their original color in `sc.data`. Else, if it is a character, it
#' should be a valid color that all un-highlighted clones are. Defaults to the
#' hexcode for gray.
#' @param scale_bg A positive numeric. Scales the brightness value of each color
#' of the non-highlighted clones by itself as a scaling factor. Defaults to 1
#' which will not alter the current brightness. Note that if
#' `color_each = FALSE` and `default_color = NULL`, this is equivalent to not
#' highlighting any clones - in this case, it may be useful to alter `scale_bg`
#' slightly so that the non-highlighted clones are darkened/brightened.
#' @param fill_legend logical indicating whether a ggplot legend of the "fill"
#' of each clonotype should be displayed.
#'
#' @details
#' Under the hood, this function simply mutates the plotting dataframe under
#' `$data` in the ggplot object, and operates on a column named `color`.
#'
#' Note that if `color_each = FALSE` and `default_color = NULL`, this is
#' equivalent to simply not highlighting anything and a warning will be shown.
#'
#' @return A ggplot object with the data modified to the highlighted colors
#' @export
#'
#' @examples
#' library(magrittr)
#' data("combined_pbmc")
#'
#' # piping the plot can be nice to read syntatically -
#' # By default, assigns unique colors to highlights and everything else is gray
#' vizAPOTC(combined_pbmc, clonecall = "aa", verbose = FALSE) %>%
#'     showCloneHighlight("CASLSGSARQLTF_CASSSTVAGEQYF")
#'
#' # one useful application is to highlight shared clones - beware that the
#' # clonotype sequences may get extremely long in the legend
#' shared_aa_clones <- names(getSharedClones(combined_pbmc, clonecall = "aa"))
#' vizAPOTC(combined_pbmc, clonecall = "aa", verbose = FALSE) %>%
#'     showCloneHighlight(shared_aa_clones)
#'
showCloneHighlight <- function(
    apotc_ggplot,
    clonotype, # could make into unlist(list(...))
    color_each = TRUE,
    default_color = "#808080",
    scale_bg = 1,
    fill_legend = TRUE
) {
    apotc_highlight_clones_error_handler()

    if (contains_duplicates(clonotype)) {
        warning(
            "`clonotype` contains duplicates - ",
            "this probably will result in erroneous outputs"
        )
    }

    if (is_false(color_each) && is.null(default_color) && scale_bg == 1) {
        warning(
            "setting `color_each` to FALSE, `default_color` to NULL without ",
            "altering `scale_bg` is equivalent to not highlighting anything"
        )
        return(apotc_ggplot)
    }

    num_seqs <- length(clonotype)
    sequence_index_map <- hash::hash(clonotype, 1:num_seqs)

    highlighted_ggplot_data <- get_ggplot_data(apotc_ggplot)

    clone_color_vector <- gen_clone_color_vector(
        color_each, clonotype, highlighted_ggplot_data
    ) # TODO be able to color by individual clone (not clonotype) for stuff like other parameters in the metadata

    num_matches <- 0

    for (i in seq_len(nrow(highlighted_ggplot_data))) {

        curr_clone <- highlighted_ggplot_data$clonotype[i]
        curr_seq_index <- sequence_index_map[[curr_clone]]

        if (is.null(curr_seq_index)) {
            if (!is.null(default_color)) {
                highlighted_ggplot_data$color[i] <- default_color
            }
            highlighted_ggplot_data$color[i] <- scale_hex_brightness(
                highlighted_ggplot_data$color[i], scale_bg
            )
            next
        }

        num_matches <- num_matches + 1
        if (identical(color_each, FALSE)) next
        highlighted_ggplot_data$color[i] <- clone_color_vector[curr_seq_index]
    }

    if (num_matches < num_seqs) {
        warning(
            ifelse(num_matches == 0, "all", "some"),
            " ", "input sequences didn't match any clone"
        )
    }

    apotc_ggplot <- set_ggplot_data(apotc_ggplot, highlighted_ggplot_data)

    # TODO technically probably possible to have seperate colors for each seq if color_each=FALSE
    if (identical(color_each, FALSE) || !fill_legend) return(apotc_ggplot)

    # FIXME fill colors doesnt show up on pkgdown though legend show up
    suppressMessages(apotc_ggplot + ggplot2::scale_fill_identity(
        guide = "legend",
        name = "clonotype",
        labels = clonotype,
        breaks = clone_color_vector
    ))

}

apotc_highlight_clones_error_handler <- function() {
    
    args <- get_parent_func_args()

    check_is_apotc_ggplot(args$apotc_ggplot)
    if (is_undetailed(args$apotc_ggplot)) {
        stop(call. = FALSE,
            "`apotc_ggplot` must be generated with `detail = TRUE`"
        )
    }

    typecheck(args$clonotype, is_character)
    typecheck(args$color_each, is_a_logical, is_a_character, is_character)
    typecheck(args$default_color, is_a_character, is.null)
    typecheck(args$scale_bg, is_a_positive_numeric)
    typecheck(args$fill_legend, is_a_logical)
}

gen_clone_color_vector <- function(color_each, sequence, plot_data) {

    if (identical(color_each, FALSE)) return(NULL)

    num_sequences <- length(sequence)

    if (is_a_character(color_each)) {
        return(rep(color_each, num_sequences))
    }

    if (is.character(color_each)) {
        if (length(color_each) != num_sequences)
            stop(call. = FALSE,
                "length of `color_each` doesn't match ",
                "the number of sequences to highlight"
            )
        return(color_each)
    }

    # assume color_each == TRUE
    gg_color_hue(num_sequences)

}

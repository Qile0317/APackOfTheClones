# APackOfTheClones's version of highlightClones

methods::setGeneric(
    "highlightClones",
    function(sc.data, ...) standardGeneric("highlightClones")
)

#' @title
#' Highlight specific clones on an APackOfTheClones ggplot
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' TODO - some combos and mention its S4 alt to scRep, and how there will be
#' a legend. 
#'
#' @param sc.data A ggplot object that is the output of [APOTCPlot] or
#' [vizAPOTC] of an APackOfTheClones plot to be highlighted
#' @param sequence character vector of the sequence(s) to highlight. Note
#' that it must be of the clonecall of the code that created the plot. A
#' warning will be shown if any of the sequences are not present.
#' @param color_each Either a logical of length 1, or a character(s). It is
#' `TRUE` by default, which assigns a unique default ggplot color to each
#' highlighted clone. If `FALSE`, each highlighted clone will retain its
#' current color in `sc.data`. It can also indicate the color of each
#' highlighted clone: if it is a character of length 1 and a valid color, all
#' highlighted clones will be of that color. Else it must be a character vector
#' of the same length as `sequence`, with each color corresponding to the
#' clone. Currently, the user must ensure `sequence` contains of unique
#' characters.
#' @param default_color A character of length 1 or `NULL` indicating the color
#' of non-highlighted clones. If `NULL`, all un-highlighted sequences will
#' retain their original color in `sc.data`. Else, if it is a character, it
#' should be a valid color that all un-highlighted clones are. Defaults to the
#' hexcode for gray.
#'
#' @details
#' TODO discuss that it modifies the ggplot data.
#'
#' @return A ggplot object with the data modified to the highlighted colors
#' @export
#'
#' @examples
#' data("combined_pbmc")
#'
#' # piping the plot can be nice to read syntatically -
#' # By default, assigns unique colors to highlights and everything else is gray
#' vizAPOTC(combined_pbmc, clonecall = "aa", verbose = FALSE) %>%
#'     APackOfTheClones::highlightClones("CASLSGSARQLTF_CASSSTVAGEQYF")
#'
#' # one useful application is to highlight shared clones - beware that the
#' # clonotype sequences may get extremely long in the legend
#' shared_aa_clones <- names(getSharedClones(combined_pbmc, clonecall = "aa"))
#' vizAPOTC(combined_pbmc, clonecall = "aa", verbose = FALSE) %>%
#'     APackOfTheClones::highlightClones(shared_aa_clones)
#'
methods::setMethod("highlightClones", "gg",
    function(
        sc.data, sequence, color_each = TRUE, default_color = "#808080"
    ) {
        apotc_highlight_clones(
            apotc_ggplot = sc.data,
            sequence = sequence,
            color_each = color_each,
            default_color = default_color
        )
    }
)

apotc_highlight_clones <- function(
    apotc_ggplot,
    sequence, # maybe allow indicies in seperate method??
    color_each, # in future allow palettes?
    default_color,
    add_legend = TRUE # pretty dumb if it was False?
) {

    apotc_highlight_clones_error_handler(hash::hash(as.list(environment())))

    # process data - probably should make sequences unique? for now assume unique

    num_seqs <- length(sequence)
    sequence_index_map <- hash::hash(sequence, 1:num_seqs)

    highlighted_ggplot_data <- get_ggplot_data(apotc_ggplot)

    clone_color_vector <- gen_clone_color_vector(
        color_each, sequence, highlighted_ggplot_data
    )

    num_matches <- 0

    for (i in seq_len(nrow(highlighted_ggplot_data))) {

        curr_clone <- highlighted_ggplot_data$clonotype[i]
        curr_seq_index <- sequence_index_map[[curr_clone]]

        if (is.null(curr_seq_index)) {
            if (!is.null(default_color)) {
                highlighted_ggplot_data$color[i] <- default_color
            }
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
    if (!add_legend || identical(color_each, FALSE)) return(apotc_ggplot)
    
    suppressMessages(
        apotc_ggplot + ggplot2::scale_fill_identity(
            guide = "legend",
            name = "clonotype",
            labels = sequence,
            breaks = clone_color_vector
        )
    )

}

apotc_highlight_clones_error_handler <- function(args) {

    check_is_apotc_ggplot(args$apotc_ggplot)

    if (!is.character(args$sequence) || (length(args$sequence) == 0)) {
        stop(call. = FALSE,
            "`sequence` must be a non-empty character vector"
        )
    }
    
    if (
        !(is_a_logical(args$color_each) || is.character(args$color_each))
            || (length(args$color_each) == 0)
    ) {
        stop(call. = FALSE,
            "`color_each` must be a length 1 logical or character, ",
            "or the character of `length(sequence)`"
        )
    }

    if (!is.null(args$default_color) && !is_a_character(args$default_color)) {
        stop(call. = FALSE,
            "`default_color` must be a character of length 1 or NULL"
        )
    }
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
                "length of `color_each` doesn't match",
                "the number of sequences to highlight"
            )
        return(color_each)
    }

    # assume color_each == TRUE
    gg_color_hue(num_sequences)

}

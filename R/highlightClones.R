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
#' TODO
#'
#' @return A ggplot object with the data modified to the highlighted colors
#' @export
#'
methods::setMethod("highlightClones", "gg",
    function(
        sc.data, sequence, color_each = FALSE, default_color = "#808080"
    ) {
        apotc_highlight_clones(
            apotc_ggplot = sc.data,
            sequence = sequence,
            color_each = color_each,
            default_color = default_color
        )
    }
)

# if color_each = FALSE, retain original color. If TRUE,
# unique colors. else custom color vector
# TODO check what happens when no sequence matches
apotc_highlight_clones <- function(
    apotc_ggplot,
    sequence,
    color_each,
    default_color,
    add_legend = TRUE
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

    if (!add_legend) return(apotc_ggplot)
    
    apotc_ggplot #%>% ggplot2::scale_fill_identity() # TODO
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
        stop(call. = FALSE, "`default_color` must be a character of length 1")
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
                "the numeber of sequences to highlight"
            )
        return(color_each)
    }

    # assume color_each == TRUE
    gg_color_hue(num_sequences)

}

#' Add NP-ROC curves to the current plot object.
#' @export
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @param x fitted NP-ROC object using \code{nproc}.
#' @param ... additional arguments.
#' @seealso \code{\link{npc}}, \code{\link{nproc}} and \code{\link{plot.nproc}}.

lines.nproc <- function(x, ...) {
    lines(x$typeI.u, 1 - x$typeII.u, type = "s", ...)
    if (x$band == TRUE) {
        lines(x$typeI.u, 1 - x$typeII.l, type = "s", ...)
    }

}


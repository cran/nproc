#' Add NP-ROC curves to the current plot object.
#' @export
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @param x fitted NP-ROC object using \code{nproc}.
#' @param col color of the added line.
#' @param ... additional arguments.
#' @seealso \code{\link{npc}}, \code{\link{nproc}} and \code{\link{plot.nproc}}.

lines.nproc <- function(x, col = "black", ...) {
    
    lines(x$typeI.u, 1 - x$typeII.u, xlab = "FPR", ylab = "TPR", type = "l", col = col)
    if (x$band == TRUE) {
        lines(x$typeI.u, 1 - x$typeII.l, col = col)
    }
    
}


#' Plot the nproc curve(s).
#' @export
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @param x fitted nproc object using \code{nproc}.
#' @param col color of the plots.
#' @param ... additional arguments.
#' @seealso \code{\link{npc}}, \code{\link{nproc}}
#' @examples
#' n = 1000
#' x = matrix(rnorm(n*2),n,2)
#' c = 1+3*x[,1]
#' y = rbinom(n,1,1/(1+exp(-c)))
#' fit = nproc(x, y, method = 'svm')
#' plot(fit)

plot.nproc <- function(x, col = "black", ...) {

plot(x$typeI.u, 1 - x$typeII.u, xlab = 'FPR', ylab = 'TPR', type = 'l', col = col)
if (x$band) {
lines(x$typeI.u, 1 - x$typeII.l, col = col)
}

}


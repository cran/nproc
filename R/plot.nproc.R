#' Plot the nproc curve(s).
#' @export
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @param x fitted nproc object using \code{nproc}.
#' @param ... additional arguments.
#' @seealso \code{\link{npc}} and \code{\link{nproc}}
#' @examples
#' n = 1000
#' x = matrix(rnorm(n*2),n,2)
#' c = 1+3*x[,1]
#' y = rbinom(n,1,1/(1+exp(-c)))
#' fit = nproc(x, y, method = 'svm')
#' plot(fit)
#' #fit = nproc(x, y, method = c('svm','logistic'))
#' #plot(fit)
#' ##Compare Confidence Curves
#' #fit = nproc(x, y, method = c('svm','logistic','lda'), conf = TRUE)
#' #plot(fit)
plot.nproc <- function(x, ...) {

    roc.lo = x$roc.lo
    roc.up = x$roc.up
    n.method = length(x$method)
    if (x$conf) {
        plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = paste("NP ROC: ",
            1 - x$delta, " Confidence Curve"))
    } else {
        plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR", ylab = "TPR", main = paste("NP ROC"))
    }

    for (i in 1:n.method) {
        range2 = which(x$cat==2)
        range1 = which(x$cat==1)
        range0 = which(x$cat==0)
        if(length(range2)>0) {
          range2 = c(range2,max(range2)+1)
        }
        if(length(range1)>0) {
          range1 = c(range1,max(range1)+1)
        }

        if(!x$adaptive){
          lines(x$alphalist, roc.lo[, 2 * i], type = "s", col = i)
        } else{
          lines(x$alphalist[range2], roc.lo[range2, 2 * i], lty = 2, type = "s", col = i)
          lines(x$alphalist[range1], roc.lo[range1, 2 * i], lty = 2, type = "s", col = i)
          lines(x$alphalist[range0], roc.lo[range0, 2 * i], type = "s", col = i)
        }
        if (x$conf) {
          if(!x$adaptive){
            lines(x$alphalist, roc.up[, 2 * i], type = "s", col = i)
          } else{
            lines(x$alphalist[range2], roc.up[range2, 2 * i], lty = 2, type = "s", col = i)
            lines(x$alphalist[range1], roc.up[range1, 2 * i], lty = 2, type = "s", col = i)
            lines(x$alphalist[range0], roc.up[range0, 2 * i], type = "s", col = i)
          }

        }
    }
    legend("bottomright", legend = x$method, col = 1:n.method, lty = rep(1, n.method))

}


#' Calculate the Neyman-Pearson ROC
#'
#' \code{nproc} calculate the Neyman-Pearson ROC
#' curve for a given sequence of type I error values.
#' @export
#' @param x n * p observation matrix. n observations, p covariates.
#' @param y n 0/1 observatons.
#' @param method classification method.
#' \itemize{
#' \item logistic: \link{glm} function with family = 'binomial'
#' \item penlog: \code{\link[glmnet]{glmnet}} in \code{glmnet} package
#' \item svm: \code{\link[e1071]{svm}} in \code{e1071} package
#' \item randomforest: \code{\link[randomForest]{randomForest}} in \code{randomForest} package
#' \item lda: \code{\link[MASS]{lda}} in \code{MASS} package
#' \item nb: \code{\link[e1071]{naiveBayes}} in \code{e1071} package
#' \item ada: \code{\link[ada]{ada}} in \code{ada} package
#' \item custom: a custom classifier. score vector needed.
#' }
#' @param score score vector corresponding to y. Required when method  = 'custom'.
#' @param alphalist the sequence of type I error values. Default = seq(from=0,to=1,by=0.01).
#' @param delta the violation rate of the type I error. Default = 0.05.
#' @param split whether the class 0 sample is split into two parts.
#' @param cv whether cross-validation is performed for calculating the roc curve.
#' @param fold number of folds for the cross-validation. Default = 5.
#' @param loc.prob the precalculated threshold locations in probability. Default = NULL.
#' @param n.cores number of cores used for parallel computing. Default = 1.
#' @seealso \code{\link{npc}}
#' @examples
#' n = 1000
#' x = matrix(rnorm(n*2),n,2)
#' c = 1+3*x[,1]
#' y = rbinom(n,1,1/(1+exp(-c)))
#' #fit = nproc(x, y, method = 'svm')
#' #fit2 = nproc(x, y, method = 'svm', cv = TRUE)
#' fit3 = nproc(x, y, method = 'penlog')
#' #fit3 = nproc(x, y, method = 'penlog',  n.cores = 2)
#' #In practice, replace 2 by the number of cores available 'detectCores()'
#' #fit4 = nproc(x, y, method = 'penlog', n.cores = detectCores())
#'
#' #Testing the custom method for nproc.
#' #fit = npc(x, y, method = 'lda', split = FALSE,  n.cores = 2) #use npc to get score list.
#' #obj = nproc(x = NULL, y = fit$y, method = 'custom', split = FALSE,
#' #score = fit$score,  n.cores = 2)




nproc <- function(x = NULL, y, method = c("logistic", "penlog", "svm", "randomforest",
                                          "lda", "nb", "ada", "custom"), score = NULL, alphalist = seq(from = 0.01, to = 0.99,
                                                                                                by = 0.01), delta = 0.05, split = TRUE, cv = FALSE, fold = 5, loc.prob = NULL,
                  n.cores = 1) {
  method = match.arg(method)
  alphalist = alphalist[alphalist > 0 & alphalist < 1]
  if (!cv) {
    fit = npc(x, y, method, score = score, alpha = alphalist, delta = delta,
              split = split, n.cores = n.cores)
    # if(is.null(loc.prob)) loc.prob = fit$loc.prob obj = npc.core(fit$y, fit$prob,
    # alpha = alphalist, delta = delta, loc.prob = loc.prob)
    v = getroc(fit$pred.y, fit$y)
  } else {
    n = length(y)
    n.fold = ceiling(n/fold)
    ind = sample(1:n)
    rocmat = array(0, c(fold, length(alphalist), 2))
    for (i in 1:fold) {
      cind = (i - 1) * n.fold + 1:n.fold
      cind = intersect(1:n, cind)
      te.ind = ind[cind]
      tr.ind = setdiff(1:n, te.ind)
      fit = npc(x[tr.ind, ], y[tr.ind], method, score = score[tr.ind], alpha = alphalist,
                delta = delta, split = split, loc.prob = loc.prob, n.cores = n.cores)
      if (is.null(loc.prob))
        loc.prob = fit$loc.prob
      pred = predict(fit, x[te.ind, ], score[te.ind])
      rocmat[i, , ] = getroc(pred$pred.label, y[te.ind])
    }
    v = apply(rocmat, c(2, 3), mean)

  }
  plot(alphalist, v[, 2], main = paste("NP ROC: ", method), xlab = "FPR", ylab = "TPR",
       type = "s")

  return(list = list(roc = v))

}

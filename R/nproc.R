#' Calculate the Neyman-Pearson Receiver Operating Characteristics
#'
#' \code{nproc} calculate the Neyman-Pearson Receiver Operating Characteristics
#' curve for a given sequence of type I error values.
#' @export
#' @param x n * p observation matrix. n observations, p covariates.
#' @param y n 0/1 observatons.
#' @param method classification method(s).
#' \itemize{
#' \item logistic: Logistic regression. \link{glm} function with family = 'binomial'
#' \item penlog: Penalized logistic regression with LASSO penalty. \code{\link[glmnet]{glmnet}} in \code{glmnet} package
#' \item svm: Support Vector Machines. \code{\link[e1071]{svm}} in \code{e1071} package
#' \item randomforest: Random Forest. \code{\link[randomForest]{randomForest}} in \code{randomForest} package
#' \item Linear Discriminant Analysis. lda: \code{\link[MASS]{lda}} in \code{MASS} package
#' \item nb: Naive Bayes. \code{\link[e1071]{naiveBayes}} in \code{e1071} package
#' \item ada: Ada-Boost. \code{\link[ada]{ada}} in \code{ada} package
#' \item custom: a custom classifier. score vector needed.
#' }
#' @param kernel kernel used in the svm method. Default = 'radial'.
#' @param score score vector corresponding to y. Required when method  = 'custom'.
#' @param pred.score score vector corresponding to the test y. Required when method  = 'custom'.
#' @param band whether to generate two np roc curves representing a confidence band. Default = FALSE.
#' @param typeI.lower whether to generate the data-driven type-I error lower bound. Default = FALSE.
#' @param delta the violation rate of the type I error. Default = 0.05.
#' @param split the number of splits for the class 0 sample. Default = 1. For ensemble
#' version, choose split > 1.  When method = 'custom',  split = 0 always.
#' @param split.ratio the ratio of splits used for the class 0 sample to train the
#' classifier. Default = 0.5.
#' @param n.cores number of cores used for parallel computing. Default = 1.
#' @param randSeed the random seed used in the algorithm.
#' @return An object with S3 class nproc.
#' \item{typeI.u}{sequence of upper bound of type I error.}
#' \item{typeII.l}{sequence of lower bound of type I error.}
#' \item{typeII.u}{sequence of upper bound of type II error.}
#' \item{auc.l}{the auc value of the lower NP-ROC curve.}
#' \item{auc.u}{the auc value of the upper NP-ROC curve.}
#' \item{band}{whether the upper NP-ROC curve is generated.}
#' \item{method}{the classification method implemented.}
#' \item{delta}{the violation rate.}
#' @seealso \code{\link{npc}}
#' @references
#' Xin Tong, Yang Feng, and Jingyi Jessica Li (2016), Neyman-Pearson (NP) classification algorithms and NP receiver operating characteristic (NP-ROC) curves, manuscript, http://arxiv.org/abs/1608.03109
#' @examples
#' n = 200
#' x = matrix(rnorm(n*2),n,2)
#' c = 1 - 3*x[,1]
#' y = rbinom(n,1,1/(1+exp(-c)))
#' #fit = nproc(x, y, method = 'svm')
#' fit2 = nproc(x, y, method = 'penlog')
#'
#' ##Plot the nproc curve
#' plot(fit2)
#' #fit3 = nproc(x, y, method = 'penlog')
#'
#' ##Plot the nproc curve
#' #plot(fit3)
#'
#' #fit3 = nproc(x, y, method = 'penlog',  n.cores = 2)
#' #In practice, replace 2 by the number of cores available 'detectCores()'
#' #fit4 = nproc(x, y, method = 'penlog', n.cores = detectCores())
#'
#' #Testing the custom method for nproc.
#' #fit = npc(x, y, method = 'lda', split = 0,  n.cores = 2) #use npc to get score list.
#' #obj = nproc(x = NULL, y = fit$y, method = 'custom', split = 0,
#' #score = fit$score,  n.cores = 2)
#'
#' #Confidence nproc curves
#' #fit6 = nproc(x, y, method = 'lda', band = TRUE)
#'
#' #nproc ensembled version
#' #fit7 = nproc(x, y, method = 'lda', split = 11)




nproc <- function(x = NULL, y, method = c("logistic", "penlog", "svm", "randomforest", "lda", "nb", "ada",
    "custom"), kernel = "radial", score = NULL, pred.score = NULL, band = FALSE, typeI.lower = FALSE, delta = 0.05, split = 1,
    split.ratio = 0.5, n.cores = 1, randSeed = 0) {
    if (!is.null(x)) {
        x = as.matrix(x)
    }
    p = ncol(x)
    if (p == 1 & method == "penlog") {
        stop("glmnet does not support the one predictor case. ")
    }
    method = match.arg(method)
    set.seed(randSeed)
    v = npc(x, y, method = method, band = band, kernel = kernel, score = score, pred.score = pred.score, delta = delta, split = split, split.ratio = split.ratio, n.cores = n.cores, randSeed = randSeed)

    split = v$split
    alphalist = seq(from = 0, to = 1, by = 0.001)
    n.alpha = length(alphalist)
    beta.u = beta.l = matrix(0,n.alpha,split)

    for(i in 1:split){
      obj = v$fits[[i]]

        beta.u[,i] = approx(obj$alpha.u, obj$beta.u, alphalist, method = 'constant', rule = 2, f = 0)$y


      if(v$band == TRUE){
        if(typeI.lower == 'FALSE'){
          beta.l[,i] = approx(obj$alpha.u, obj$beta.l, alphalist, method = 'constant', rule = 2, f = 1)$y
        }else{
            beta.u[,i] = approx(obj$alpha.l, obj$beta.l, alphalist, method = 'constant', rule = 2, f = 0)$y
          }
      }
      }
    beta.u.m = apply(beta.u, 1, mean)
    auc.l = sum(diff(alphalist) * (1-beta.u.m[-1]))
    auc.u = NULL
    beta.l.m = NULL
    if (v$band == TRUE){
    beta.l.m = apply(beta.l, 1, mean)
    auc.u = sum(diff(alphalist) * (1-beta.l.m[-1]))
    }



    object = list(typeII.u = beta.u.m, typeII.l = beta.l.m, auc.l = auc.l, auc.u = auc.u, band = band, method = method,
        typeI.u = alphalist, delta = delta)
    class(object) = "nproc"
    return(object)
}

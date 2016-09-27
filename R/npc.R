#' Calculate the Neyman-Pearson Classifier from a sample of class 0 and class 1.
#'
#' Give a type I error upper bound alpha and the violation upper bound delta, \code{npc} calculate the Neyman-Pearson Classifier
#' which control the type I error with in alpha with probability at least 1-delta
#' @export
#' @importFrom e1071 svm
#' @importFrom e1071 naiveBayes
#' @importFrom glmnet cv.glmnet
#' @importFrom MASS lda
#' @importFrom randomForest randomForest
#' @importFrom ada ada
#' @importFrom parallel mcmapply
#' @importFrom parallel mclapply
#' @importFrom graphics plot
#' @importFrom stats approx
#' @importFrom stats glm
#' @importFrom stats predict
#' @importFrom stats pbinom
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @importFrom graphics polygon
#' @param x n * p observation matrix. n observations, p covariates.
#' @param y n 0/1 observatons.
#' @param method classification method.
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
#' @param band  whether to generate the type II error lower and upper bounds. If TRUE, the class 1 data needs to be splitted into two parts like class 0. Default = FALSE.
#' @param kernel kernel used in the svm method. Default = 'radial'.
#' @param score score vector corresponding to y. Required when method  = 'custom'.
#' @param pred.score predicted score vector for the test sample. Optional when method  = 'custom'.
#' @param alpha the desirable control on type I error. Default = 0.05.
#' @param delta the violation rate of the type I error. Default = 0.05.
#' @param split the number of splits for the class 0 sample. Default = 1. For ensemble
#' version, choose split > 1.  When method = 'custom',  split = 0 always.
#' @param split.ratio the ratio of splits used for the class 0 sample to train the
#' classifier. Default = 0.5.
#' @param n.cores number of cores used for parallel computing. Default = 1. WARNING:
#' windows machine is not supported.
#' @param randSeed the random seed used in the algorithm.
#' @return An object with S3 class npc.
#' \item{fits}{a list of length split, represents the fit during each split.}
#'  \item{band}{whether the lower bound of type II erro is calculated.}
#'  \item{method}{the classification method.}
#'   \item{split}{the number of splits used.}
#' @seealso \code{\link{nproc}} and \code{\link{predict.npc}}
#' @references
#' Xin Tong, Yang Feng, and Jingyi Jessica Li (2016), Neyman-Pearson (NP) classification algorithms and NP receiver operating characteristic (NP-ROC) curves, manuscript, http://arxiv.org/abs/1608.03109.
#' @examples
#' set.seed(1)
#' n = 1000
#' x = matrix(rnorm(n*2),n,2)
#' c = 1+3*x[,1]
#' y = rbinom(n,1,1/(1+exp(-c)))
#' xtest = matrix(rnorm(n*2),n,2)
#' ctest = 1+3*xtest[,1]
#' ytest = rbinom(n,1,1/(1+exp(-ctest)))
#'
#' ##Use svm classifier and the default type I error control with alpha=0.05
#' fit = npc(x, y, method = 'svm')
#' pred = predict(fit,xtest)
#' fit.score = predict(fit,x)
#' accuracy = mean(pred$pred.label==ytest)
#' cat('Overall Accuracy: ',  accuracy,'\n')
#' ind0 = which(ytest==0)
#' typeI = mean(pred$pred.label[ind0]!=ytest[ind0]) #type I error on test set
#' cat('Type I error: ', typeI, '\n')
#'
#' ##Ensembled svm classifier with split = 11,  alpha=0.05
#' #fit = npc(x, y, method = 'svm', split = 11)
#' #pred = predict(fit,xtest)
#' #accuracy = mean(pred$pred.label==ytest)
#' #cat('Overall Accuracy: ',  accuracy,'\n')
#' #ind0 = which(ytest==0)
#' #typeI = mean(pred$pred.label[ind0]!=ytest[ind0]) #type I error on test set
#' #cat('Type I error: ', typeI, '\n')
#'
#' ##Now, change the method to logistic regression and change alpha to 0.1
#' fit = npc(x, y, method = 'logistic', alpha = 0.1)
#' pred = predict(fit,xtest)
#' accuracy = mean(pred$pred.label==ytest)
#' cat('Overall Accuracy: ',  accuracy,'\n')
#' ind0 = which(ytest==0)
#' typeI = mean(pred$pred.label[ind0]!=ytest[ind0]) #type I error on test set
#' cat('Type I error: ', typeI, '\n')
#'
#' ##Now, change the method to adaboost
#' #fit = npc(x, y, method = 'ada', alpha = 0.1)
#' #pred = predict(fit,xtest)
#' #accuracy = mean(pred$pred.label==ytest)
#' #cat('Overall Accuracy: ',  accuracy,'\n')
#' #ind0 = which(ytest==0)
#' #typeI = mean(pred$pred.label[ind0]!=ytest[ind0]) #type I error on test set
#' #cat('Type I error: ', typeI, '\n')
#'
#' ##A 'custom' npc classifier with y and score.
#' #fit2 = npc(y = y, score = fit.score$pred.score,
#' #pred.score = pred$pred.score, method = 'custom')

npc <- function(x = NULL, y, method = c("logistic", "penlog", "svm", "randomforest", "lda", "nb", "ada",
    "custom"), band = FALSE, kernel = "radial", score = NULL, pred.score = NULL, alpha = 0.05, delta = 0.05, split = 1,
    split.ratio = 0.5, n.cores = 1, randSeed = 0) {
    if (!is.null(x)) {
        x = as.matrix(x)
    }
    method = match.arg(method)
    set.seed(randSeed)
    if (method == "custom" & is.null(score)) {
        stop("score is needed when specifying method = \"custom\"")
    }
    if (!is.null(score)) {
        ## custom method, user specifed score vector
        test.score = score
        test.y = y
        pred.score = pred.score
        obj = npc.core(test.y, test.score, pred.score = pred.score, alpha = alpha, delta = delta)
        object = list(pred.y = obj$pred.y, y = test.y, score = test.score, cutoff = obj$cutoff, sign = obj$sign,
            method = method)
        class(object) = "npc"
        return(object)
    } else {
        ## not custom method
        object = NULL
        n = nrow(x)
        p = ncol(x)
        if (p == 1 & method == "penlog") {
            stop("glmnet does not support the one predictor case. ")
        }
        ind0 = which(y == 0)  ##indices for class 0
        ind1 = which(y == 1)  ##indices for class 1
        n0 = length(ind0)
        n1 = length(ind1)
        if (split == 0) {
            ## no split, use all class 0 obs for training and for calculating the cutoff
          fits = npc.split(x, y, p, kernel, alpha, delta, ind0, ind0, ind1, ind1, method, pred.score, n.cores = n.cores)
        } else {
            ## with split

                n0.1 = round(n0 * split.ratio)  ##default size for calculating the classifier
                n1.1 = round(n1 * split.ratio)

                  ind01.mat = sapply(1:split, f <- function(i) {
                    sample(ind0, n0.1)})
                  ind11.mat = sapply(1:split, f <- function(i) {
                      sample(ind1, n1.1)
                  })
            ind02.mat = sapply(1:split, f <- function(i) {
              setdiff(ind0, ind01.mat[, i])})
            ind12.mat = sapply(1:split, f <- function(i) {
                setdiff(ind1, ind11.mat[, i])
            })
            n0.cores = max(1,floor(n.cores/split))
            if (band == TRUE) {
              fits = mclapply(1:split, f <- function(i) {
                npc.split(x, y, p, kernel, alpha, delta, ind01.mat[, i], ind02.mat[, i], ind11.mat[,i], ind12.mat[,i], method,
                          pred.score, n.cores = n0.cores)
              }, mc.cores = n.cores)
            } else{
              fits = mclapply(1:split, f <- function(i) {
                npc.split(x, y, p, kernel, alpha, delta, ind01.mat[, i], ind02.mat[, i], ind1, ind1, method,
                          pred.score, n.cores = n0.cores)
              }, mc.cores = n.cores)
            }

        }
        object$fits = fits
        object$band = band
        object$split = split
        object$method = method
        class(object) = "npc"
        return(object)

    }


}

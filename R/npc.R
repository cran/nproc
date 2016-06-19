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
#' @param kernel kernel used in the svm method. Default = 'radial'.
#' @param score score vector corresponding to y. Required when method  = 'custom'.
#' @param pred.score predicted score vector for the test sample. Optional when method  = 'custom'.
#' @param alpha the desirable control on type I error. Default = 0.05.
#' @param delta the violation rate of the type I error. Default = 0.05.
#' @param split the number of splits for the class 0 sample. Default = 1. For ensemble
#' version, choose split > 1.  When method = 'custom',  split = 0 always.
#' @param split.ratio the ratio of splits used for the class 0 sample to train the
#' classifier. Default = 0.5.
#' @param n01.min the minimum number of class 0 observations used for training classifier. Default = 10.
#' @param adaptive whether to adaptively choose the number of observations used for
#' determining the cutoff when 1/2 n0 < lower bound < n0. Default = F.
#' @param n.cores number of cores used for parallel computing. Default = 1. WARNING:
#' windows machine is not supported.
#' @param randSeed the random seed used in the algorithm.
#' @return An object with S3 class npc.
#' \item{fit}{the fit from the specified classifier.}
#' \item{score}{the score vector for each observation.}
#'  \item{cutoff}{thecutoff determined via bootstrap to achieve the specified type I error
#'   control.}
#'  \item{sign}{whether class 1 has a larger average score than the
#'   class 0.}
#'  \item{method}{the classification method.}
#'  \item{loc.prob}{the
#'   percentile used to determine the cutoff for the specified type I error
#'   control.}
#'   \item{split}{the number of splits used.}
#'   \item{n0toosmall}{whether the class 0 sample size is too small for the desirable type
#'   I error control}
#' @seealso \code{\link{nproc}} and \code{\link{predict.npc}}
#' @examples
#' set.seed(0)
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
#' fit = npc(x, y, method = 'svm', split = 11)
#' pred = predict(fit,xtest)
#' accuracy = mean(pred$pred.label==ytest)
#' cat('Overall Accuracy: ',  accuracy,'\n')
#' ind0 = which(ytest==0)
#' typeI = mean(pred$pred.label[ind0]!=ytest[ind0]) #type I error on test set
#' cat('Type I error: ', typeI, '\n')
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
#' fit = npc(x, y, method = 'ada', alpha = 0.1)
#' pred = predict(fit,xtest)
#' accuracy = mean(pred$pred.label==ytest)
#' cat('Overall Accuracy: ',  accuracy,'\n')
#' ind0 = which(ytest==0)
#' typeI = mean(pred$pred.label[ind0]!=ytest[ind0]) #type I error on test set
#' cat('Type I error: ', typeI, '\n')
#'
#' ##A 'custom' npc classifier with y and score.
#' fit2 = npc(y = y, score = fit.score$pred.score,
#' pred.score = pred$pred.score, method = 'custom')

npc <- function(x = NULL, y, method = c("logistic", "penlog", "svm", "randomforest",
    "lda", "nb", "ada", "custom"), kernel = "radial", score = NULL, pred.score = NULL,
    alpha = 0.05, delta = 0.05, split = 1, split.ratio = 0.5, n01.min = 10, adaptive = FALSE, n.cores = 1, randSeed = 0) {
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
        obj = npc.core(test.y, test.score, pred.score = pred.score, alpha = alpha,
            delta = delta)
        object = list(pred.y = obj$pred.y, y = test.y, score = test.score, cutoff = obj$cutoff,
            sign = obj$sign, method = method)
        class(object) = "npc"
        object$n0toosmall = 'FALSE'
        return(object)
    } else { ##not custom method
        object = NULL
        n = nrow(x)
        p = ncol(x)
        ind0 = which(y == 0)  ##indices for class 0
        ind1 = which(y == 1)  ##indices for class 1
        n0 = length(ind0)

        n02.min = ceiling(log(delta)/log(1-min(alpha)))
        if(n0 < n01.min * (split>0) + n02.min) { ##too few class 0 obs, note the
          ##difference when split=0 and >0
          object$n0toosmall = TRUE
          warning('too few class 0 observations to ensure a type I error control.')
          class(object) = 'npc'
          return(object)
        } else if(split == 0){ ##no split, use all class 0 obs for training and for calculating the cutoff
            object = npc.split(x, y, p, kernel, alpha, delta, ind0, ind0, ind1, method,
                pred.score)
        } else { ##with split

            if(adaptive == F){
              n0.1 = round(n0*split.ratio) ##default size for calculating the classifier
              if(n0 < n0.1+n02.min){
              object$n0toosmall = TRUE
              warning('too few class 0 observations to ensure a type I error control.')
              class(object) = 'npc'
              object$n0toosmall = 'TRUE'
              return(object)
            } else{
              ind01.mat = sapply(1:split, f <- function(i) {
                sample(ind0, n0.1)
            })
            }
            }else {
              ind01.mat = sapply(1:split, f <- function(i) {
                sample(ind0, n0-n02.min)
              })
            }


            ind02.mat = sapply(1:split, f <- function(i) {
                setdiff(ind0, ind01.mat[, i])
            })
            object = mclapply(1:split, f <- function(i) {
                npc.split(x, y, p, kernel, alpha, delta, ind01.mat[, i], ind02.mat[,
                  i], ind1, method, pred.score)
            }, mc.cores = n.cores)
        }
        object$split = split
        object$method = method
        class(object) = "npc"
        object$n0toosmall = 'FALSE'
        return(object)

    }


}

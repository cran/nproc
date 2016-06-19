# Calculate the Neyman-Pearson Classifier from a vector of scores. In general, the
# classifier works as follows. The training sample is divided into two parts, one
# part is a mixed sample of class 0 and class 1, the other part is just a sample of
# class 0.  @param y the true response for the training observations. A vector of
# 0/1 values, with value 0/1 representing that the observation is from class 0/1.
# @param score the score list for a set of training observations. An example would
# be the predicted success probabilities for a group of patients.  @param
# pred.score prediction score list for a set of testing observations.  @param alpha
# the desirable control on type I error. Default = 0.05.  @param delta the
# violation rate of the type I error. Default = 0.05.  Default = NULL.  @param
# n.cores number of cores used for parallel computing. Default = 1.  @seealso
# \code{\link{nproc}}, \code{\link{predict.npc}}, and \code{\link{nproc}}.
npc.core <- function(y, score, pred.score = NULL, alpha = 0.05, delta = 0.05) {

    ind0 = which(y == 0)
    ind1 = which(y == 1)
    n0 = length(ind0)

    # loc = bs(1:n0, method = 'bs', alpha, delta, n.cores = n.cores)
    loc = find.order(1:n0, alpha, delta)
    if(length(ind0)==0 || length(ind1)==0){
      stop('both class 0 and class 1 responses are needed to decide the cutoff.')
    }
    sig = mean(score[ind1]) > mean(score[ind0])  ##whether the class 0 has a larger average score than class 1
    if (sig == FALSE) {
        loc = n0 + 1 - loc
    }
    cutoff = (sort(score[ind0]))[loc]
    names(cutoff) = alpha
    if (is.null(pred.score)) {
        pred.score = score
    }
    if (sig == TRUE) {
        pred.y = outer(pred.score, cutoff, ">")
    } else {
        pred.y = outer(pred.score, cutoff, "<=")
    }
    return(list = list(cutoff = cutoff, sign = sig, pred.y = pred.y))
}

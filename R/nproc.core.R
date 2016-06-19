nproc.core <- function(x = NULL, y, method = c("logistic", "penlog", "svm", "randomforest",
    "lda", "nb", "ada", "custom"), score = NULL, conf = FALSE, alphalist = seq(from = 0,
    to = 1, by = 0.01), delta = 0.05, split = 1, split.ratio = 0.5, n01.min = 10, adaptive = FALSE, cv = FALSE,
    fold = 5, n.cores = 1) {
    method = match.arg(method)
        cat = rep(0,length(alphalist))

    if (!cv) {
        roc = matrix(0,length(alphalist),2)
        n0 = sum(y==0)
        n02.min = n0-n01.min*(split>0)
        n02.min.split = ceiling(n0*(1-split.ratio*(split>0)))
        range1 = which(alphalist < 1 - delta^(1/n02.min))
        range2 = which(1 - delta^(1/n02.min)<=alphalist & alphalist < 1-delta^(1/n02.min.split))
        range3 = which(alphalist >1-delta^(1/n02.min.split))
        cat[range1] = 2 ##cannot control type I error
        cat[range2] = 1 ##can control with more than the default split of obs
        if(adaptive){
        for(alpha.ind in range2){
          fit = npc(x, y, method, score = score, alpha = alphalist[alpha.ind], delta = delta,
                    split = split, split.ratio = split.ratio, n01.min = n01.min, adaptive = adaptive, n.cores = n.cores)
          pred = predict(fit, x, score)
          roc[alpha.ind,] = getroc(pred$pred.label, y)
          }
        }
        fit = npc(x, y, method, score = score, alpha = alphalist[range3], delta = delta,
                 split = split, split.ratio = split.ratio, n01.min = n01.min, adaptive = adaptive, n.cores = n.cores)
        pred = predict(fit, x, score)
        roc[range3,] = getroc(pred$pred.label, y)
    } else {
        n = length(y)
        ind0 = which(y==0)
        ind1 = which(y==1)
        ind = c(ind0,ind1)
        x = x[ind,]
        y = y[ind] ##reorder data according to y values
        n0 = sum(y==0)
        n1 = sum(y==1)
        n0.fold = ceiling(n0/fold) ##stratified cross-validation
        n1.fold = ceiling(n1/fold)
        n.fold = ceiling(n/fold)
        ind0 = sample(1:n0)
        ind1 = sample(1:n1)
        n0.tr = floor(n0/fold*(fold-1))
        n02.min = n0.tr-n01.min*(split>0)
        n02.min.split = ceiling(n0.tr*(1-split.ratio*(split>0)))
        range1 = which(alphalist < 1 - delta^(1/n02.min))
        range2 = which(1 - delta^(1/n02.min)<=alphalist & alphalist < 1-delta^(1/n02.min.split))
        range3 = which(alphalist >1-delta^(1/n02.min.split))
        cat[range1] = 2 ##cannot control type I error
        cat[range2] = 1 ##can control with more than the default split of obs
        rocmat = mcmapply(1:fold, FUN = function(i) {
          cind0 = (i - 1) * n0.fold + 1:n0.fold
          cind0 = intersect(1:n0, cind0)
          cind1 = (i - 1) * n1.fold + 1:n1.fold
          cind1 = intersect(1:n1, cind1)
            te.ind = c(ind0[cind0],ind1[cind1]+n0)
            tr.ind = setdiff(1:n, te.ind)
            a.cores = max(1, floor(n.cores/fold))
            roc = matrix(0,length(alphalist),2)
               if(adaptive){
              for(alpha.ind in range2){
                fit = npc(x[tr.ind,], y[tr.ind], method, score = score, alpha = alphalist[alpha.ind], delta = delta,
                          split = split, split.ratio = split.ratio, n01.min = n01.min, adaptive = adaptive, n.cores = a.cores)
                pred = predict(fit, x[te.ind,], score[te.ind])
                roc[alpha.ind,] = getroc(pred$pred.label, y[te.ind])
              }
            }
            fit = npc(x[tr.ind,], y[tr.ind], method, score = score, alpha = alphalist[range3], delta = delta,
                      split = split, split.ratio = split.ratio, n01.min = n01.min, adaptive = FALSE, n.cores = a.cores)
            # if(is.null(loc.prob)) loc.prob = fit$loc.prob obj = npc.core(fit$y, fit$prob,
            # alpha = alphalist, delta = delta, loc.prob = loc.prob)
            pred = predict(fit, x[te.ind,], score[te.ind])
            roc[range3,] = getroc(pred$pred.label, y[te.ind])
            roc
        }, mc.cores = fold)
        rocmat = array(rocmat, c(length(alphalist), 2, 5))

        # for(i in 1:fold) { cind = (i - 1) * n.fold + 1:n.fold cind = intersect(1:n, cind)
        # te.ind = ind[cind] tr.ind = setdiff(1:n, te.ind) fit = npc(x[tr.ind, ],
        # y[tr.ind], method, score = score[tr.ind], alpha = alphalist, delta = delta, split
        # = split, loc.prob = loc.prob, n.cores = n.cores) pred = predict(fit, x[te.ind, ],
        # score[te.ind]) rocmat[i, , ] = getroc(pred$pred.label, y[te.ind]) }
        roc = apply(rocmat, c(1, 2), mean)

    }
     return(list(v = roc, cat = cat))
}


pred.npc.core <- function(object, newx = NULL) {
    if (object$method == "logistic") {
        pred.score = predict(object$fit, data.frame(newx), type = "response")
    } else if (object$method == "penlog") {
        pred.score = predict(object$fit$glmnet.fit, newx = newx, type = "response",
            s = object$fit$lambda.min)
    } else if (object$method == "svm") {
        pred.score = attr(predict(object$fit, newx, decision.values = TRUE), "decision.values")[,
            1]
    } else if (object$method == "randomforest") {
        pred.score = predict(object$fit, newx, type = "prob")[, 2]
    } else if (object$method == "lda") {
        pred.score = predict(object$fit, data.frame(newx))$posterior[, 2]
    } else if (object$method == "nb") {
        pred.score = predict(object$fit, data.frame(newx), type = "raw")[, 2]
    } else if (object$method == "ada") {
        pred.score = predict(object$fit, data.frame(newx), type = "probs")[, 2]
    }
    pred.score = as.vector(pred.score)
    if (object$sign == TRUE) {
        pred.label = outer(pred.score, object$cutoff, ">")
    } else {
        pred.label = outer(pred.score, object$cutoff, "<=")
    }
    list(pred.label = pred.label, pred.score = pred.score)

}

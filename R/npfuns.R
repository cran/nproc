npc.split <- function(x, y, p, kernel, alpha, delta, ind01, ind02, ind1, method, pred.score = pred.score) {
    train.ind = c(ind01, ind1)
    test.ind = c(ind02, ind1)
    train.x = as.matrix(x[train.ind, ])
    train.y = y[train.ind]
    test.x = as.matrix(x[test.ind, ])
    test.y = y[test.ind]
    colnames(train.x) = paste("x", 1:p, sep = "")
    colnames(test.x) = paste("x", 1:p, sep = "")
    if (method == "logistic") {
        train.data = data.frame(train.x, y = train.y)
        fit = glm(y ~ ., data = train.data, family = "binomial")
        test.score = predict(fit, data.frame(test.x), type = "response")
    } else if (method == "penlog") {
        fit = cv.glmnet(train.x, train.y, family = "binomial")
        test.score = predict(fit$glmnet.fit, newx = test.x, type = "response", s = fit$lambda.min)
        test.score = as.vector(test.score)
    } else if (method == "svm") {
        train.y = as.factor(train.y)
        fit = svm(train.x, train.y, kernel = kernel)
        test.score = attr(predict(fit, test.x, decision.values = TRUE), "decision.values")[,
            1]
    } else if (method == "randomforest") {
        train.y = as.factor(train.y)
        fit = randomForest(train.x, train.y)
        test.score = predict(fit, test.x, type = "prob")[, 2]
    } else if (method == "lda") {
        fit = lda(train.x, train.y)
        test.score = predict(fit, data.frame(test.x))$posterior[, 2]
    } else if (method == "nb") {
        fit = naiveBayes(train.x, train.y)
        test.score = predict(fit, data.frame(test.x), type = "raw")[, 2]
    } else if (method == "ada") {
        train.data = data.frame(train.x, y = train.y)
        fit = ada(y ~ ., data = train.data)
        test.score = predict(fit, data.frame(test.x), type = "probs")[, 2]
    }

    obj = npc.core(test.y, test.score, pred.score = pred.score, alpha = alpha, delta = delta)

    object = list(fit = fit, pred.y = obj$pred.y, y = test.y, score = test.score, cutoff = obj$cutoff,
        sign = obj$sign, method = method)
    class(object) = "npc"
    return(object)
}
getroc <- function(yhat, y) {
    TPR = apply(yhat, 2, f <- function(x) {
        sum((x == y) & (y == 1))/sum(y == 1)
    })
    FPR = apply(yhat, 2, f <- function(x) {
        sum((x != y) & (y == 0))/sum(y == 0)
    })
    cbind(FPR, TPR)
}

get.vio <- function(roc.th, roc.em, nproc.em, alphalist) {
    alphaall = sort(c(roc.th[, 1], roc.em[, 1]))
    t1 = approx(roc.th[, 1], roc.th[, 2], alphaall, rule = 2)$y
    t2 = approx(roc.em[, 1], roc.em[, 2], alphaall, rule = 2)$y
    loc = which(t2 > t1)
    vio1 = sum(diff(alphaall)[loc[loc < length(alphaall)]])
    alphaall2 = sort(c(roc.th[, 1], alphalist))
    t3 = approx(roc.th[, 1], roc.th[, 2], alphaall2, rule = 2)$y
    t4 = approx(alphalist, nproc.em[, 2], alphaall2, rule = 2)$y
    loc = which(t4 > t3)
    vio2 = sum(diff(alphaall)[loc[loc < length(alphaall2)]])
    return(c(vio1, vio2))
}

getalpha <- function(roc1, roc2) {
    alphaall = sort(c(roc1[, 1], roc2[, 1]))
    t1 = approx(roc1[, 1], roc1[, 2], alphaall, rule = 2)$y
    t2 = approx(roc2[, 1], roc2[, 2], alphaall, rule = 2)$y
    alphaall[which.min(diff(sign(t1 - t2) != 0)) + 1]
}






## find the order such that the type-I error bound is satisfied with probability at
## least 1-delta
find.order <- function(prob, alpha = 0.05, delta = 0.05) {
    ### prob: probability of assigning to 1 method: alpha: a vector of type-I error
    ### delta: violation rate control
    n = length(prob)
    n.alpha = length(alpha)
    order = rep(n, n.alpha)
    if (sum(alpha > 1 - delta^(1/n)) == 0) {
        return(order)
    } else {
        im = min(which(alpha > 1 - delta^(1/n)))
        for (i in im:n.alpha) {
            order[i] = min(which(1 - pbinom(0:(n - 1), n, 1 - alpha[i]) < delta))
        }
        # min.r.alpha = 1-log(delta)/log(n) r.alpha = alpha[alpha>min.r.alpha]
        return(order)

    }
}



## try replace bootstrap with subsampling interpolation
bs <- function(prob, method = c("bs", "ss"), alpha = c(0.01, 0.05), delta = 0.05, B = 1000,
    n.cores = 1) {
    ### prob: probability of assigning to 1 method: bs: boostrap, ss: subsampling n/2 of
    ### n alpha: a vector of type-I error delta: violation rate control B: number of
    ### bootstraps set.seed(0)
    n = length(prob)
    cutoff.list = as.vector(sort(prob))
    if (method == "bs") {
        cutoff.bs = matrix(cutoff.list[sample(x = 1:n, size = n * B, replace = TRUE)],
            ncol = B)  # store the bootstrapped lr.Hat into a matrix. Every column is a bootstrap sample
    } else if (method == "ss") {
        cutoff.bs = matrix(0, round(n/2), B)
        for (i in 1:B) {
            cutoff.bs[, i] = cutoff.list[sample(1:n, size = round(n/2), replace = FALSE)]
        }
    }
    if (n.cores > 1) {
        percents <- mcmapply(cutoff.list, FUN = function(cutoff) {
            # for every possible cutoff, find the percentage of empirical type I errors > alpha
            temp1 = cutoff.bs - cutoff
            temp2 = matrix(as.numeric(temp1 > 0), ncol = B)
            temp3 = colSums(temp2)/n
            apply(outer(temp3, alpha, ">"), 2, sum)/B
        }, mc.cores = n.cores)
    } else {
        percents <- sapply(cutoff.list, FUN = function(cutoff) {
            # for every possible cutoff, find the percentage of empirical type I errors > alpha
            temp1 = cutoff.bs - cutoff
            temp2 = matrix(as.numeric(temp1 > 0), ncol = B)
            temp3 = colSums(temp2)/n
            apply(outer(temp3, alpha, ">"), 2, sum)/B
        })
    }

    if (is.matrix(percents) == FALSE) {
        percents = matrix(percents, 1, length(percents))
    }
    ind1vec = apply(percents, 1, f <- function(x) {
        min(which(x <= delta))
    })
    # ind2vec = apply(percents, 1, f <- function(x) { max(which(x > delta)) })
    w1 = abs(percents[ind1vec] - delta)
    # w2 = abs(percents[ind2vec] - delta) cutoff = (w1* cutoff.list[ind2] + w2*
    # cutoff.list[ind1])/(w1+w2)

    cutoff = cutoff.list[ind1vec]

    return(cutoff)

}

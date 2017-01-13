npc.core <- function(y, score, alpha = NULL, delta = 0.05, n.cores = 1) {

    ind0 = which(y == 0)
    ind1 = which(y == 1)
    if (length(ind0) == 0 || length(ind1) == 0) {
        stop("both class 0 and class 1 responses are needed to decide the cutoff.")
    }
    sig = mean(score[ind1]) > mean(score[ind0])  ##whether the class 0 has a larger average score than class 1
    if (sig == FALSE) {
        score = -score
    }
    # loc = bs(1:n0, method = 'bs', alpha, delta, n.cores = n.cores)

    score0 = sort(score[ind0])
    score1 = sort(score[ind1])
    obj = find.order(score0, score1, delta, n.cores = n.cores)

    cutoff.list = obj$scores


    cutoff = NULL
    if (!is.null(alpha)) {
        loc = max(which(obj$alpha.u <= alpha))
        if (loc > 1) {
            cutoff = cutoff.list[loc]
        } else {
            cutoff = Inf
        }
    }
    return(list = list(cutoff = cutoff, loc = loc, sign = sig, alpha.l = obj$alpha.l, alpha.u = obj$alpha.u,
        beta.l = obj$beta.l, beta.u = obj$beta.u))
}

pred.npc.core <- function(object, newx = NULL) {
    if (object$method == "logistic") {
        pred.score = predict(object$fit, data.frame(newx), type = "response")
    } else if (object$method == "penlog") {
        pred.score = predict(object$fit$glmnet.fit, newx = newx, type = "response", s = object$fit$lambda.min)
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
    if (object$sign == FALSE) {
        pred.score = -pred.score
    }
    pred.label = outer(pred.score, object$cutoff, ">")
    list(pred.label = pred.label, pred.score = pred.score)

}


npc.split <- function(x, y, p, kernel, alpha, delta, ind01, ind02, ind11, ind12, method,
    n.cores = 1) {
    train.ind = c(ind01, ind11)
    test.ind = c(ind02, ind12)
    train.x = x[train.ind, , drop = FALSE]
    train.y = y[train.ind]
    test.x = x[test.ind, , drop = FALSE]
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
    obj = npc.core(test.y, test.score, alpha = alpha, delta = delta, n.cores = n.cores)

    object = list(fit = fit, y = test.y, score = test.score, cutoff = obj$cutoff,
        sign = obj$sign, method = method, beta.l = obj$beta.l, beta.u = obj$beta.u, alpha.l = obj$alpha.l,
        alpha.u = obj$alpha.u)
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






## find the order such that the type-I error bound is satisfied with probability at
## least 1-delta
find.order <- function(score0, score1 = NULL, delta = 0.05, n.cores = 1) {


    score0 = sort(score0, decreasing = TRUE)
    score1 = sort(score1, decreasing = TRUE)
    n0 = length(score0)
    n1 = length(score1)
    sd0 = sd(score0)
    sd1 = sd(score1)

    score0 = score0 + rnorm(n0, 0, sd0/100)
    score1 = score1 + rnorm(n1, 0, sd1/100)

    # scores = sort(c(score0,score1), decreasing = TRUE)
    scores = score0
    # n = n0 + n1 nc = n + 1 ###number of different classifiers
    nc = n0 + 1  ###number of different classifiers
    alpha.l = alpha.u = beta.l = beta.u = rep(0, nc)
    alpha.l[1] = alpha.u[1] = 0
    beta.l[1] = beta.u[1] = 1
    alpha.l[nc] = alpha.u[nc] = 1
    beta.l[nc] = beta.u[nc] = 0
    v.list = seq(from = 0, to = 1, by = 0.001)
   # ru0 = apply(outer(score0, scores, "<="), 2, sum)
   # rl1 = apply(outer(score1, scores, "<="), 2, sum)
   # rl0 = n0 + 1 - apply(outer(score0, scores, ">="), 2, sum)
   # ru1 = n1 + 1 - apply(outer(score1, scores, ">="), 2, sum)
    ru0 = mcmapply(function(s){sum(score0 <= s)}, scores, mc.cores = n.cores)
    rl1 = mcmapply(function(s){sum(score1 <= s)}, scores, mc.cores = n.cores  )
    rl0 = n0 + 1 - mcmapply(function(s){sum(score0 >= s)}, scores, mc.cores = n.cores)
    ru1 = n1 + 1 - mcmapply(function(s){sum(score1 >= s)}, scores, mc.cores = n.cores)

    alpha.mat = mcmapply(f <- function(i) {
        if (ru0[i] == 0) {
            alpha.u = 1
            alpha.l = 1
        } else if (ru0[i] == n0) {
            alpha.l = 0
            alpha.u = v.list[min(which(pbinom(n0 - ru0[i], n0, v.list) <= delta))]
        } else {
            alpha.u = v.list[min(which(pbinom(n0 - ru0[i], n0, v.list) <= delta))]
            alpha.l = v.list[max(which(pbinom(n0 - rl0[i], n0, v.list) >= 1 - delta))]
        }
        if (rl1[i] == n1) {
            beta.u = 1
            beta.l = v.list[max(which(pbinom(rl1[i] - 1, n1, v.list) <= delta))]
        } else if (rl1[i] == 0) {
            beta.l = 0
            beta.u = 0
        } else {
            # beta.u[i] = v.list[which.max(1-pbinom(ru1[i], n1, v.list) >= 1-delta)]
            beta.u = v.list[min(which(pbinom(ru1[i] - 1, n1, v.list) <= delta))]
            beta.l = v.list[max(which(pbinom(rl1[i] - 1, n1, v.list) >= 1 - delta))]
        }
        c(alpha.l, alpha.u, beta.l, beta.u)
    }, 2:n0, mc.cores = n.cores)
    alpha.l[2:n0] = alpha.mat[1, ]
    alpha.u[2:n0] = alpha.mat[2, ]
    beta.l[2:n0] = alpha.mat[3, ]
    beta.u[2:n0] = alpha.mat[4, ]


    return(list(scores = scores, beta.l = beta.l, beta.u = beta.u, alpha.l = alpha.l, alpha.u = alpha.u))

}




rm(list=ls())
library(glmnet)
library(arm)
library(lattice)
# library(numDeriv)
# library(h2o)
source("FISTA.R")

# generate test data
seed = 1234
set.seed(seed)

nsims = 1000
scale.x.true = FALSE

# shape.true = exp(rnorm(1))
shape.true = 1000
m = 1000 # num of examples
n = 15 # num of variables
numzero = 10
alpha.EN = 0.5
ABSTOL = 1e-8
maxiter = 1e4
smallest.idx = 5

# true solution with sparsity
if(scale.x.true){
        x.true = rnorm(n, sd = sqrt(1/max(m,n)))
} else {
        x.true = rnorm(n)
}
# x.true = rnorm(n)
# x.true = log(x.true + min(c(x.true, 0)) + 0.0001)
zero.idx = sample(1:n, numzero)
x.true[zero.idx] = 0

# generate A matrix and b

A = matrix(rnorm(m * n, sqrt(1/max(m,n))), nrow = m) # independent variables
A = scale(A, center = TRUE)
# A = log(A + min(c(A, 0)) + 0.0001)
# A = log(A)
# scales = exp(A %*% x.true) / shape.true
# b = sapply(scales, function(x) rgamma(1, shape = shape.true, scale = x))  # dependent variable, A %*% x.true + noise
# b = c(scale(b, center = FALSE))

b_true <- exp(A %*% x.true)

reslist = list()

for(sim.iter in 1:nsims){       

b <- rgamma(m, rate = shape.true / b_true, shape = shape.true)

# nll_gamma_optim(c(x.true, shape.true), A, b)
# f(x.true)

# set alpha.EN and lambda.EN

# lambda.EN.max = norm(t(A) %*% b, "i")
# lambda.EN.max = find_lambda_max_Gamma_GLM_EN(A, b, alpha.EN, shape.true)
# lambda_vec = generate_lambda_grid(A, b, alpha.EN, shape.true)
# lambda.EN = lambda_vec[50]
# lambda.EN = 0

# lambda.EN = 0

# use FISTA implementation

# f = function(x){
#         nll_gamma_glm(x, A, b, shape.true)
# }
# 
# f_grad = function(x){
#         nll_gamma_glm_grad(x, A, b, shape.true)
# }
# 
# f_grad_num = function(x){
#         grad(f, x)
# }
# 
# g = function(x){
#         regularizer_EN(x, lambda.EN = lambda.EN, alpha.EN = alpha.EN)
# }
# 
# g_prox = function(x, lambda){
#         prox_EN(x, lambda = lambda, lambda.EN = lambda.EN, alpha.EN = alpha.EN)
# }
# 
# 
# res.FISTA = FISTA(rep(0, n), f, g, f_grad, g_prox, ABSTOL = ABSTOL, maxiter = maxiter)
# res.df = cbind(x.true, res.FISTA$x)
# apply(res.df, 2, function(x) sum(x==0))
# res.df


# res.glm.optim = fit_glm_gamma(c(rnorm(n), exp(rnorm(1))), A, b)
# c(tail(res.glm.optim$par,1), shape.true)
# res.glmnet = glmnet(A, b, alpha = alpha.EN, lambda = lambda.EN, standardize = FALSE, intercept = FALSE)
# res.glmnet$beta
# res.FISTA.numgrad = FISTA(rep(0, n), f, g, f_grad_num, g_prox, ABSTOL = ABSTOL, maxiter = maxiter)
# res.FISTA.numgrad
m_glm <- glm(b ~ A - 1, family = Gamma(link = "log"))
coef(m_glm)
# (res.optim = optim(rep(0, ncol(A)), function(x) f(x) + g(x), method = "BFGS"))
# cbind(x.true, res.FISTA$x, res.optim$par)
# # f_grad_num(rep(0, ncol(A)))
# # f_grad(rep(0, ncol(A)))
# res.optim$value
# tail(res.FISTA$objvals,1)
# res.glmGammaNet = glmGammaNet(rep(0, ncol(A)), A, b, lambda.EN, shape.true)
# cbind(x.true, res.glmGammaNet$x)

res.cv.glmGammaNet = cv.glmGammaNet(A, b, alpha.EN = alpha.EN, nlambda = 100)

# fit glmGammaNet using the max lambda with NLL that is less than smallest.idx 

best.idx.percentile = max(sort(res.cv.glmGammaNet$NLL_vec, index.return = TRUE)$ix[1:smallest.idx])

best.percentile.lambda.EN = res.cv.glmGammaNet$lambda.EN.vec[best.idx.percentile]

best.percentile.fit = glmGammaNet(res.cv.glmGammaNet$x, A, b,
                       lambda.EN = best.percentile.lambda.EN,
                       shape0 = res.cv.glmGammaNet$shape,
                       alpha.EN = alpha.EN,
                       ABSTOL = ABSTOL,
                       maxiter = maxiter)

# fit glmGammaNet using max lambda with NLL less than NLL_vec[best.idx] + NLL_sds[best.idx]

best.idx = res.cv.glmGammaNet$best.idx
NLL_threshold = res.cv.glmGammaNet$NLL_vec[best.idx] + res.cv.glmGammaNet$NLL_sds[best.idx]
best.idx.1sd = max(which(res.cv.glmGammaNet$NLL_vec < NLL_threshold))

best.1sd.fit = glmGammaNet(res.cv.glmGammaNet$x, A, b,
                                  lambda.EN = res.cv.glmGammaNet$lambda.EN.vec[best.idx.1sd],
                                  shape0 = res.cv.glmGammaNet$shape,
                                  alpha.EN = alpha.EN,
                                  ABSTOL = ABSTOL,
                                  maxiter = maxiter)

#plot(res.cv.glmGammaNet$lambda.EN.vec, res.cv.glmGammaNet$NLL_vec, type = "l")

res.cv.glmnet = cv.glmnet(A, b, alpha = alpha.EN, standardize = FALSE, intercept = FALSE)
res.cv.glmnet$cvm
###plot(res.cv.glmnet$cvm)
###plot(res.cv.glmnet)
lambda.glmnet = res.cv.glmnet$lambda[which.min(res.cv.glmnet$cvm)]
res.glmnet = glmnet(A, b, alpha = alpha.EN, lambda = lambda.glmnet, standardize = FALSE, intercept = FALSE)
# res.glmnet$beta
# x.true
# (res.coeffs = data.frame(x.true = x.true, x.glm = coef(m_glm), x.glmGammaNet = res.cv.glmGammaNet$x, x.glmnet = as.vector(res.glmnet$beta)))
# res.cv.glmGammaNet$best.idx
# nll_gamma_optim(c(x.true, shape.true), A, b)
# nll_gamma_optim(c(res.cv.glmGammaNet$x, res.cv.glmGammaNet$shape), A, b)
# apply(res.coeffs, 2, function(x) sum( abs(x) <= 1e-8))
# res.cv.glmGammaNet$NLL_vec[res.cv.glmGammaNet$best.idx]
# NLL_vec = res.cv.glmGammaNet$NLL_vec
# best.idx = res.cv.glmGammaNet$best.idx
# NLL_sds = res.cv.glmGammaNet$NLL_sds
# which(NLL_vec <= NLL_vec[best.idx] + NLL_sds[best.idx])
#plot(res.cv.glmGammaNet$NLL_vec, type = "l")
# x.tmp = glmGammaNet(rep(0, ncol(A)), A, b, lambda.EN = res.cv.glmGammaNet$lambda.EN.vec[82], res.cv.glmGammaNet$shape)
# cbind(x.true, x.tmp$x)
# res.cv.glmGammaNet$x
# best.idx
(res.coeffs = data.frame(x.true = x.true, 
                         x.glm = coef(m_glm), 
                         x.glmGammaNet = res.cv.glmGammaNet$x, 
                         x.glmGammaNet.percentile = best.percentile.fit$x,
                         x.glmGammaNet.1sd = best.1sd.fit$x,
                         x.glmnet = as.vector(res.glmnet$beta)))
reslist[[sim.iter]] = res.coeffs
}
reslist

# compute the L1 norm difference

L1Norm = function(x){
        sum(abs(x))
}

computeL1Error = function(x, percentage = FALSE){
        res = NULL
        for(i in 2:ncol(x)){
                res = c(res, sum(abs((x[,i] - x[,1]))))
        }
        if(percentage)
                res = res / sum(abs(x[,1]))
        res
}

error.L1 = sapply(reslist, computeL1Error, percentage = FALSE)
str(error.L1)
apply(error.L1, 1, mean)
histogram(error.L1[1,])
histogram(error.L1[2,])
histogram(error.L1[3,])

# compute the number of correct zero coefficients

computeCorrectZeros = function(x, true.zero.idx, tol = 1e-8){
        res = NULL
        for(i in 2:ncol(x)){
                idx = which(abs(x[,i]) <= tol)
                nzero = 0
                for(i in idx){
                        if (i %in% true.zero.idx)
                                nzero = nzero + 1
                }
                res = c(res, nzero)
        }
        res
}

# reslist[[1]]
# computeCorrectZeros(reslist[[1]], zero.idx)


zeros.correct = sapply(reslist, computeCorrectZeros, true.zero.idx = zero.idx)
# histogram(zeros.correct[1,])
# histogram(zeros.correct[2,])
# histogram(zeros.correct[3,])

# compute the number of incorrect zero coefficients

computeIncorrectZeros = function(x, true.zero.idx, tol = 1e-8){
        res = NULL
        for(i in 2:ncol(x)){
                idx = which(abs(x[,i]) <= tol)
                nzero = 0
                for(i in idx){
                        if (!i %in% true.zero.idx)
                                nzero = nzero + 1
                }
                res = c(res, nzero)
        }
        res
}

# reslist[[1]]
# computeIncorrectZeros(reslist[[1]], zero.idx)


zeros.incorrect = sapply(reslist, computeIncorrectZeros, true.zero.idx = zero.idx)
# histogram(zeros.incorrect[1,])
# histogram(zeros.incorrect[2,])
# histogram(zeros.incorrect[3,])

# compute the zero score (correct zeros - incorrect zeros)

zeros.score = zeros.correct - zeros.incorrect
# apply(zeros.score, 1 , mean)
# histogram(zeros.score[1,])
# histogram(zeros.score[2,])
# histogram(zeros.score[3,])

# performance metrics
res.df = data.frame(error.L1 = apply(error.L1, 1, mean),
                    zeros.correct = apply(zeros.correct, 1 , mean),
                    zeros.incorrect = apply(zeros.incorrect, 1, mean),
                    zeros.score = apply(zeros.score, 1 , mean)
)
rownames(res.df) = c("glmGamma", "glmGammaNet", "glmGammaNet.percentile", "glmGammaNet.1sd", "glmnet")

filename = paste0( seed, "_", m, "_", n, "_", smallest.idx, "_", numzero, "_", nsims, "_", shape.true, "_", scale.x.true)

num_nonzeros = sapply(reslist, function(x) apply(x, 2, function(x) sum(abs(x) > 1e-8)))
discrete.histogram(num_nonzeros["x.glmGammaNet",])
discrete.histogram(num_nonzeros["x.glmGammaNet.percentile",])
discrete.histogram(num_nonzeros["x.glmGammaNet.1sd",])
reslist[[1]]
rownames(num_nonzeros)
write.csv(res.df, paste0("./data/",filename, ".csv"))
save.image(paste0("./data/", filename, ".Rdata"))

# histogram(zeros.score[2,], breaks = -3:3)
# hist(zeros.score[2,], breaks = -3:3)
# barplot(zeros.score[2,])
# discrete.histogram(zeros.score[2,], bar.width = 0.6)
# ?discrete.histogram
# png(paste0("./plots/", filename, ".png"), 1024, 768)
# discrete.histogram(zeros.score[2,], bar.width = 0.6)
# dev.off()


# Histogram of Zeros
reslist[[1]]
res.df


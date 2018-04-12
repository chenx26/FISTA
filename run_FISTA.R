rm(list=ls())
library(glmnet)
source("FISTA.R")

# generate test data

set.seed(123)

m = 100 # num of examples
n = 5  # num of variables


# true solution with sparsity
x.true = rnorm(n)
x.true[2] = 0
x.true[4] = 0

# generate A matrix and b

A = matrix(rnorm(m * n), nrow = m) # independent variables
b = A %*% x.true + rnorm(m, 0, 0.01)    # dependent variable, A %*% x.true + noise

# set alpha.EN and lambda.EN
alpha.EN = 0.5
lambda.EN.max = norm(t(A) %*% b, "i")
lambda.EN = 0.0001 * lambda.EN.max
ABSTOL = 1e-5
maxiter = 1e3
# lambda.EN = 0

# use FISTA implementation

f = function(x){
        mse(x, A, b)/length(b)
}

f_grad = function(x){
        mse_grad(x, A, b)/length(b)
}

g = function(x){
        regularizer_EN(x, lambda.EN = lambda.EN, alpha.EN = alpha.EN)
}

g_prox = function(x, lambda){
        prox_EN(x, lambda = lambda, lambda.EN = lambda.EN, alpha.EN = alpha.EN)
}

res.FISTA = FISTA(rep(0, n), f, g, f_grad, g_prox, ABSTOL = ABSTOL, maxiter = maxiter)
res.FISTA
res.glmnet = glmnet(A, b, alpha = alpha.EN, lambda = lambda.EN, standardize = FALSE, intercept = FALSE)
res.glmnet$beta
cbind(x.true, res.FISTA$x, res.glmnet$beta)

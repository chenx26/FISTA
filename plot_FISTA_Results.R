rm(list=ls())
library(arm)
library(knitr)
load("./data/1234_200_6_1000_1_FALSE.Rdata")
x.true
zero.idx.true = which(x.true == 0)
tmp = reslist[[1]]

# plot histogram of number of zeros

count_zeros = function(x){
        apply(x, 2, function(x) sum(x == 0))
}
zeros.df = sapply(reslist, count_zeros)
dim(zeros.df)
discrete.histogram(zeros.df[3,], bar.width = 0.6)

# create table of number of times subsets of 0 up to 3 zero coefficients identified
i = 1
subset.counts = NULL
for(i in 1:length(reslist)){
        tmp = reslist[[i]]
        subset.counts = rbind(subset.counts, apply(tmp, 2, function(x) sum(which(x == 0) %in% zero.idx.true)))
}
subset.counts.mat = matrix(0, nrow = 4, ncol = ncol(subset.counts))
rownames(subset.counts.mat) = colnames(subset.counts)
for(i in 1:ncol(subset.counts)){
        for(j in 1:nrow(subset.counts)){
                subset.counts.mat[i,subset.counts[j, i]+1] =  subset.counts.mat[i,subset.counts[j, i]+1] + 1
        }
}
colnames(subset.counts.mat) = 0:3
kable(subset.counts.mat)

# create table of number of times each zero is correctly identified

zero.variable.counts = matrix(0, nrow = 4, ncol = 3)
colnames(zero.variable.counts) = c("z1", "z2", "z3")
rownames(zero.variable.counts) = colnames(subset.counts)

i = 1
for(i in 1:length(reslist)){
        tmp = reslist[[i]]
        for(j in 1:length(zero.idx.true)){
                for(k in 1:ncol(tmp)){
                        if(tmp[zero.idx.true[j], k] == 0){
                                zero.variable.counts[k, j] = zero.variable.counts[k, j] + 1
                        }
                }
        }
}
kable(subset.counts.mat)
kable(zero.variable.counts)

x = runif(5)
x
sort(x)
sort(x, index.return = TRUE)$ix[3]

## demonstrates that nipals is better than BPCA for the metabolite data???
library(pcaMethods)

data(metaboliteData)

met <- scale(t(metaboliteData), center=TRUE, scale=FALSE)

gone <- sample(1:7280, 7280 * 0.05)

goal <- met[gone]

met[gone] <- NA

pcNip <- pca(met, scale="none", center=TRUE, method="nipa", nPcs=5)
pcBpc <- pca(met, scale="none", center=TRUE, method="bpca", nPcs=5)

nipPred <- (pcNip@scores %*% t(pcNip@loadings))[gone]
bpcPred <- (pcBpc@scores %*% t(pcBpc@loadings))[gone]

## squared error
sum((goal - nipPred)^2)
sum((goal - bpcPred)^2)

## visible in plot as well
par(mfrow=c(1,2))
plot(goal, nipPred)
plot(goal, bpcPred)

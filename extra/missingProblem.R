## demonstrates that nipals is better than BPCA for the metabolite data
library(pcaMethods)

data(metaboliteData)
met <- t(metaboliteData)
met[sample(1:7280, 7280 * 0.15)] <- NA

goal <- t(metaboliteData)[is.na(met)]

pcNip <- pca(met, scale="uv", center=TRUE, method="nipa", nPcs=5)
pcBpc <- pca(met, scale="uv", center=TRUE, method="bpca", nPcs=5)

nipPred <- fitted(pcNip)[wasna(pcNip)]
bpcPred <- fitted(pcBpc)[wasna(pcBpc)]

## squared error
sum((goal - nipPred)^2, na.rm=TRUE)
sum((goal - bpcPred)^2, na.rm=TRUE)

## visible in plot as well
par(mfrow=c(1,2))
plot(goal, nipPred)
abline(c(0,1), lwd=2, col="blue")
plot(goal, bpcPred)
abline(c(0,1), lwd=2, col="blue")

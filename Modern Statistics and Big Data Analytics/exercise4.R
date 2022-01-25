library(MASS)
library(pdfCluster)
library(fpc)
library(cluster)
library(mclust)
library(flexclust)
library(EMMIXskew)


### 1
p05 <- bundestag(2005)

p05_mvn <- EmSkew(p05, g=7, distr="mvn")
EmSkew.contours(p05, p05_mvn)

p05_mvt <- EmSkew(p05, g=7, distr="mvt")

p05_msn <- EmSkew(p05, g=7, distr="msn")

p05_mst <- EmSkew(p05 ,g=7, distr="mst")


### 2
data(oliveoil)
olive <- oliveoil[, 3:10]

kolive <- kmeans(olive, 3, nstart=100)
molive <- Mclust(olive, G=3, "EII")

adjustedRandIndex(molive$classification, kolive$cluster)

data(geyser)

kgeyser <- kmeans(geyser, 3, nstart=100)
mgeyser <- Mclust(geyser, G=3, "EII")

adjustedRandIndex(mgeyser$classification, kgeyser$cluster)



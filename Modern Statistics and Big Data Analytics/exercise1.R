library(MASS)
library(pdfCluster)
library(fpc)

data(oliveoil)
str(oliveoil)
olive <- oliveoil[,3:10]
solive <- scale(olive)

clusterdata2 <- as.matrix(read.table("clusterdata2.txt"))
str(clusterdata2)
sclusterdata <- scale(clusterdata2)

### K means
set.seed(12345)

olivesk3 <- kmeans(solive, 3, nstart = 100)
plot(solive, col=olivesk3$cluster, pch=clusym[olivesk3$cluster])

olivesk9 <- kmeans(solive, 9, nstart = 100)
plot(solive, col=olivesk9$cluster, pch=clusym[olivesk9$cluster])

### gap statistic method
nco <- clusGap(solive, kmeans, 20, B = 100, d.power = 2, spaceH0 = "original", nstart = 100)
print(nco, method = "Tibs2001SEmax")

ncc <- clusGap(clusterdata2, kmeans, 20, B = 100, d.power = 2, spaceH0 = "original", nstart = 100)
print(ncc, method = "Tibs2001SEmax")
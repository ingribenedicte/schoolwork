library(MASS)
library(pdfCluster)
library(fpc)
library(cluster)
library(mclust)
library(StatMatch)


veronica <- read.table("veronica.txt")
veronicam <- as.matrix(veronica)
p <- nrow(veronicam)

data(geyser)
str(geyser)

### Multidimensional Scaling Graphs

jveronica <- dist(veronica, method="binary")
smveronica <- dist(veronica, method="manhattan")/p

mdsjveronica <- cmdscale(jveronica, k = 2)
mdssmveronica <- cmdscale(smveronica, k = 2)

plot(mdsjveronica)
plot(mdssmveronica)

dgeyser2 <- dist(geyser, method="euclidean")
dgeyser1 <- dist(geyser, method="manhattan")
dgeyserm <- matrix(0, ncol = 299, nrow = 299)
geysercov <- cov(geyser)
mgeyser <- as.matrix(geyser)
for (i in 1:299)
  dgeyserm[i,] <- mahalanobis(mgeyser, mgeyser[i,], geysercov)

mdsdgeyser2 <- cmdscale(dgeyser2, k = 2)
mdsdgeyser1 <- cmdscale(dgeyser1, k = 2)
mdsdgeyserm <- cmdscale(dgeyserm, k = 2)

plot(mdsdgeyser2)
plot(mdsdgeyser1)
plot(mdsdgeyserm)

### Gower dissimilarity

x1 <- c('blue', 1, 1, 12)
x2 <- c('red', 0, NA, NA)
x3 <- c('red', 1, 0, 17)
x <- cbind(x1, x2, x3)

dx <- gower.dist(x)
head(dx)

### K-means

data(oliveoil)
olive <- oliveoil[, 3:10]

olive1 <- olive[c(0,2:8)]
olive2 <- olive[c(1,3:8)]
olive3 <- olive[c(1:2,4:8)]
olive4 <- olive[c(1:3,5:8)]
olive5 <- olive[c(1:4,6:8)]
olive6 <- olive[c(1:5,7:8)]
olive7 <- olive[c(1:6,8)]
olive8 <- olive[c(1:7)]

ARIs_us <- vector(, 8)

olivek <- kmeans(olive, 9, nstart = 100)

olive1k <- kmeans(olive1, 9, nstart = 100)
ARIs_us[1] <-adjustedRandIndex(olive1k$cluster,olivek$cluster)
olive2k <- kmeans(olive2, 9, nstart = 100)
ARIs_us[2] <-adjustedRandIndex(olive2k$cluster,olivek$cluster)
olive3k <- kmeans(olive3, 9, nstart = 100)
ARIs_us[3] <-adjustedRandIndex(olive3k$cluster,olivek$cluster)
olive4k <- kmeans(olive4, 9, nstart = 100)
ARIs_us[4] <-adjustedRandIndex(olive4k$cluster,olivek$cluster)
olive5k <- kmeans(olive5, 9, nstart = 100)
ARIs_us[5] <-adjustedRandIndex(olive5k$cluster,olivek$cluster)
olive6k <- kmeans(olive6, 9, nstart = 100)
ARIs_us[6] <-adjustedRandIndex(olive6k$cluster,olivek$cluster)
olive7k <- kmeans(olive7, 9, nstart = 100)
ARIs_us[7] <-adjustedRandIndex(olive7k$cluster,olivek$cluster)
olive8k <- kmeans(olive8, 9, nstart = 100)
ARIs_us[8] <-adjustedRandIndex(olive8k$cluster,olivek$cluster)

solive <- scale(olive)

solive1 <- solive[c(0, 2:8)]
solive2 <- solive[c(1, 3:8)]
solive3 <- solive[c(1:2, 4:8)]
solive4 <- solive[c(1:3, 5:8)]
solive5 <- solive[c(1:4, 6:8)]
solive6 <- solive[c(1:5, 7:8)]
solive7 <- solive[c(1:6, 8)]
solive8 <- solive[c(1:7)]

ARIs_s <- vector(, 8)

solivek <- kmeans(solive, 9, nstart = 100)

solive1k <- kmeans(solive1, 9, nstart = 100)
ARIs_s[1] <-adjustedRandIndex(solive1k$cluster,solivek$cluster)
solive2k <- kmeans(solive2, 9, nstart = 100)
ARIs_s[2] <-adjustedRandIndex(solive2k$cluster,solivek$cluster)
solive3k <- kmeans(solive3, 9, nstart = 100)
ARIs_s[3] <-adjustedRandIndex(solive3k$cluster,solivek$cluster)
solive4k <- kmeans(solive4, 9, nstart = 100)
ARIs_s[4] <-adjustedRandIndex(solive4k$cluster,solivek$cluster)
solive5k <- kmeans(solive5, 9, nstart = 100)
ARIs_s[5] <-adjustedRandIndex(solive5k$cluster,solivek$cluster)
solive6k <- kmeans(solive6, 9, nstart = 100)
ARIs_s[6] <-adjustedRandIndex(solive6k$cluster,solivek$cluster)
solive7k <- kmeans(solive7, 9, nstart = 100)
ARIs_s[7] <-adjustedRandIndex(solive7k$cluster,solivek$cluster)
solive8k <- kmeans(solive8, 9, nstart = 100)
ARIs_s[8] <-adjustedRandIndex(solive8k$cluster,solivek$cluster)
library(MASS)
library(flexclust)
library(fpc)
library(pdfCluster)
library(mclust)
library(cluster)
library(dbscan)

### Hierarchical clustering
veronica <- read.table("veronica.txt")
jveronica <- dist(veronica, method="binary")

plantclusts <- hclust(jveronica, method="single")
plantclustc <- hclust(jveronica, method="complete")
plantclusta <- hclust(jveronica, method="average")

ARI_avg <- vector(, 10)

for (i in 1:10) {
  K = 2*i
  plantclusts_K <- cutree(plantclusts, K)
  plantclustc_K <- cutree(plantclustc, K)
  plantclusta_K <- cutree(plantclusta, K)
  
  ARIsc <- adjustedRandIndex(plantclusts_K, plantclustc_K)
  ARIsa <- adjustedRandIndex(plantclusts_K, plantclusta_K)
  ARIca <- adjustedRandIndex(plantclustc_K, plantclusta_K)
  
  ARI_avg[i] = (ARIsc + ARIsa + ARIca)/3
}


### PAM clustering

# PAM
pamclust <- pam(jveronica,8)
pamclust[4]

### DBSCAN
p05 <- bundestag(2005)

eps <- 1
minPts <- 25

dbscan(p05, eps, minPts, weights = NULL, borderPoints = TRUE)
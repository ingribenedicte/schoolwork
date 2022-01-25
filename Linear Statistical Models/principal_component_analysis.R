#Principle component analysis

#a
X <- na.omit(USArrests)
X

# using prcomp to get loadings:
pca <- prcomp(X,scale = TRUE)
pca$rotation

# manually: normalized eigenvectors are loadings
s <- cor(X)
eigen(s)


#b
#sample variances are eigenvalues
eigen(s)$values

#standard deviation squared using prcomp:
pca$sdev^2


#c
#The i'th principal component of X is the i'th loading times X
t(pca$rotation) %*% t(scale(USArrests))[, "Alabama"]


#d
plot(pca$x[,1], pca$x[,2], type="n")
text(pca$x[,1], pca$x[,2], rownames(USArrests), cex=0.6)
biplot(pca, scale=0, cex=0.6)
plot(scale(USArrests)[, 2], scale(USArrests)[, 4])
abline(0, pca$rotation[2, 1] / pca$rotation[4, 1], col = "red")

#e
eigen(s)$values
#M = 2
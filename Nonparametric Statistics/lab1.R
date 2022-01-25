### Exercise 1
x = 8
n = 10
p_0 = 0.5

binom.test(x, n, p_0, alternative = "g")

### Exercise 2
x = 80
n = 90
p_0 = 0.75

binom.test(x, n, p_0, alternative = "t")

### Exercise 3
nh = c(20,22,13,22,10,13)
n = sum(nh)
H = length(nh)
ph_star = rep(1/H,H)
nh_star = n*ph_star

chi = sum((nh-nh_star)^2/(nh_star))

p_value <- pchisq(chi, H-1, ncp= 0, lower.tail = FALSE)
p_value

### Exercise 4
nh = c(13,9,6,5,7,10)
n = sum(nh)
H = length(nh)

lambda_hat = mean(nh)
x = 0:11
phstar = dpois(x,lambda_hat)

phstar_new = matrix(NA,6,1)
phstar_new [1] <- sum(phstar[1:4])
phstar_new [2] <- phstar[5]
phstar_new [3] <- phstar[6]
phstar_new [4] <- phstar[7]
phstar_new [5] <- phstar[8]
phstar_new [6] <- 1- sum(phstar_new[1:5])

chi = sum((nh-n*phstar_new)^2/n*phstar_new)

p_value = pchisq(chi,H-1, ncp = 0, lower.tail = TRUE)
p_value

### Exercise 5
x <- matrix(1:12, nrow = 3, dimnames = list(c("null","partial","total"), c("1","2","3","4")))
x

C1 = c(5,12,6)
C2 = c(4,10,11)
C3 = c(7,15,7)
C4 = c(8,14,9)

X2 = cbind(C1,C2,C3,C4)

chi = chisq.test(X2)
chi

### Exercise 6
C = matrix(nrow = 2, c(2,7,8,7))
fisher.test(C)

### Exercise 7
C = matrix(nrow = 2, c(28,7,13,27))
mcnemar.test(C)
### Exercise 2
x <- data.frame(appreciation=c(1,3,4,5,4,2,4,3,1,5,4,1,2,4,3),
                song=rep(c(1,2,3),each=5),
                subject=rep(c(1,2,3,4,5),times=3))

friedman.test(appreciation ~ song | subject, data = x)

### Exercise 3
theta_0 = 300

x = c(235, 230, 180, 250, 280, 330, 440, 430, 260, 225, 240, 235, 215)
x.sel <- x[x!=theta_0]

sign.pos <- x.sel-theta_0 > 0

binom.test(sum(sign.pos),length(sign.pos),0.5,alternative="t")

### Exercise 4
x <- matrix(nrow=2,c(5,2,1,7))
fisher.test(x,alternative="greater")

### Exercise 5
n=10
K=4

x <- data.frame(concentration=c(82,89,95,88,81,79,68,67,58,66,10,28,39,41,68,38,25,59,27,44,
                                19,26,50,60,65,74,51,49,51,57,84,95,99,85,92,78,83,79,87,64),
                condition=rep(c(1,2,3,4),each=n),
                subject=rep(c(1,2,3,4,5,6,7,8,9,10),times=K))

friedman.test(concentration ~ condition | subject, data=x)

### Exercise 6
x <- matrix(nrow=2,c(50,5,48,16,56,8))

chisq.test(x)

### Exercise 8
nh = c(2,2,3,24,69)
n = sum(nh)

phstar = matrix(NA,5,1)
phstar[1] = dbinom(0,4,0.95)
phstar[2] = dbinom(1,4,0.95)
phstar[3] = dbinom(2,4,0.95)
phstar[4] = dbinom(3,4,0.95)
phstar[5] = dbinom(4,4,0.95)

nhstar = n*phstar

chisq.test(nh,phstar)

### Exercise 9
A = c(92, 114, 82, 164, 167, 110, 135)
B = c(156, 123, 198, 83, 242, 176, 185, 217)

theta_hat = median(c(A,B))

x = matrix(nrow=2,c(sum(A>theta_hat),sum(A<theta_hat),sum(B>theta_hat),sum(B<theta_hat)))

fisher.test(x,alternative="less")

### Exercise 10
x = c(0,1,2,4)
nh = c(42,14,10,4)
lambda = sum(x*nh)/sum(nh)

phstar = dpois(x,lambda=lambda)
phstar[4] = 1 - sum(phstar[1:3])

phstar*sum(nh)

# have to merge last two categories
phstar_new = c(phstar[1:2],sum(phstar[3:4]))
nh_new = c(42,14,14)


chisq.test(nh_new,p=phstar_new)

### Exercise 11
library(RVAideMemoire)

cars <-read.table("Exercise11.txt",header=TRUE,sep='')

x <- data.frame(pref=c(cars[,2],cars[,3],cars[,4]),
             brand=rep(c('C','G','H'),each=12),sub=rep(1:12,3))

cochran.qtest(pref ~brand | sub, data=x)

### Exercise 12
data <- read.table("Exercise12.txt",header=TRUE)
data

kruskal.test(data)

### Exercise 13
x = matrix(nrow=2,c(523,230,345,554))

mcnemar.test(x)

### Exercise 14
nh = c(118,71,62,69)
phstar = c(0.32,0.27,0.16,0.25)
chisq.test(nh,p=phstar)

### Exercise 15
x =c(4,2,3,2,2,2,-2,-1,5,3,4,1,2,-1,5,-1,-2,-2,2)

wilcox.test(x, alternative="t")

### Exercise 16
binom.test(8,10,0.5,alternative = "t")

### Exercise 17
x = c(49,141,2,13,64,8,175,30,7,179,9,7)

data = data.frame(pounds=x,county=rep(c('A','B','C','D'),each=3))

kruskal.test(data$pounds,data$county)

### Exercise 18
cockt<-matrix(c(3,4,2,6,6,4,4,3,2,4,3,5,6,8,5,7,5,8,7,9,
                7,6,6,9,8,10,7,6,4,9),15,2,byrow=TRUE)
cockt.data<-data.frame(A=cockt[,1],B=cockt[,2])

wilcox.test(x=cockt.data$A,y=cockt.data$B,paired=TRUE,alternative="two.sided")

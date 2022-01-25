
### Exercise 8
library(RVAideMemoire)
data <- read.table("probl_8.txt")
data.test <- data.frame (outcome = c(data[,1], data [,2], data [,3], data[,4]), questions = rep
                              (colnames(data), each=15), student= rep(1:15,4))

cochran.qtest(outcome~questions|student,data.test)

### Exercise 9
x = c(24, 26, 32, 40, 22, 25, 30, 44, 21, 18, 26, 27, 28, 28, 27)
x.sel <- x[x!=25]
sign.pos = x.sel-25 > 0

binom.test(sum(sign.pos),length(sign.pos),0.5,alternative = "g")

### Exercise 10
data <- read.table("probl_10.txt",header=TRUE)
D = data[,2]-data[,1]
D.sel <- D[D!=0]

sign.pos = D.sel > 0

binom.test(sum(sign.pos),length(sign.pos),0.5,alternative="g")

### Exercise 11
data <- read.table("probl_11.txt",header=TRUE)
theta_hat = median(data[,1])
data.sel = data[data[,1]!=theta_hat,]
y <- table(data.sel[,1]-theta_hat>0, data[,2])

# homogeneity chisq. test:
chisq.test(y)

### Exercise 12
data <- read.table("probl_12.txt",header=TRUE)

n = 22
K = 3

x <- data.frame(time=c(data$Round_out, data$Narrow_angle, data$Wide_angle),
                    method=rep(c('Round_out', 'Narrow_angle', 'Wide_angle'), each=n),
                    player=rep(data$player,times=K))

friedman.test(time ~ method | player, data = x)

### Exercise 13
data <- data <- read.table("probl_13.csv",header=TRUE,sep =";")

n = 12
K = 3

x <- data.frame(measure=c(data$N, data$R, data$A),
                condition=rep(c('N','R','A'),each=n),
                subject=rep(data$Subject,times=K))

friedman.test(measure ~ condition | subject, data = x)

### Exercise 3
x = matrix(nrow = 2, c(9,1,9,1))
mcnemar.test(x)

x = 9
n = 10
p_0 = 0.5
binom.test(x,n,p_0,alternative="g")

### Exercise 4
x = c(646, 241, 472, 385, 297, 505, 224, 187)
theta_0 = 430

# two-sided:
wilcox.test(x, alternative = "t", mu= theta_0)

# one-sided:
theta_hat = mean(x)
theta_hat
wilcox.test(x, alternative = "l", mu= theta_0)

# for fun:
wilcox.test(x, alternative = "g", mu= theta_0)

### Exercise 5
data <- data <- read.table("exercise5.txt",header=TRUE)
kruskal.test(data[,1],data[,2])

# H_0 rejected, try multiple testing methods:
library(coin)
library(NSM3)

cSDCFlig(alpha=0.01,n=rep(10,4),method="Exact")

cNDWol(alpha=0.01,n=rep(10,4),method="Exact")

### Exercise 6
data <- read.table("exercise6.txt",header=TRUE)
kruskal.test(data[,1],data[,2])

# H_0 rejected, try multiple testing methods:
kruskal_test(quantity ~ line, data=data, distribution = asymptotic)

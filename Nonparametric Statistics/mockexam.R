### Exercise 1
nh = c(20,14,18,17,22,29)
phstar = rep(1/6,6)
phstar*sum(nh)

# all values over 5, so we can use the goodness-of-fit chisquare

chisq.test(nh,p=phstar)

### Exercise 2
data <- matrix(nrow=3,c(11.1,15.3,22.7,23.8,25.6,31.2,25.8,
                        32.6,40.8,52.1,52.8,63.1,59.5,55.3,
                        63.3,65.0,58.8,61.4,41.1,78.1,60.2), byrow=TRUE)

friedman.test(data)

### Exercise 3
library(boot)

data = motor

theta.r <- function(data,i) {
  cor(data[i,1],data[i,2])
}

set.seed(11)
r.boot <- boot(motor,theta.r,R=10000)

sample.r<-cor(data[,1],data[,2])

r.boot.ci <- boot.ci(r.boot,conf=0.95,type=c("norm","basic","perc","bca"))
r.boot.ci

### Exercise 4
library(gss)

n = length(buffalo)
IQR = 1

# Optimal number of bins according to Sturge's rule: 7
k = ceiling(1 + log(n,2))
k

# Optimal number of bins according to Scott's rule: 5
h.sc = 3.5*sd(buffalo)*n^(-1/3)
k.sc = ceiling((max(buffalo)-min(buffalo))/h.sc)
k.sc

# Optimal number of bins according to Freedman-Diaconis' rule: 7
h.fd = 2*IQR(buffalo)*n^(-1/3)
k.fd = ceiling((max(buffalo)-min(buffalo))/h.fd)
k.fd

breaks.fd<- seq(from=min(buffalo),by=h.fd,length.out = ceiling(k.fd)+1)
hist(buffalo,breaks = breaks.fd, freq=FALSE)

plot(density(buffalo,bw="SJ",adjust=1,kernel="g"))
lines(density(buffalo,bw="SJ",adjust=1,kernel="e"))

### Exercise 5
data(motor)
plot(motor$times,motor$accel)

library(np)
motor.cv.ls <- npreg(bwmethod='cv.ls', txdat=motor$times, tydat=motor$accel)
lines(motor$times,predict(motor.cv.ls),col=2)

span.cv=supsmu(x=motor$times, y=motor$accel,span="cv")
lines(span.cv,col=3)
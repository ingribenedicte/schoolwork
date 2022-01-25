library(TSA)

### 1b
n = 1000

# WN process:
e <- rnorm(n)

# theoretical:
ARMAspec(model = list(c(0, 0, 0)), plot = T)
         
# estimated:
spectrum(e)

# AR(1) process:
phi=0.8
e <- rnorm(n + 1)
y = 0
y[1] = 0
for(t in 1:(n - 1)){
  y[t + 1] = phi * y[t] + e[t + 1] 
}

# theoretical:
ARMAspec(model = list(ar = c(phi)), plot = T)

# estimated:
spectrum(y)

### 2c
sigma2 <- 1

vp_func <- function(p, d) {
  if (p * d < 0.5) {
    sigma2*(gamma(1 - 2 * p * d) / (gamma(1 - p * d))^2)
  }
}

# d = 0.499:
vp_1 <- rep(NA, 199)
vp_1[100] <- sigma2
for (p in seq(-0.99, 0.99, 0.01)) {
  vp_1[(p + 1) * 100] = vp_func(p, 0.499)
}

# d = 0.4:
vp_2 <- rep(NA, 199)
vp_2[100] <- sigma2
for (p in seq(-0.99, 0.99, 0.01)) {
  vp_2[(p + 1) * 100] = vp_func(p, 0.4)
}

# d = 0.1:
vp_3 <- rep(NA, 199)
vp_3[100] <- sigma2
for (p in seq(-0.99, 0.99, 0.01)) {
  vp_3[(p + 1) * 100] = vp_func(p, 0.1)
}

plot(vp_1, col = "red", type = 'l', xlab = "p", ylab = "", ylim=c(0,10))
lines(vp_2, col="blue", type='l')
lines(vp_3, col="green", type='l')


### 5
phi_1 <- function(x) {
  1.8*cos(1.5-cos(4*pi*x))
}

phi_2 = -0.81
n = 1000
e <-rnorm(n)

y = 0
y[1]= 0  
y[2] = 0
  
for(t in 3:n){
  y[t] = phi_1(t / n) * y[t - 1] + phi_2 * y[t - 2] + e[t]
}

ts.plot(y, lwd = 1, col="red", main= "ARMA(2,1) process")

# 1)Derive the spectral density function of a zero mean ARMA(2,1) process 
# and plot its graph in R both using the function of omega and using the ARMAspec function
# Check they are equal for the parameter values ?1=0.8, ?2=0.3, ?=0.4, s2=0.2
 
library(TSA)

phi1 <- 0.8
phi2 <- -0.3
theta <- 0.4
sigma2 <- 0.2
omega <- seq(0, pi, pi/500)
f_omega <- sigma2 / (2 * pi) * (1.16 - 0.8 *cos(omega)) / (1.13 - 2.08 * cos(omega) + 1.2*(cos(omega))^2)
par(mfrow = c(1,2))

plot(omega, f_omega, type="l", ylab="f(omega)", lwd=2, col="blue", ylim=c(0,0.08))
points(omega[f_omega==max(f_omega)], max(f_omega))
abline(h=0)
abline(v=0)
abline(h=max(f_omega), lty=2)
abline(v=omega[f_omega==max(f_omega)], lty=2)
legend(x=omega[f_omega==max(f_omega)]-0.3, y=max(f_omega)+0.005, legend=paste("(",round(omega[f_omega==max(f_omega)],3),",",round(max(f_omega),3),")",sep=""),bty="n",cex=0.65)

As<-ARMAspec(model=list(ar=c(phi1,phi2),ma=-theta),plot=F)
plot(As$freq*2*pi,As$spec*sigma2/(2*pi),type="l",ylab="f(omega)",xlab="omega",main="ARMAspec",lwd=2,col="green",ylim=c(0,0.08))
points(omega[f_omega==max(f_omega)],max(f_omega))
abline(h=0)
abline(v=0)
abline(h=max(As$spec*sigma2/(2*pi)),lty=2)
abline(v=(As$freq*2*pi)[As$spec*sigma2/(2*pi)==max(As$spec*sigma2/(2*pi))],lty=2)
legend(x=(As$freq*2*pi)[As$spec*sigma2/(2*pi)==max(As$spec*sigma2/(2*pi))]-0.3, y = max(As$spec*sigma2/(2*pi))+0.005, legend=paste("(",round((As$freq*2*pi)[As$spec*sigma2/(2*pi)==max(As$spec*sigma2/(2*pi))],3),",",round(max(As$spec*sigma2/(2*pi)),3),")",sep=""),bty="n",cex=0.65)

# 2) For an MA(1) process: derive the variance profile 
# and plot (in R) its graph 
# for a range of value of ? in the invertibility region and setting  s2=1. 

sigma2<-1
vp <- rep(NA, 199)
vp[100] <- sigma2
pal = rainbow(5)
for (p in seq(-0.99, 0.99, 0.01)){
	if (p != 0){
		som = 0
		for (k in seq(0, 1000, 1)){	
			add<-((choose(p,k))^2)*(0.9)^(2*k)
			som<-som+add
			}
	vp[(p+1)*100] <- (som)^(1/p)*sigma2
	}
}
plot(seq(-0.99,0.99,0.01), vp, type="l", xlab="p", xlim=c(-1,1.4), ylim=c(0,2), main="MA(1)", col=2, lwd=2)
legend(x=0.95, y=vp[199]+0.15, legend="|Theta|=0.9", bty="n", cex=0.75, text.col=pal[1])
abline(h=0)
abline(v=0)
abline(h=1, lty=2)
c = 2
for (theta in seq(0.7, 0.1, -0.2)){
	vp <- rep(NA,199)
	vp[100] <- sigma2
	for (p in seq(-0.99, 0.99, 0.01)){
		if (p!=0){
			som=0
			for (k in seq(0,1000,1)){	
				add <- ((choose(p,k))^2)*(theta)^(2*k)
				som <- som + add
				}
		vp[(p+1)*100] <- (som)^(1/p)*sigma2
		}
	}
	lines(seq(-0.99,0.99,0.01),vp,col=pal[c],lwd=2)
	legend(x=0.95,y=vp[199]+0.15,legend=paste("|Theta|=",theta),bty="n",cex=0.75,text.col=pal[c])
	c<-c+1
}

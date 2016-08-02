truncbinormrnd2 <- function(m=100,mu1=0,mu2=0,sigma1=1,sigma2=1,rho=0.9,c1=0,c2=0) {
	x <- matrix(0,m,2)
	x1 <- c1
	x2 <- c2
	sd1 <- sqrt(sigma1*(1-rho^2))
	sd2 <- sqrt(sigma2*(1-rho^2))
	for (i in 1:m) {
		mn1 <- mu1+rho*sigma1/sigma2*(x2-mu2)
		a1 <- (c1-mn1)/sd1
		x1 <- mn1+sd1*qnorm(runif(1,pnorm(a1),1))
		mn2 <- mu2+rho*sigma2/sigma1*(x1-mu1)
		a2 <- (c2-mn2)/sd2
		x2 <- mn2+sd2*qnorm(runif(1,pnorm(a2),1))
		x[i,1] <- x1
		x[i,2] <- x2
		}
	return(x)
	}

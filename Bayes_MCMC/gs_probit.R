gs_probit <- function(y,X,beta0,sigma0,iter=1000,burnin=100) {
	dmvnorm0 <- function(mu,sigma) exp(-.5*(log(det(as.matrix(sigma)))+t(mu)%*%solve(sigma)%*%mu))
	n <- nrow(X)
	k <- ncol(X)
	XX <- crossprod(X)
	tilde_sigma <- solve(XX+solve(sigma0))
	sd_beta <- sqrt(diag(tilde_sigma))
	chol_sigma <- chol(tilde_sigma)
	A0beta0 <- solve(sigma0,beta0)
	index_one <- y == 1 
	index_zero <- !index_one
	num_one <- sum(index_one)
	num_zero <- sum(index_zero)
	z <- matrix(0,n,1)
	mcmc <- matrix(0,iter,k)
	colnames(mcmc) <- colnames(X)
	beta_i <- solve(XX,crossprod(X,2*y-1))
	sddr_T <- 0
	sddr_F <- 0
	for (i in -burnin:iter) {
		mean_z <- X%*%beta_i
		prob_z <- pnorm(0,mean_z)
		z[index_one] <- qnorm(runif(num_one,min=prob_z[index_one],max=1),mean_z[index_one])
		z[index_zero] <- qnorm(runif(num_zero,min=0,max=prob_z[index_zero]),mean_z[index_zero])
		tilde_beta <- tilde_sigma%*%(crossprod(X,z)+A0beta0)
		beta_i <- tilde_beta+crossprod(chol_sigma,matrix(rnorm(k),k,1))
		if (i>0) {
			mcmc[i,1:k] <- t(beta_i)		
			sddr_T <- sddr_T+dnorm(0,mean=tilde_beta,sd=sd_beta)
			sddr_F <- sddr_F+dmvnorm0(tilde_beta[2:k],tilde_sigma[2:k,2:k])
		}
	}
	sddr_T <- log10(sddr_T/iter)-log10(dnorm(0,mean=beta0,sd=sqrt(diag(sigma0))))
	sddr_F <- log10(sddr_F/iter)-log10(dmvnorm0(beta0[2:k],sigma0[2:k,2:k]))
	LogBayesFactor <- t(as.matrix(sddr_T))
	Posterior_Mean <- apply(mcmc,2,mean)
	Posterior_SD <- apply(mcmc,2,sd)
	Posterior_QT <- apply(mcmc,2,quantile,c(0.01,0.025,0.05,0.5,0.95,0.975,0.99))
	stats <- rbind(Posterior_Mean,Posterior_SD,Posterior_QT,t(sddr_T))
	rownames(stats)[nrow(stats)] <- "LogBayesFactor"
	results <- list(sample=mcmc,statistics=stats,BF=as.double(sddr_F))
	return(results)
}
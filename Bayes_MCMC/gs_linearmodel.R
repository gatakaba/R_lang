gs_linearmodel <- function(y,X,beta0,sigma0,nu0,lambda0,iter=1000,burnin=100) {
	dmvnorm0 <- function(mu,sigma) exp(-.5*(log(det(as.matrix(sigma)))+t(mu)%*%solve(sigma)%*%mu))
	n <- nrow(X)
	k <- ncol(X)
	Xy <- crossprod(X,y)
	XX <- crossprod(X)
	hat_beta <- solve(XX,Xy)
	rss <- as.double(crossprod(y-X%*%hat_beta))
	tilde_nu <- .5*(nu0+n)
	hat_lambda <- lambda0+rss
	A0 <- solve(sigma0)
	A0beta0 <- solve(sigma0,beta0) 
	mcmc <- matrix(0,iter,k+1)
	colnames(mcmc) <- c(colnames(X),"Sigma")
	beta_i <- hat_beta
	sigma_i <- rss/(n-k)
	sddr_T <- 0
	sddr_F <- 0
	for (i in -burnin:iter) {
		tilde_lambda <- .5*(hat_lambda+as.double(t(beta_i-hat_beta)%*%XX%*%(beta_i-hat_beta)))
		sigma_i <- 1/rgamma(1,tilde_nu,tilde_lambda)
		tilde_sigma <- solve(XX/sigma_i+A0)
		tilde_beta <- tilde_sigma%*%(Xy/sigma_i+A0beta0)
		beta_i <- tilde_beta+crossprod(chol(tilde_sigma),matrix(rnorm(k),k,1))
		if (i>0) {
      ls
			mcmc[i,1:k] <- t(beta_i)
			mcmc[i,k+1] <- sqrt(sigma_i)
			sddr_T <- sddr_T+dnorm(0,mean=tilde_beta,sd=sqrt(diag(tilde_sigma)))
			sddr_F <- sddr_F+dmvnorm0(tilde_beta[2:k],tilde_sigma[2:k,2:k])
		}
	}
	sddr_T <- log10(sddr_T/iter)-log10(dnorm(0,mean=beta0,sd=sqrt(diag(sigma0))))
	sddr_F <- log10(sddr_F/iter)-log10(dmvnorm0(beta0[2:k],sigma0[2:k,2:k]))
	Posterior_Mean <- apply(mcmc,2,mean)
	Posterior_SD <- apply(mcmc,2,sd)
	Posterior_QT <- apply(mcmc,2,quantile,c(0.01,0.025,0.05,0.5,0.95,0.975,0.99))
	stats <- rbind(Posterior_Mean,Posterior_SD,Posterior_QT,LogBayesFactor=c(sddr_T,sddr_F))
	results <- list(sample=mcmc,statistics=stats)
	return(results)
}

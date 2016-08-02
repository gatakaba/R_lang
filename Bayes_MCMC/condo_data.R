data <- read.csv("condo.txt")
n <- nrow(data)
k <- ncol(data)
y <- as.matrix(data[,1])
X0 <- as.matrix(data[,2:k])
X <- cbind(1,X0)
colnames(X)[1] <- "Constant"
# hyper-paraemters
beta0 <- matrix(0,k,1)
sigma0 <- diag(100,k,k)
nu0 <- 0.1
lambda0 <- 0.1

source("gs_linearmodel.R")
result <-gs_linearmodel(y,X,beta0,sigma0,nu0,lambda0)

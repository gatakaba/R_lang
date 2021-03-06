data <- read.csv("credit.txt")
n <- nrow(data)
k <- ncol(data)
y <- as.matrix(data[,1])
X0 <- as.matrix(data[,2:(k-1)])
CA_Dummy <- as.double(data$CA == "California")
X <- cbind(1,X0,CA_Dummy)
colnames(X)[1] <- "Constant"
# hyper-paraemters
beta0 <- matrix(0,k,1)
sigma0 <- diag(100,k,k)
nu0 <- 0.1
lambda0 <- 0.1




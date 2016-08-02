library(MCMCpack)


n_iter=1100

n = 50.0
x_mean = 5.0
mu0 = 0.0
tau0_square = 100.0
n0 = 10.0
lambda0 = 8.0
mu = 0.0
sigma = 1.0

mu_list=numeric(n_iter)
sigma_list=numeric(n_iter)



for (i in 1:n_iter){
  p1 = (n * sigma ** -1 * x_mean + tau0_square ** -1 * mu0) /(n * sigma ** -1 + tau0_square ** -1)
  p2 = 1 / (n * sigma ** -1 + tau0_square ** -1)
  p3 = (n + n0) / 2.0
  p4 = (n * (mu - x_mean) ** 2 + 50 + lambda0) / 2.0 
  
  mu = rnorm(1,p1, p2)
  sigma = rinvgamma(1,shape=p3,scale=p4)
  
  mu_list[i]=mu
  sigma_list[i]=sigma
  
}
png("mu.png")
plot(mu_list, type = "l")
dev.off()
png("sigma.png")
plot(sigma_list, type = "l")
dev.off()



mu_list=mu_list[100:1100]
sigma_list=sigma_list[100:1100]


png("mu_hist.png")
hist(mu_list, type = "l")
dev.off()
png("sigma_hist.png")
hist(sigma_list, type = "l")
dev.off()

print(mean(mu_list))
print(mean(sigma_list))

print(mean(mu_list)+qt(0.975,n_iter-1)*sd(mu_list)/sqrt(length(mu_list)))
print(mean(mu_list)-qt(0.975,n_iter-1)*sd(mu_list)/sqrt(length(mu_list)))

print(mean(sigma_list)+qt(0.975,n_iter-1)*sd(sigma_list)/sqrt(length(sigma_list)))
print(mean(sigma_list)-qt(0.975,n_iter-1)*sd(sigma_list)/sqrt(length(sigma_list)))



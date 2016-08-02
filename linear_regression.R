#データの作成
N <- 100
weight<- 60+ rnorm(N)*3
height <-weight *3+rnorm(N)*10-10
#height=a*weight+const というモデルにフィティング
lm <- lm(height~weight)
print(summary(lm))
#結果の表示

plot(weight,height)
abline(lm , col="blue")

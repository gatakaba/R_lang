library(ggplot2)

N=1000
point=rnorm(N)
english.exam = data.frame(
  Reading = 2*point,
  Writing = 2*point+rnorm(N),
  Speaking= point+rnorm(N)*5
)
g = ggplot(english.exam, aes(Reading, Writing))
#g + geom_point()
g + geom_point() + stat_smooth(method = "lm")

library(ggplot2)
library(reshape2)

set.seed(1)
n_values = c(100, 250, 1000)
posterior = matrix(NA, 500, 3)
for(i in 1:length(n_values)){
  data = generate_from_simple_logistic_model(n_values[i])
  obj = gibbs_sampler(data$y, data$X, lambda=100, n_iter=500)
  posterior[, i] = obj$beta[,1]
}
df = data.frame(posterior[-c(1:100), ])
colnames(df) = n_values
df.m = melt(df)
ggplot(df.m, aes(value, col=variable)) + geom_density() + theme_bw() + 
  scale_color_discrete("Sample size") + xlab(expression(beta[1])) + ylab("Posterior distribution")

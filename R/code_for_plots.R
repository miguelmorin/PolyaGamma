# Prior, likelihood and posterior plot

if(FALSE){
  library(mvtnorm)
  library(ggplot2)
  library(gridExtra)

  n = 100
  x1 = rnorm(n)
  x2 = rnorm(n)
  psi = 1*x1 + 1*x2
  y = rbinom(n, 1, 1 / (1 + exp(-psi)))
  plot(x1, x2, col=y+1, pch=16)


  Xnorm1 = rmvnorm(20000, mean=c(0, 0), sigma = 0.2*diag(2))
  dfnorm1 = data.frame(Xnorm1)
  Xnorm2 = rmvnorm(20000, mean=c(0, 0), sigma = 0.03*diag(2))
  dfnorm2 = data.frame(Xnorm2)

  X = cbind(x1, x2)
  obj = gibbs_sampler(y, X, lambda=0, n_iter_total=1050, burn_in=50)
  df1 = data.frame(obj$beta)
  obj = gibbs_sampler(y, X, lambda=1, n_iter_total=1050, burn_in=50)
  df2 = data.frame(obj$beta)
  obj = gibbs_sampler(y, X, lambda=10, n_iter_total=1050, burn_in=50)
  df3 = data.frame(obj$beta)


  xlim = c(-1, 2)
  ylim = c(-1, 2)

  p1 = ggplot() +
    stat_density2d(aes(X1, X2), dfnorm1, col="black") +
    stat_density2d(aes(X1, X2), df1) +
    stat_density2d(aes(X1, X2), df2, col="red") +
    theme_classic() +
    geom_vline(xintercept=0) + geom_hline(yintercept=0) +
    coord_cartesian(xlim=xlim, ylim=ylim) +
    ggtitle("(a)") + xlab(expression(beta[1])) + ylab(expression(beta[2]))

  p2 = ggplot() +
    stat_density2d(aes(X1, X2), dfnorm2, col="black") +
    stat_density2d(aes(X1, X2), df1) +
    stat_density2d(aes(X1, X2), df3, col="red") +
    theme_classic() +
    geom_vline(xintercept=0) + geom_hline(yintercept=0) +
    coord_cartesian(xlim=xlim, ylim=ylim) +
    ggtitle("(b)") + xlab(expression(beta[1])) + ylab(expression(beta[2]))

  grid.arrange(p1, p2, ncol=2)
}

logistic = function(x)1/(1+exp(-x))

#' Simulate data from a simple logistic model
#' @export
generate_from_simple_logistic_model = function(n = 1000){
  x = rnorm(n)
  X = cbind(1, x)
  y = rbinom(n, 1, logistic(1 + 1*x))
  return(list(X=X, y=y))
}

# data = generate_from_simple_logistic_model(1000)
# obj = gibbs_sampler(data$y, data$X, b=rep(0, 2), B=10*diag(2), 500)
# plot(obj)

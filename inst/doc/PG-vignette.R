## ----fig.height=2--------------------------------------------------------
library(PolyaGamma)
data = generate_from_simple_logistic_model(n=100)
obj = gibbs_sampler(data$y, data$X, lambda=0.001, n_iter_total=150, burn_in=50)
obj


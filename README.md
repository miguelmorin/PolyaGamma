# PolyaGamma

An R package for Bayesian logistic regression. 
The posterior distribution of the parameters is obtained via Gibbs sampling using Polya-Gamma latent variables (see paper ["Bayesian Inference for Logistic Models Using Pólya–Gamma Latent Variables"](http://www.tandfonline.com/doi/abs/10.1080/01621459.2013.829001) for details). 

This package can be installed as follows:

```R
devtools::install_github("kasparmartens/PolyaGamma")
```

A small example, how to use the main function `gibbs_sampler`

```R
# load library
library(PolyaGamma)
# generate data from the logistic model 1 / (1 + exp(-(1 + x)))
data = generate_from_simple_logistic_model(n=100)
# obtain a sample from the posterior distribution of beta
obj = gibbs_sampler(data$y, data$X, lambda=0.001, n_iter_total=200, burn_in=50)
obj
plot(obj)

```

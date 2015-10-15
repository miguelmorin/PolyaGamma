# Here will be some functions for Polya-Gamma Gibbs sampling

# Naive sampler for PG(1, z)
# (based on the finite approximation of infinite sum)
single_sample_naive_polyagamma = function(z, max_k = 100){
  g = rexp(max_k, 1)
  out = 1 / (2*pi**2) * sum(g / (0.25*(1:max_k - 1)**2 + z**2 / (4*pi**2)))
  return(out)
}#end function

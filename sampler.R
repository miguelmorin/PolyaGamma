library(mvtnorm)

# Here will be some functions for Polya-Gamma Gibbs sampling

# Naive sampler for PG(1, z)
# (based on the finite approximation of infinite sum)
rpolyagamma_naive = function(z, max_k = 100){
  g = rexp(max_k, 1)
  out = 1 / (2*pi**2) * sum(g / ((1:max_k - 1/2)**2 + z**2 / (4*pi**2)))
  return(out)
}#end function


#FUNCTION: cdf(x) of inverse Normal distribution with parameters mu and lambda
#AUTHOR: SHERMAN IP
#DATE: 15/10/15
pinversen = function(x,mu,lambda){
  #check if x, mu and lambda are positive real numbers
  if ((x<=0)|(mu<=0)|(lambda<=0)){
    stop("Parameters in pinversen() are not of the correct type");
  }#end if
  
  #assign variables
  sqrt_lambda_over_x = sqrt(lambda/x);
  x_over_mu = x/mu;
  
  #work out the cdf and return it
  cdf = pnorm(sqrt_lambda_over_x*(x_over_mu-1));
  cdf = cdf + exp(2*lambda/mu)*pnorm(-sqrt_lambda_over_x*(x_over_mu+1));
  return (cdf);
}#end pinversen


#FUNCTION: Piecewise coefficients
#PARAMETERS:
#Parameter x, nth term, truncate at t
#AUTHOR: SHERMAN IP
#DATE: 15/10/15
a_n = function(x,n,t){
  #check if n is positive integer, x is positive real, t is positive real
  if ((x<=0)|(t<=0)|(round(n)!=n)|(n<0)){
    stop("Parameters in a_n are not of the correct type");
  }#end if
  
  #set a for x<=t
  if (x<=t){
    a = pi*(n+0.5)*(sqrt(2/(pi*x)))^3*exp(-2*(n+0.5)^2/x);
  }#end if
  #else, set a for x>t
  else{
    a = pi*(n+0.5)*exp(-0.5*(n+0.5)^2*pi^2*x);
  }#end else
  
  #return a
  return(a);
}#end a_n

#FUNCTION: Generate sample from inverse Normal with parameters (mu, 1) truncated with a max of t
#AUTHOR: SHERMAN IP
#DATE: 15/10/15
rinversen = function(mu,t){
  #check if mu and t are positive real
  if ((mu<=0)|(t<=0)){
    stop("Parameters in rinversen are not of the correct type");
  }#end if
  
  #x is the random sample
  x = t;
  
  #while x is equal or more than t, sample
  while(x>=t){
    
    #for large mu, use chi-squared approximation
    if (mu>t){
      repeat{
        #sample exponential (full details in paper)
        repeat{
          E = rexp(2);
          if (E[1]^2<=2*E[2]/t){
            E = E[1];
            break;
          }#end if
        }#end repeat
        
        #assign x, alpha
        x = t/(1+t*E)^2;
        alpha = exp(-0.5*x/mu^2);
        
        #accept if U(0,1) <= alpha
        if (runif(1)<=alpha){
          break;
        }#end if
        
      }#end repeat
    }#end if
    
    #else for small mu...
    else{
      y = (rnorm(1))^2;
      x = mu+0.5*mu^2*y-0.5*mu*sqrt(4*mu*y+(mu*y)^2);
      if (runif(1)>mu/(mu+x)){
        x = mu*mu/x;
      }#end if
    }#else
    
  }#end while
  
  #return the sample
  return(x);
}#end rinverse

#DEBUG FUNCTION: pdf(x) of inverse normal with parameter mu, lambda
#AUTHOR: SHERMAN IP
#DATE: 15/10/15
dinversen = function(x,mu,lambda){
  return(sqrt(lambda/(2*pi*x^3))*exp(-lambda*(x-mu)^2/(2*mu^2*x)));
}#end dinversen

#DEBUG FUNCTION: pdf(x|x<t) where x~IG(mu,lambda)
dtruncatedinversen = function(x,mu,lambda,t){
  return(dinversen(x,mu,lambda)/pinversen(t,mu,lambda));
}#end dtruncatedinversen

#DEBUG FUNCTION: plot histogram and pdf of truncated IG(mu,1) at t
testinversen = function(mu,t){
  x_plot = seq(from=0.001,to=t,by=0.001);
  f_plot = sapply(x_plot,dtruncatedinversen,mu=mu,lambda=1,t=t);
  x = replicate(10000,rinversen(mu,t));
  hist(x,freq=FALSE);
  lines(x_plot,f_plot);
}#end testinversen


# Generate parameter vector beta from a multivariate normal
# beta ~ N(m, V), see details in the paper
generate_mv_normal = function(w, y, X, b, B){
  Binv = solve(B)
  temp = t(X) %*% diag(w) %*% X + Binv
  V = solve(temp)
  kappa = y - 0.5
  m = V %*% (t(X) %*% kappa + Binv %*% b)
  beta = as.numeric(rmvnorm(1, mean=m, sigma=V))
  return(beta)
}

# Check whether the user has provided arguments 
# with matching dimensionalities
check_dimensions = function(y, X, b, B){
  if(!is.vector(y)) stop("y must be a vector")
  if(!is.vector(b)) stop("b must be a vector")
  if(!is.matrix(X)) stop("X must be a matrix")
  if(!is.matrix(B)) stop("B must be a matrix")
  if(length(y) != nrow(X)) stop("nrow(X) must equal length(y)")
  if(length(b) != ncol(X)) stop("ncol(X) must equal length(b)")
}

# Gibbs two-step sampling procedure 
# for parameter vector beta and the latent Polya-Gamma variables
gibbs_sampler = function(y, X, b, B, n_iter=100){
  # Check if everything is OK with dimensions
  check_dimensions(y, X, b, B)
  
  # number of parameters
  m = ncol(X)
  # number of data points
  n = nrow(X)
  # Starting values for beta; initialise w
  beta = b
  w = rep(NA, n)
  
  # Store the values of all betas and all w
  beta_all = matrix(0, n_iter, m)
  w_all = matrix(0, n_iter, n)
  
  for(k in 1:n_iter){
    # draw elements of w from PG
    for(i in 1:n){
      psi = as.numeric(X[i, ] %*% beta)
      w[i] = rpolyagamma_naive(psi)
    }
    # draw beta from a multivariate normal
    beta = generate_mv_normal(w, y, X, b, B)
    beta_all[k, ] = beta
    w_all[k, ] = w
  }
  return(list("beta" = beta_all, "w" = w_all))
}

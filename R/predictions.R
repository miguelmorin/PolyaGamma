#FUNCTION: return the logistic regression for multiple betas
#PARAMETERS: X and beta_design are design matrices
logisticRegression = function(X,beta_design){
  #check if X and beta_design are matrices
  if ( (!(is.matrix(X)|is.vector(X))) | (!(is.matrix(beta_design)|is.vector(beta_design))) ){
    stop("Parameters in logisticRegression are not of the correct type");
  }#end if
  #get a matrix of etas
  eta = beta_design %*% t(X);
  #return logistic regression
  p = 1/(1+exp(-eta));
  return(p);
}#end logisticRegression

#' Predict the y value for all data points,
#' given  and the sample from posterior of betas
#'
#' @param X A design matrix.
#' @param beta Sample from the posterior distribution of betas. Alternatively, could be a vector with estimated betas.
#'
#' @export
get_predictions = function(X, beta){
  samples = logisticRegression(X, beta)
  posterior_mean = colMeans(samples)
  ypred = round(posterior_mean)
  return(ypred)
}

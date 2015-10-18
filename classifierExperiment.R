#AUTHOR: SHERMAN IP
#DATE: 16/10/15

classifierExperiment = function(){
  #set random seed
  set.seed(92486);
  
  #read data
  data = read.table("SAheart.data",sep=",",head=T,row.names=1);
  #get the response vector
  y = data$chd;
  #select numeric features and save it as a matrix
  data = data[c("sbp","tobacco","ldl","famhist","alcohol","age")];
  data$famhist = (data$famhist=="Present");
  X = as.matrix(data);
  
  #get the sample size and the number of features
  n_total = length(y);
  p = ncol(X)+1;
  
  #get size of training and testing set
  n_train = round(0.5*n_total);
  n_test = n_total - n_train;
  
  #assign data to training/testing set
  data_pointer = sample(1:n_total,n_total);
  train_pointer = data_pointer[1:n_train];
  test_pointer = data_pointer[-(1:n_train)];
  
  #assign variables for training/testing set
  X_train = X[train_pointer,];
  y_train = y[train_pointer];
  X_test = X[test_pointer,];
  y_test = y[test_pointer];
  
  #normalise the data
  x_center = colMeans(X_train);
  x_scale = apply(X_train,2,sd);
  X_train = scale(X_train,x_center,x_scale);
  X_test = scale(X_test,x_center,x_scale);
  
  #fit using logistic model
  model = glm(y_train~X_train,family=binomial(link="logit"));
  beta_logistic = model$coefficients;
  
  #add constant term to X_train and X_test
  X_train = cbind(1,X_train);
  X_test = cbind(1,X_test);
  
  #get error of logistic regression
  logistic_train_error = getTestError(y_train,as.vector(logisticRegression(X_train,beta_logistic)));
  logistic_test_error = getTestError(y_test,as.vector(logisticRegression(X_test,beta_logistic)));
  
  #get range of lambda to investigate
  lambda_exp_vector = seq(3,10);
  
  n_error = 20; #number of times to repeat the experiment
  n_samples = 100; #number of betas to sample from the chain
  n_chain = 150; #the length of the chain
  
  #create matrix to stroe training and testing error
  test_error = matrix(0,ncol=length(lambda_exp_vector),nrow=n_error);
  train_error = test_error;
  
  #for every lambda
  for (i in 1:length(lambda_exp_vector)){
    #get lambda
    lambda = lambda_exp_vector[i];
    #repeat n_error times
    for (j in 1:n_error){
      #get a chain
      chain = gibbs_sampler(y_train, X_train, b = replicate(p,0), B = diag(replicate(p,lambda)),n_chain)$beta;
      #take the last part of the chain
      beta_posterior = chain[(n_chain-n_samples+1):n_chain,];
      #average the logistic regression, round it and use it for prediction
      train_error[j,i] = getTestError(y_train,(colMeans(logisticRegression(X_train,beta_posterior))));
      test_error[j,i] = getTestError(y_test,(colMeans(logisticRegression(X_test,beta_posterior))));
    }#end for
  }#end for
  
  #plot the training and testing error
  par(mfrow=c(1,2));
  boxplot(train_error,names=paste("10E",sapply(lambda_exp_vector,toString),sep=""),xlab="Prior variance",ylab="Training error");
  abline(logistic_train_error,0,col="red")
  boxplot(test_error,names=paste("10E",sapply(lambda_exp_vector,toString),sep=""),xlab="Prior variance",ylab="Testing error");
  abline(logistic_test_error,0,col="red");
}#end classifierExperiment

#FUNCTION: get test error of classifing y using prediction yHat
getTestError = function(y,yHat){
  
  #check y is a binary vector, yHat is a vector between 0 and 1
  if ( (!is.vector(y))|(!is.vector(yHat))| (!all((y==0)|(y==1))) | (!all((yHat>=0)|(yHat<=1))) ){
    stop("Parameters in getTestError are not of the correct type");
  }#end if
  
  #check y and yHat have the correct dimenstions
  n = length(y);
  if ( (n!=length(yHat)) ){
    stop("Dimensions are not of the correct size in getTestError");
  }#end if
  
  #return error
  return(sum(y!=round(yHat))/n);
  
}#end getTestError


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
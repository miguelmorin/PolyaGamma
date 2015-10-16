#AUTHOR: SHERMAN IP
#DATE: 16/10/15

classifierExperiment = function(){
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
  n_train = round(0.75*n_total);
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
  
  #fit using logistic model
  model = glm(y_train~X_train,family=binomial(link="logit"));
  beta_logistic = model$coefficients;
  
  #add constant term to X_train and X_test
  X_train = cbind(1,X_train);
  X_test = cbind(1,X_test);
  
  #get error of logistic regression
  logistic_train_error = getTestError(X_train,y_train,beta_logistic);
  logistic_test_error = getTestError(X_test,y_test,beta_logistic);
  
  #get range of lambda to investigate
  lambda_exp_vector = seq(1,5);
  
  n_error = 20; #number of chains to run
  n_samples = 50; #number of betas to sample from the chain
  n_chain = 100; #the length of the chain
  
  #create matrix to stroe training and testing error
  test_error = matrix(0,nrow=length(lambda_exp_vector),ncol=n_error);
  train_error = test_error;
  
  #for every lambda
  for (i in 1:length(lambda_exp_vector)){
    #get lambda
    lambda = lambda_exp_vector[i];
    
    #for n_error times
    for (j in 1:n_error){
      #get a chain
      chain = gibbs_sampler(y_train, X_train, b = replicate(p,0), B = diag(replicate(p,lambda)),n_chain)$beta;
      #take the last part of the chain and take the sample mean
      beta_posterior = colMeans(chain[(n_chain-n_samples):n_chain,]);
      #get the traing and testing error
      train_error[i,j] = getTestError(X_train,y_train,beta_posterior);
      test_error[i,j] = getTestError(X_test,y_test,beta_posterior);
    }#end for
  }#end for
  
  #plot the training and testing error
  par(mfrow=c(1,2));
  boxplot(t(train_error),names=paste("10E",sapply(lambda_exp_vector,toString),sep=""),xlab="Prior variance",ylab="Training error");
  abline(logistic_train_error,0,col="red")
  boxplot(t(test_error),names=paste("10E",sapply(lambda_exp_vector,toString),sep=""),xlab="Prior variance",ylab="Testing error");
  abline(logistic_test_error,0,col="red");
}#end classifierExperiment

#FUNCTION: get test error of classifing y with design matrix X using the logistic regression with parameter beta
getTestError = function(X,y,beta){
  
  #check X is a vector, y is a binary vector, beta is a vector
  if ((!is.matrix(X))|(!is.vector(y))|(!is.vector(beta))|!all((y==0)|(y==1))){
    stop("Parameters in getTestError are not of the correct type");
  }#end if
  
  #check X, y and beta have the correct dimenstions
  n = length(y);
  p = length(beta);
  if ( (n!=nrow(X)) | (p!=ncol(X)) ){
    stop("Dimensions are not of the correct size in getTestError");
  }#end if
  
  #get vector of symanatic components
  eta = X %*% beta;
  #get vector of logistic regression
  p = 1/(1+exp(-eta));
  #return error
  return(sum(y!=round(p))/n);
  
}#end getTestError
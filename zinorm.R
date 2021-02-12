
## Predictions is a vector of predicted values
zi_norm_estimator <- function(predictions,tol=0.00001){
  p0_init <- mean(predictions==0) # Estimate initial value for the zero-iflation factor
  mu_norm_init <- mean(predictions) # Estimate initial value for mean of the normal distribution
  sd_norm_init <- sd(predictions) # Estimate initial value for sd of normal distribution
  ll_old <- sum(dnorm(predictions,mu_norm_init,sd_norm_init,log=T)) # Estimate loglikelihood of normal with no zero inflation (initial)
  
  p0 <- rep(p0_init/(p0_init + (1 - p0_init) * dnorm(0, mu_norm_init)),length(predictions)) # Estimate probability of zero given the likelihood of zero from original normal distribution
  p0[predictions>0] <- 0 # Values greater than zero have probability=0 of being 0
  weights <-  (1 - p0)/mean(1-p0) # Create weights based on probability of being zero, normalize weights
  mu_norm <- mean(weights*predictions) # Estimate new mean of normal filtering out excess zeroes
  sd_norm <- sqrt(sum(weights)/(sum(weights)-1)*mean(weights*(predictions-mu_norm)^2)) # Estimate new sd of normal filtering out excess zeroes
  ll_new <-  sum(weights*dnorm(predictions,mu_norm,sd_norm,log=T))/mean(weights) # Estimate new loglikelihood of normal filtering out excess zeroes
  
  ## Repeat the block above until convergence
  while (abs((ll_old - ll_new)/ll_old) > tol) {
    ll_old <- ll_new
    p0 <- p0/(p0 + (1 - p0) * dnorm(0, mu_norm))
    p0[predictions>0] <- 0
    weights <-  (1 - p0)/mean(1-p0)
    mu_norm <- mean(weights*predictions)
    sd_norm <- sqrt(sum(weights)/(sum(weights)-1)*mean(weights*(predictions-mu_norm)^2))
    ll_new <-  sum(weights*dnorm(predictions,mu_norm,sd_norm,log=T))/mean(weights)
  }
  
  ## Cleaning the results and returning
  out <- c(mean(p0),mu_norm,sd_norm)
  names(out) <- c('ZI-factor','mu-norm','sd-norm')
  return(out)
}











# .libPaths(c(.libPaths(), getwd()))

# Likelihood Function
l <- function(x, z, y, bet, sig, ga, phi){
  
  x <- as.matrix(x)
  z <- as.matrix(x)
  y <- as.numeric(y)
  bet <- as.numeric(bet)
  sig <- as.numeric(sig)
  ga  <- as.numeric(ga)
  
  eta <- x %*% bet + z %*% ga
  mu  <- exp(eta)/(1+exp(eta))
  mu  <- pmin(mu, 1-.Machine$double.eps)
  mu  <- pmax(mu, .Machine$double.eps )
  mu[is.nan(mu)] <- 1-.Machine$double.eps
  
  return(
    mu*phi*log(y/(1-y)) + phi*log(1-y) + 
      log(gamma(phi)) - log(gamma(mu*phi)) - log(gamma((1-mu)*phi))-
      sum(0.5*ga^2*sig^(-1)) - 0.5*log(prod(sig)) 
  )
  
}

# Link Derivative
d_mu   <- function(eta){
  
  D_mu <- pmin(exp(eta), .Machine$double.xmax)/(1+pmin(exp(eta), .Machine$double.xmax))^2
  D_mu <- pmin(D_mu, 1-.Machine$double.eps)
  D_mu <- pmax(D_mu, .Machine$double.eps )
  D_mu[is.nan(D_mu)]  <- .Machine$double.eps
  
  
  return(D_mu)
}
d_mu_2 <- function(eta){
  
  D_mu_2 <- exp(eta)*(1-exp(eta))/(1+exp(eta))^3
  D_mu_2 <- pmin(D_mu_2, 1-.Machine$double.eps)
  D_mu_2 <- pmax(D_mu_2, .Machine$double.eps )
  D_mu_2[is.nan(D_mu_2)]  <- .Machine$double.eps
  
  return(D_mu_2)
  
}

# First order gradient of gamma
g_gamma   <- function(x, z, y, bet, sig, ga, phi){
  
  x <- as.matrix(x)
  z <- as.matrix(z)
  y <- as.numeric(y)
  bet <- as.numeric(bet)
  sig <- as.numeric(sig)
  ga <- as.numeric(ga)
  
  eta <- x %*% bet + z %*% ga
  mu  <- pmin(exp(eta), .Machine$double.xmax)/(1+pmin(exp(eta), .Machine$double.xmax))
  mu[mu <= .Machine$double.eps]   <- .Machine$double.eps
  mu[mu >= 1-.Machine$double.eps] <- 1-.Machine$double.eps
  
  g <- t(t(phi*log(y/(1-y))*d_mu(eta) - 
                 phi*digamma(mu*phi)*d_mu(eta) + 
                 phi*digamma((1-mu)*phi)*d_mu(eta)) %*% z) - diag(1/sig)%*%ga 
  
  return(g)
  
}
# Second order gradient of gamma
g_gamma_2 <- function(x, z, y, bet, sig, ga, phi){
  
  x <- as.matrix(x)
  z <- as.matrix(z)
  y <- as.numeric(y)
  bet <- as.numeric(bet)
  sig <- as.numeric(sig)
  ga  <- as.numeric(ga)
  
  eta <- x %*% bet + z %*% ga
  mu  <- pmin(exp(eta), .Machine$double.xmax)/(1+pmin(exp(eta), .Machine$double.xmax))
  mu[mu <= .Machine$double.eps]   <- .Machine$double.eps
  mu[mu >= 1-.Machine$double.eps] <- 1-.Machine$double.eps
  
  V   <- (phi*log(y/(1-y))-phi*digamma(mu*phi)+phi*digamma((1-mu)*phi))*d_mu_2(eta) -
    phi^2*(trigamma(mu*phi)+trigamma((1-mu)*phi))*(d_mu(eta))^2
  ZVZ <- t(z * as.numeric(V)) %*% z
  
  g <- ZVZ - diag(1/sig)
  
  return(g)
  
}
# First order gradient of beta
g_beta    <- function(x, z, y, bet, sig, ga, phi){
  
  x <- as.matrix(x)
  z <- as.matrix(z)
  y <- as.numeric(y)
  sig <- as.numeric(sig)
  ga  <- as.numeric(ga)
  bet <- as.numeric(bet)
  
  eta <- x %*% bet + z %*% ga
  mu   <- exp(eta)/(1+exp(eta))
  mu   <- pmin(mu, 1-.Machine$double.eps)
  mu   <- pmax(mu, .Machine$double.eps )
  mu[is.nan(mu)] <- 1-.Machine$double.eps
  
  G_bet <- apply(x * phi*as.numeric((log(y/(1-y)) - digamma(mu*phi) + digamma((1-mu)*phi))*d_mu(eta)),2,sum)
  
  return(G_bet)
  
}
# First order gradient of sigma
g_sigma   <- function(ga, sig){
  
  sig   <- as.numeric(sig)
  ga    <- as.numeric(ga)
  
  return(0.5*ga^2/sig^2-0.5/sig)
  
}
# First order gradient of phi
g_phi     <- function(x, z, y, bet, sig, ga, phi){ 
  
  x <- as.matrix(x)
  z <- as.matrix(z)
  y <- as.numeric(y)
  sig <- as.numeric(sig)
  ga  <- as.numeric(ga)
  bet <- as.numeric(bet)
  
  eta  <- x %*% bet + z %*% ga
  mu   <- exp(eta)/(1+exp(eta))
  mu   <- pmin(mu, 1-.Machine$double.eps)
  mu   <- pmax(mu, .Machine$double.eps )
  mu[is.nan(mu)] <- 1-.Machine$double.eps
  
  G_phi <- sum(mu*log(y/(1-y)) + log(1-y) + digamma(phi) - mu*digamma(mu*phi) - (1-mu)*digamma((1-mu)*phi))
  
  return(G_phi)
  
}  



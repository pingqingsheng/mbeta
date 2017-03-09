# -----------------------------------------------------------------------------------------
# define some function to make the SGD neat

epsilon <- 1e-4

g_h <- function(x, y, ga, bet){
  
  # x <- X
  # y <- Y
  # ga  <- gamma_temp[i,]
  # bet <- beta_old_sgd
  
  x <- as.matrix(x)
  y <- as.matrix(y) ; y[y==1] <- 1-.Machine$double.eps; y[y==0] <- .Machine$double.eps
  ga  <- as.matrix(ga)
  bet <- as.matrix(bet)
  
  eta <- x[,1:5] %*% bet + x[,6:10] %*% ga
  mu  <- (exp(eta)^(-1)+1)^(-1) 
  mu[mu > 1-.Machine$double.eps] <- 1-.Machine$double.eps
  mu[mu < .Machine$double.eps]   <- .Machine$double.eps
  
  numer  <- exp(eta)
  denom  <- (1+exp(eta))^2
  d_mu   <- as.vector(numer/denom)
  d_mu[mu < .Machine$double.eps   | mu == .Machine$double.eps]   <- .Machine$double.eps
  d_mu[mu > 1-.Machine$double.eps | mu == 1-.Machine$double.eps] <- .Machine$double.eps
  
  g_beta <- phi_old_sgd*t(x[,1:5]) %*% diag(d_mu) %*% (log(y/(1-y)) - digamma(mu*phi_old_sgd) + digamma((1-mu)*phi_old_sgd))
  g_phi  <- sum(mu*(log(y/(1-y)) - digamma(mu*phi_old_sgd) + digamma((1-mu)*phi_old_sgd)) + log(1-y) - digamma((1-mu)*phi_old_sgd) + digamma(phi_old_sgd))
  g_tau  <- -1/2 * outer(as.numeric(ga), as.numeric(ga))
  
  return(list(g_beta = g_beta, g_phi = g_phi, g_tau = g_tau))
} 

# before calculating that cumbersome matrix derivative wrt to vector, define h_2 in advance

# input x and y should be vector within a cluster

d_h_1 <- function(x, y, ga, bet, tau){
  
  x <- as.matrix(x)
  y <- as.matrix(y) ; y[y==1] <- 1-.Machine$double.eps; y[y==0] <- .Machine$double.eps
  ga  <- as.matrix(ga)
  bet <- as.matrix(bet)
  
  eta <- x[,1:5] %*% bet + x[,6:10] %*% ga
  mu  <- (exp(eta)^(-1)+1)^(-1)
  mu[mu > 1-.Machine$double.eps] <- 1-.Machine$double.eps
  mu[mu < .Machine$double.eps]   <- .Machine$double.eps
  
  numer  <- exp(eta)
  denom  <- (1+exp(eta))^2
  d_mu   <- as.matrix(numer/denom)
  d_mu[mu < .Machine$double.eps   | mu == .Machine$double.eps]   <- .Machine$double.eps
  d_mu[mu > 1-.Machine$double.eps | mu == 1-.Machine$double.eps] <- .Machine$double.eps
  
  d_h_1 <- phi_old_sgd * t(x[,6:10]) %*% diag(as.vector(d_mu)) %*% (log(y/(1-y)) - digamma(mu*phi_old_sgd) + digamma((1-mu)*phi_old_sgd)) - tau %*% ga
  
  return(d_h_1)
}

d_h_2 <- function(x, y, ga, bet, tau){
  
  # x <- X
  # y <- Y
  # ga  <- gamma_temp[i,]
  # bet <- beta_old_sgd
  # tau <- tau_old_sgd
  
  x <- as.matrix(x)
  y <- as.matrix(y) ; y[y==1] <- 1-.Machine$double.eps; y[y==0] <- .Machine$double.eps
  ga  <- as.matrix(ga)
  bet <- as.matrix(bet)
  
  eta <- as.vector(x[,c(1:5)] %*% bet + x[,c(6:10)] %*% ga)
  mu  <- (exp(eta)^(-1)+1)^(-1) 
  mu[mu > 1-.Machine$double.eps] <- 1-.Machine$double.eps
  mu[mu < .Machine$double.eps]   <- .Machine$double.eps
  
  numer  <- exp(eta)
  denom  <- (1+exp(eta))^2
  d_mu   <- as.numeric(numer/denom)
  d_mu[mu < .Machine$double.eps   | mu == .Machine$double.eps]   <- .Machine$double.eps
  d_mu[mu > 1-.Machine$double.eps | mu == 1-.Machine$double.eps] <- .Machine$double.eps
  
  numer  <- exp(eta)*(1-exp(eta))
  denom  <- (1+exp(eta))^3
  d_2_mu <- as.numeric(numer/denom)
  d_2_mu[mu < .Machine$double.eps   | mu == .Machine$double.eps]   <- .Machine$double.eps
  d_2_mu[mu > 1-.Machine$double.eps | mu == 1-.Machine$double.eps] <- .Machine$double.eps
  
  A <- as.numeric((log(y/(1-y)) - digamma(mu*phi_old_sgd) + digamma((1-mu)*phi_old_sgd))*d_2_mu)
  B <- -as.numeric(phi_old_sgd*(trigamma(mu*phi_old_sgd)+trigamma((1-mu)*phi_old_sgd))*d_mu^2)
  w <- A+B
  d_h_2 <- phi_old_sgd * t(x[,6:10]) %*% diag(w) %*% x[,6:10] - tau 
  
  return(d_h_2)
}

# now let's deal with the trace of a rank 3 tensor 

g_d_h_2 <- function(x, y, ga, bet, tau){
  
  # x <- X
  # y <- Y
  # ga  <- gamma_temp[i,]
  # bet <- beta_old_sgd
  # tau <- tau_old_sgd

  x <- as.matrix(x)
  y <- as.matrix(y); y[y==1] <- 1-.Machine$double.eps ; y[y==0] <- .Machine$double.eps
  ga  <- as.matrix(ga)
  bet <- as.matrix(bet)
  
  eta <- as.vector(x[,c(1:5)] %*% bet + x[,c(6:10)] %*% ga)
  mu  <- (exp(eta)^(-1)+1)^(-1) 
  mu[mu > 1-.Machine$double.eps] <- 1-.Machine$double.eps
  mu[mu < .Machine$double.eps]   <- .Machine$double.eps
  
  numer  <- exp(eta)
  denom  <- (1+exp(eta))^2
  d_mu   <- as.vector(numer/denom)
  d_mu[mu < .Machine$double.eps   | mu == .Machine$double.eps]   <- .Machine$double.eps
  d_mu[mu > 1-.Machine$double.eps | mu == 1-.Machine$double.eps] <- .Machine$double.eps
  
  numer  <- exp(eta)*(1-exp(eta))
  denom  <- (1+exp(eta))^3
  d_2_mu <- numer/denom
  d_2_mu[mu < .Machine$double.eps   | mu == .Machine$double.eps]      <- .Machine$double.eps
  d_2_mu[mu > 1-.Machine$double.eps | mu == 1-.Machine$double.eps]    <- .Machine$double.eps
  
  numer  <- exp(eta)*(exp(eta)-(2+sqrt(3)))*(exp(eta)-(2-sqrt(3)))
  denom  <- (1+exp(eta))^3
  d_3_mu <- as.numeric(numer/denom)
  d_3_mu[mu < .Machine$double.eps   | mu == .Machine$double.eps]      <- .Machine$double.eps
  d_3_mu[mu > 1-.Machine$double.eps | mu == 1-.Machine$double.eps]    <- .Machine$double.eps
  
  g_beta_diag <- as.numeric((log(y/(1-y))-digamma(mu*phi_old_sgd)+digamma((1-mu)*phi_old_sgd))*d_3_mu -
    3*(phi_old_sgd*trigamma(mu*phi_old_sgd) + phi_old_sgd*trigamma((1-mu)*phi_old_sgd))*d_mu*d_2_mu-
    (phi_old_sgd^2*psigamma(mu*phi_old_sgd, deriv = 2) - phi_old_sgd^2 * psigamma((1-mu)*phi_old_sgd, deriv=2))*(d_mu)^3)
  g_beta <- apply(x[,1:5], 2, 
                   FUN = function(w) sum(diag( phi_old_sgd * self_solve(d_h_2(x,y,ga,bet,tau)) %*% t(x[,6:10]) %*% diag(g_beta_diag*w) %*% x[,6:10])) ) 
  
  p1 <- as.numeric((log(y/(1-y)) - phi_old_sgd*mu*trigamma(mu*phi_old_sgd) + (1-mu)*phi_old_sgd*trigamma((1-mu)*phi_old_sgd) - digamma(mu*phi_old_sgd) + digamma((1-mu)*phi_old_sgd))*d_2_mu)
  p2 <- as.numeric((phi_old_sgd*trigamma(mu*phi_old_sgd) + phi_old_sgd*trigamma((1-mu)*phi_old_sgd) + phi*trigamma(mu*phi_old_sgd) + mu*phi_old_sgd^2*psigamma(mu*phi_old_sgd, deriv=2) + phi_old_sgd^2*(1-mu)*psigamma((1-mu)*phi_old_sgd, deriv=2) + phi_old_sgd*trigamma((1-mu)*phi_old_sgd))*d_mu^2)
  g_phi  <- sum(diag( self_solve(d_h_2(x,y,ga,bet,tau)) %*% t(x[,6:10]) %*% diag(as.numeric(p1-p2)) %*% x[,6:10] ))
  
  return(list( g_beta = g_beta, g_phi = g_phi ))
} 

# define a function to give working vector for IRLS algo 

working_vec <- function(x, y, ga, bet, tau){ 
  
  # x <- Dat[,3:12]
  # y <- Dat[,1]
  # ga  <- gamma_new
  # bet <- beta_new_sgd
  # tau <- tau_new_sgd
  # 
  # ga  <- gamma_temp[i,]
  # bet <- beta_old_sgd
  # tau <- tau_old_sgd
  
  x <- as.matrix(x)
  y <- as.matrix(y) ; y[y==1] <- 1-.Machine$double.eps; y[y==0] <- .Machine$double.eps
  ga  <- as.matrix(ga)
  bet <- as.matrix(bet)
  tau   <- as.matrix(tau)
  
  d_h <- list(NULL)
  H   <- list(NULL)
  
  for (i in 1:length(unique(Dat$group))){
    
    Inx <- which(Dat$group==i)
    
    eta <- as.vector(x[Inx,1:5] %*% bet + x[Inx,6:10] %*% ga[i,])
    mu  <- (exp(eta)^(-1)+1)^(-1)
    mu[mu > 1-.Machine$double.eps] <- 1-.Machine$double.eps
    mu[mu < .Machine$double.eps]   <- .Machine$double.eps
    
    numer  <- exp(eta)
    denom  <- (1+exp(eta))^2
    d_mu   <- as.numeric(numer/denom)
    d_mu[mu < .Machine$double.eps   | mu == .Machine$double.eps]   <- .Machine$double.eps
    d_mu[mu > 1-.Machine$double.eps | mu == 1-.Machine$double.eps] <- .Machine$double.eps
    
    numer  <- exp(eta)*(1-exp(eta))
    denom  <- (1+exp(eta))^3
    d_2_mu <- as.numeric(numer/denom)
    d_2_mu[mu < .Machine$double.eps   | mu == .Machine$double.eps]   <- .Machine$double.eps
    d_2_mu[mu > 1-.Machine$double.eps | mu == 1-.Machine$double.eps] <- .Machine$double.eps
    
    A <- (log(y[Inx]/(1-y[Inx]))-digamma(mu*phi_new_sgd)+digamma((1-mu)*phi_new_sgd))*d_2_mu
    B <- (-phi_new_sgd*psigamma(mu*phi_new_sgd, deriv = 2)-phi_new_sgd*psigamma((1-mu)*phi_new_sgd))*d_mu^2
    w <- as.vector(A+B)
    
    J   <- phi_new_sgd * t(x[Inx,6:10]) %*% diag(d_mu) %*% (log(y[Inx]/(1-y[Inx])) - digamma(mu*phi_new_sgd) + digamma((1-mu)*phi_new_sgd)) - tau %*% ga[i,]
    d_J <- phi_new_sgd * t(x[Inx,6:10]) %*% diag(w) %*% x[Inx,6:10] - tau
    
    d_h[[i]] <- J
    H[[i]]   <- d_J
    
    # Fisher Scoring
    # I <- (t(J*x[Dat$group==i, 6:10]) - as.numeric(tau %*% ga[i,])) %*% t(t(J*x[Dat$group==i, 6:10]) - as.numeric(tau %*% ga[i,]))
    # H[[i]]   <- I/length(Inx)
    
  }
  
  return(list(d_h=d_h, H=H))
}

self_solve <- function(M){
  return(eigen(M)$vectors %*% diag(1/eigen(M)$values) %*% t(eigen(M)$vectors))
}
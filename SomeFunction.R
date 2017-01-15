# -----------------------------------------------------------------------------------------
# define some function to make the SGD neat

epsilon <- 1e-4

g_h <- function(x, y, ga, bet){
  
  x <- as.matrix(x)
  y <- as.matrix(y) ; y[y==1] <- 1-epsilon; y[y==0] <- epsilon
  ga  <- as.matrix(ga)
  bet <- as.matrix(bet)
  
  eta <- x[,1:5] %*% bet + x[,6:10] %*% ga
  mu  <- (exp(eta)^(-1)+1)^(-1) 
  
  numer  <- exp(eta)
  denom  <- (1+exp(eta))^2
  d_mu   <- numer/denom
  d_mu[mu < .Machine$double.eps]   <- .Machine$double.eps
  d_mu[mu > 1-.Machine$double.eps] <- .Machine$double.eps
  
  gradient <- diag(as.vector((phi*log(y/(1+y)) - phi*log(mu*phi) + (1/2-mu*phi)/(mu) + 
                                phi*log((1-mu)*phi) - (mu*phi-phi+1/2)/(1-mu))*d_mu)) %*% x[,1:5]
  gradient[mu < .Machine$double.eps]   <- .Machine$double.eps
  gradient[mu > 1-.Machine$double.eps] <- 1/2
  
  return(gradient)
} 

# before calculating that cumbersome matrix derivative wrt to vector, define h_2 in advance

# input x and y should be vector within a cluster

d_h_2 <- function(x, y, ga, bet, D){
  
  x <- as.matrix(x)
  y <- as.matrix(y) ; y[y==1] <- 1-epsilon; y[y==0] <- epsilon
  ga  <- as.matrix(ga)
  bet <- as.matrix(bet)
  
  eta <- as.vector(x[,c(1:5)] %*% bet + x[,c(6:10)] %*% ga)
  mu  <- (exp(eta)^(-1)+1)^(-1) 
  
  numer  <- exp(eta)
  denom  <- (1+exp(eta))^2
  d_mu   <- numer/denom
  d_mu[mu < .Machine$double.eps]   <- .Machine$double.eps
  d_mu[mu > 1-.Machine$double.eps] <- .Machine$double.eps
  
  numer  <- exp(eta)*(1-exp(eta))
  denom  <- (1+exp(eta))^3
  d_2_mu <- numer/denom
  d_2_mu[numer < .Machine$double.eps & denom < .Machine$double.eps] <- .Machine$double.eps
  d_2_mu[numer==-Inf & denom==Inf] <- .Machine$double.eps
  
  d_h_2 <- (phi*log(y/(1-y)) - phi*log(mu*phi) - (1/2-mu*phi)/mu + phi*log((1-mu)*phi) - (mu*phi + 1/2 - phi)/(1-mu))*d_2_mu +
    (1/(2*mu^2) - phi/mu - phi/(1-mu) - 1/(2*(1-mu)^2))*(d_mu)^2 
  d_h_2[mu < .Machine$double.eps]   <- .Machine$double.eps
  d_h_2[mu > 1-.Machine$double.eps] <- .Machine$double.eps
  
  return(t(x[,c(6:10)]) %*% diag(as.vector(d_h_2)) %*% x[,c(6:10)] - self_solve(D))
}

# now let's deal with the trace of a rank 3 tensor 

g_d_h_2 <- function(x, y, ga, bet, D){
  
  x <- as.matrix(x)
  y <- as.matrix(y); y[y==1] <- 1-epsilon ; y[y==0] <- epsilon
  ga  <- as.matrix(ga)
  bet <- as.matrix(bet)
  
  eta <- as.vector(x[,c(1:5)] %*% bet + x[,c(6:10)] %*% ga)
  mu  <- (exp(eta)^(-1)+1)^(-1) 
  
  numer  <- exp(eta)
  denom  <- (1+exp(eta))^2
  d_mu   <- numer/denom
  d_mu[mu < .Machine$double.eps]   <- .Machine$double.eps
  d_mu[mu > 1-.Machine$double.eps] <- .Machine$double.eps
  
  numer  <- exp(eta)*(1-exp(eta))
  denom  <- (1+exp(eta))^3
  d_2_mu <- numer/denom
  d_2_mu[numer==0 & denom==0]      <- .Machine$double.eps
  d_2_mu[numer==-Inf & denom==Inf] <- .Machine$double.eps
  
  numer  <- exp(eta)*(exp(eta)-(2-sqrt(3)))*(exp(eta)-(2+sqrt(3)))
  denom  <- (exp(eta)+1)^4
  d_3_mu <- numer/denom
  d_3_mu[numer==0 & denom==0]     <- .Machine$double.eps
  d_3_mu[numer==Inf & denom==Inf] <- .Machine$double.eps
  
  G_d_h_2 <- ((phi*log(y/(1-y))) - phi*log(mu*phi) - (1/2-mu*phi)/mu + phi*log((1-mu)*phi) - (mu*phi +1/2 - phi)/(1-mu))*d_3_mu+
    (-phi/(1-mu) + (1/2-mu*phi)/mu^2 - (1/2-phi+mu*phi)/(1-mu)^2 + 1/mu^2 - 2*phi/mu - 2*phi/(1-mu) - 1/(1-mu)^2)*d_mu*d_2_mu+
    (-1/(1-mu)^3 - 1/(2*mu^2) - phi/(1-mu)^2 + phi/mu^2)*(d_mu)^3
    
  G_d_h_2[mu < .Machine$double.eps]     <- 1.5
  G_d_h_2[mu > 1 - .Machine$double.eps] <- -2.5
  
  tensor_tr <- -apply(x[,1:5], 2, FUN = function(w) sum(diag( self_solve(d_h_2(x, y, ga, bet, D)) * t(x[,6:10]) %*% diag(as.vector(G_d_h_2*w)) %*% x[,6:10])))
  gradient  <- tensor_tr
    
  return(gradient)
  
}

g_psi <- function(x, y, ga, bet, D){ 
  
  x <- as.matrix(x)
  y <- as.matrix(y) ; y[y==1] <- 1-epsilon; y[y==0] <- epsilon
  ga  <- as.vector(ga)
  bet <- as.matrix(bet)
  
  eta <- as.vector(x[,c(1:5)] %*% bet + x[,c(6:10)] %*% ga)
  mu  <- (exp(eta)^(-1)+1)^(-1)  
  
  numer <- exp(eta)*(1-exp(eta))
  denom <- (1+exp(eta))^3
  d_mu  <- numer/denom
  d_mu[numer==-Inf & denom==Inf] <- .Machine$double.eps
  d_mu[numer==0 & denom==0] <- .Machine$double.eps
  
  d_h_tau   <- -1/2*ga %*% t(ga)
  d_D_tau   <- -t(self_solve(tau))
  d_h_2_tau <- -t(self_solve(d_h_2(x, y, ga, bet, tau)))
  
  gradient <- d_h_tau - (1/2)*d_D_tau + (1/2)*d_h_2_tau
  
  return(gradient)
}

# define a function to give working vector for IRLS algo 

working_vec <- function(x, y, ga, bet, D){ 
  
  x <- as.matrix(x)
  y <- as.matrix(y) ; y[y==1] <- 1-epsilon; y[y==0] <- epsilon
  ga  <- as.matrix(ga)
  bet <- as.matrix(bet)
  D   <- as.matrix(D)
  
  d_h <- list(NULL)
  H   <- list(NULL)
  
  for (i in 1:length(unique(Dat$group))){
    
    Inx <- which(Dat$group==i)
    
    eta <- as.vector(x[which(Dat$group==i),1:5] %*% bet + x[which(Dat$group==i),6:10] %*% ga[i,])
    mu  <- (exp(eta)^(-1)+1)^(-1)
    
    numer  <- exp(eta)
    denom  <- (1+exp(eta))^2
    d_mu   <- numer/denom
    d_mu[mu < .Machine$double.eps]   <- .Machine$double.eps
    d_mu[mu > 1-.Machine$double.eps] <- .Machine$double.eps
    
    numer  <- exp(eta)*(1-exp(eta))
    denom  <- (1+exp(eta))^3
    d_2_mu <- numer/denom
    d_2_mu[numer==0 & denom==0] <- .Machine$double.eps
    d_2_mu[numer==-Inf & denom==Inf] <- .Machine$double.eps
    
    numer  <- exp(eta)*(exp(eta)-(2-sqrt(3)))*(exp(eta)-(2+sqrt(3)))
    denom  <- (exp(eta)+1)^4
    d_3_mu <- numer/denom
    d_3_mu[numer==0 & denom==0]     <- .Machine$double.eps
    d_3_mu[numer==Inf & denom==Inf] <- .Machine$double.eps
    
    J <- (phi*log(y[Inx]/(1-y[Inx]))-phi*log(mu*phi)+(1/2-mu*phi)/mu+
            phi*log((1-mu)*phi)-(mu*phi+1/2-phi)/(1-mu))*d_mu
    J[mu < .Machine$double.eps]   <- .Machine$double.eps
    J[mu > 1-.Machine$double.eps] <- 1/2
    
    d_J <- (phi*log(y[Inx]/(1-y[Inx]))-phi*log(mu*phi)-(1/2-mu*phi)/mu + phi*log((1-mu)*phi) - (mu*phi+1/2-phi)/(1-mu))*d_2_mu+
      (1/(2*mu)-phi/mu-phi/(1-mu)-1/(2*(1-mu)^2))*d_mu^2
    d_J[mu < .Machine$double.eps]   <- 0
    d_J[mu > 1-.Machine$double.eps] <- 0
  
    d_h[[i]] <- apply(J * x[Dat$group==i,6:10],2,sum) - self_solve(D) %*% ga[i,]
    H[[i]]   <- t(x[Dat$group==i,6:10]) %*% diag(as.vector(d_J)) %*% x[Dat$group==i,6:10] - self_solve(D)
    
  }
  
  return(list( d_h=d_h, H=H))
}

self_solve <- function(M){
  return(eigen(M)$vectors %*% diag(1/eigen(M)$values) %*% t(eigen(M)$vectors))
}
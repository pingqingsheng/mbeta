library(MASS)
# library(doParallel)
library(doMC) #if you use Mac or Linux
library(bigmemory)
library(bit64)
set.seed(123)

rm(list=ls())

# Simulate Some Data y ~ beta(mu*phi, (1-mu)phi) | g(mu) = logit(X*Beta + Z*gamma)
fixed   <- mvrnorm(1,m=rep(0.0001,5), Sigma=diag(rep(0.1,5)))
D       <- diag(rep(10, length(fixed)))
random  <- as.matrix(mvrnorm(100, m=rep(0.001,5), Sigma=D))
random  <- matrix(rep(random, each=100),ncol=5, byrow=FALSE)

phi     <- 1

X <- matrix(0,10000,4)
for(i in 1:100){
  X[((i-1)*100+1):(i*100),] <- mvrnorm(100, m=rep(i/100,4), Sigma=diag(rep(1,4)))
}
X <- cbind(1,X)
Z <- X

group <- rep(seq_len(100), each=100)
Dat   <- as.data.frame(cbind(X,Z))
Dat   <- cbind(group, Dat)
# Let X=Z at present time, but of course we could sim different Z
eta <- X %*% fixed + apply(Z*random, 1, sum)
mu  <- exp(eta)/(1+exp(eta))
mu  <- pmin(mu, 1-.Machine$double.eps)
mu  <- pmax(mu, .Machine$double.eps )
mu[is.nan(mu)] <- 1-.Machine$double.eps
Y   <- rbeta(length(group), shape1 = mu*phi, shape2 = (1-mu)*phi)
Y[Y==1] <- 1-.Machine$double.eps # deal with NaN and Inf problem
Y[Y==0] <- .Machine$double.eps   # deal with NaN and Inf problem
names(Dat) <- c("group",paste0("x",seq_len(5)),paste0("z",seq_len(5)))
Dat <- cbind(Y=Y, Dat)

# Initialization
# Use mini batch (possiblly find a extra package do SGD to accelerate or rewrite with C)
terminate   <- FALSE
tolerance_1 <- 1 # tolerance for fixed effect
tolerance_2 <- 1 # tolerance for random effect
tolerance_3 <- 5 # tolerance for covariance component

# Starting Value (Start with the real value)
# gamma_new <- unique(random)
# beta_new  <- fixed
# sig_new   <- as.numeric(diag(D))
# phi_new   <- 5

gamma_new <- matrix(rep(100, prod(dim(unique(random)))), dim(unique(random))[1], dim(unique(random))[2])
beta_new  <- rep(100, dim(X)[2])
sig_new   <- as.numeric(rep(100, dim(gamma_new)[2]))
phi_new   <- 5


rate   <- 0.01 # Learning Rate
maxIte <- 30
cond_1 <- FALSE
cond_2 <- FALSE
cond_3 <- FALSE

marker <- 1
verbose  <- TRUE # For monitor purpose
parallel <- FALSE # For parallel

while(!terminate){

  # Laplace approximate random effect

  # Use IRWLS here
  steps <- 1 # Initial Steps
  if(!parallel){

    for(i in 1:length(unique(group))){

      x   <- X[group==i, ]
      z   <- Z[group==i, ]
      y   <- Y[group==i]

      gamma_temp <- gamma_new[unique(group) == i,]

      marker_inner <- 0
      while(max(steps) > 0.1 & marker_inner < 100){

        Gradient <- g_gamma(x, z, y,  beta_new, sig_new, gamma_temp, phi_new)
        Hessian  <- g_gamma_2(x, z, y,  beta_new, sig_new, gamma_temp, phi_new)

        # Transform Newtown-Rapson to a IRWLS problem and use QR decomposite to speed the algo
        qrHessian   <- qr(Hessian)
        steps       <- backsolve(qr.R(qrHessian), diag(rep(1,dim(Hessian)[1]))) %*% t(qr.Q(qrHessian)) %*% Gradient
        (gamma_temp <- gamma_temp - t(steps))

        marker_inner <- marker_inner + 1

      }

      gamma_new[i,] <- gamma_temp

    }

    if(verbose){

      cat("Max group (IRWLS) random effect:", max(gamma_temp), "\n")

    }

  }else{

    registerDoMC(20)
    # registerDoMC(length(unique(group)))

    gamma_new <- foreach(i = seq_len(max(group))) %dopar% {

      x   <- X[group==i, ]
      z   <- Z[group==i, ]
      y   <- Y[group==i]

      gamma_temp <- gamma_new[unique(group) == i,]

      marker_inner <- 0
      while(max(steps) > 0.1 & marker_inner < 200){

        Gradient <- g_gamma(x, z, y,  beta_new, sig_new, gamma_temp, phi_new)
        Hessian  <- g_gamma_2(x, z, y,  beta_new, sig_new, gamma_temp, phi_new)

        # Transform Newtown-Rapson to a IRWLS problem and use QR decomposite to speed the algo
        qrHessian <- qr(Hessian)
        steps <- solve(qr.R(qrHessian)) %*% t(qr.Q(qrHessian)) %*% Gradient
        gamma_temp <- gamma_temp - t(steps)

        marker_inner <- marker_inner + 1

      }

      return(gamma_temp)

    }

    registerDoSEQ()

    gamma_new <- matrix(unlist(gamma_new), ncol=dim(Z)[2], byrow=FALSE)
  } # Use DoParallel

  # Stochastic gradient descent to get fixed effect | phi | covariance compoent

  # Use ADAM here to avoid explosive gradient
  # ADAM parameter setting
  a_beta <- 0.9; a_sig  <- 0.9; a_phi  <- 0.9
  m_beta <- 0;   m_sig  <- 0;   m_phi  <- 0
  v_beta <- 0;   v_sig  <- 0;   v_phi  <- 0

  beta_sgd <- beta_old <- beta_new
  phi_sgd  <- phi_old  <- phi_new
  sig_sgd  <- sig_old  <- sig_new

  marker_inner <- 0
  converge <- FALSE
  while(!converge){

    Inx <- sample(unique(group),1)
    x   <- X[group==Inx, ]
    z   <- Z[group==Inx, ]
    y   <- Y[group==Inx]

    beta_gradient  <- g_beta(x, z, y, beta_sgd, sig_sgd, gamma_new[Inx,], phi_sgd)
    sigma_gradient <- g_sigma(gamma_new[Inx,], sig_sgd)
    phi_gradient   <- g_phi(x, z, y, beta_sgd, sig_sgd, gamma_new[Inx,], phi_sgd)

    m_beta <- (a_beta*m_beta + (1-a_beta)*beta_gradient)/(1-a_beta^(marker_inner+1))
    v_beta <- (a_beta*v_beta + (1-a_beta)*beta_gradient^2)/(1-a_beta^(marker_inner+1))
    m_sig  <- (a_sig*m_sig + (1-a_sig)*sigma_gradient)/(1-a_sig^(marker_inner+1))
    v_sig  <- (a_sig*v_sig + (1-a_sig)*sigma_gradient^2)/(1-a_sig^(marker_inner+1))
    m_phi  <- (a_phi*m_phi + (1-a_phi)*phi_gradient)/(1-a_phi^(marker_inner+1))
    v_phi  <- (a_phi*v_phi + (1-a_phi)*phi_gradient^2)/(1-a_phi^(marker_inner+1))

    step_beta <- rate*m_beta/(sqrt(v_beta))
    step_sig  <- rate*m_sig/(sqrt(v_sig))
    step_phi  <- rate*m_phi/(sqrt(v_phi))

    beta_sgd <- beta_sgd -  step_beta
    sig_sgd  <- sig_sgd  -  step_sig
    phi_sgd  <- phi_sgd  -  step_phi

    if(max(abs(step_beta)) < 0.01 & max(abs(step_sig)) < 0.01 & step_phi < 0.01 | marker_inner > 500) {

      beta_new <- beta_sgd
      phi_new  <- phi_sgd
      sig_new  <- sig_sgd

      converge <- TRUE
    }

    marker_inner <- marker_inner + 1

  }

  if(verbose){

    cat("Final descent result (maxsg_beta|max s_phi|max s_sigma)",  max(step_beta), max(step_sig), step_phi, "\n")

  }

  #  Test convergence
  cond_1 <- ( beta_new - beta_old < tolerance_1 )
  cond_2 <- ( phi_new  - phi_old  < tolerance_2 )
  cond_3 <- ( sig_new  - sig_old  < tolerance_3 )

  if(!any(c(cond_1,cond_2,cond_3)==FALSE)) {

    terminate <- TRUE

    cat(" Estimation Finished, Final Output:")

    (beta_hat <- as.numeric(beta_new))
    (phi_hat  <- phi_new)
    (sig_hat  <- diag(sig_new))
     gamma_hat <- as.numeric(gamma_new)

  }

  if( marker > maxIte ){

    terminate <- TRUE

    cat(" Maximum iteration reached: fail to converge ----- final output ")

    (beta_hat <- as.numeric(beta_new))
    (phi_hat  <- phi_new)
    (sig_hat  <- diag(sig_new))
    (gamma_hat <- as.numeric(gamma_new))

  }

  marker <- marker + 1

}

gamma_hat  <- matrix(rep(gamma_new,100), ncol=5, byrow=FALSE)
eta_hat    <- as.matrix(Dat[,3:7]) %*% beta_hat + apply(as.matrix(Dat[,8:12])*gamma_hat, 1, sum)
mu_hat   <- 1/(1+exp(eta_hat)^(-1))

y[y > 1-.Machine$double.eps] <- 1-.Machine$double.eps
y[y < .Machine$double.eps]    <- .Machine$double.eps

mu_real <- y

RMSE <- sqrt(sum((mu_hat - mu_real)^2))

# dev.cur()
# png('trialserver.png', width=1000, height=500, res=100, units="px")
plot(density(y), type="l", lwd=2)
lines(density(mu_hat, bw=density(y)$bw), lwd=2, lty=2, main="Line: Real vs Fitted", col="red")
legend("top", legend=c("Fitted","Real"), lty=c(2,1), lwd=2, col=c(2,1))
# dev.off()



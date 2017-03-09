
# Simulate some data
library(MASS)
set.seed(234)
fixed   <- mvrnorm(1,m=rep(5,5), Sigma=diag(rep(1,5)))
D <- outer(rep(sqrt(3),5), rep(sqrt(3),5))
diag(D) <- 9
random <- as.matrix(mvrnorm(100, m=rep(0,5), Sigma=D))
random <- matrix(rep(random, each=100),ncol=5, byrow=FALSE)
phi <- 5

X <- matrix(0,10000,4)
for(i in 1:100){
  X[((i-1)*100+1):(i*100),] <- mvrnorm(100, m=rep(i/100,4), Sigma=diag(rep(1,4)))
}
X <- cbind(1,X)
Z <- X

group <- rep(seq_len(100), each=100)
Dat <- as.data.frame(cbind(X,Z))
Dat <- cbind(group, Dat)

#  Let X=Z at present time, but of course we could sim different Z
eta <- X %*% fixed + apply(Z*random, 1, sum)
mu  <- as.vector(pmax(1/(1+exp(eta)^(-1)), .Machine$double.eps))
y   <- rbeta(length(group), shape1 = mu*phi, shape2 = (1-mu)*phi)
names(Dat) <- c("group",paste0("x",seq_len(5)),paste0("z",seq_len(5)))
Dat <- cbind(y=y, Dat)

# Estimation with Laplace Approximation 

# initialization 

# Use mini batch (possiblly find a extra package do SGD to accelerate or rewrite with C)
terminate   <- FALSE 
tolerance_1 <- 1
tolerance_2 <- 1
tolerance_3 <- 5
# gamma_new   <- gamma_0 
gamma_new    <- unique(random)
beta_new_sgd <- fixed
# D_new       <- matrix(15,5,5)
# diag(D_new) <- 20
D_new  <- D
phi_new_sgd <- 5
  
rate   <- 0.01
cond_1 <- FALSE
cond_2 <- FALSE
cond_3 <- FALSE

marker <- 1
verbose <- TRUE

# SGD then Newtown-Rapson iteratively
while(!terminate){
  
  beta_old <- beta_new_sgd
  phi_old  <- phi_new_sgd
  D_old    <- D_new
  tau_new_sgd  <- self_solve(D_old)
  
  gamma_old  <- matrix(rep(gamma_new, each=100), ncol=5, byrow  = FALSE)
  gamma_temp <- aggregate(gamma_old, by=list(Dat$group), unique, simplify=TRUE)
  gamma_temp <- as.matrix(gamma_temp[,-1]) 
  
  # Use ADAM updating scheme 
  b_1 <- 0.9
  b_2 <- 0.9
  
  m_beta  <- 0
  nu_beta <- 0
  m_D  <- 0
  nu_D <- 0
  m_phi  <- 0
  nu_phi <- 0
  
  # m_beta_old  <- 0
  # nu_beta_old <- 0
  # m_D_old  <- 0
  # nu_D_old <- 0
  
  # Give gamma_i (gamma_old), use SGD to optimize beta and D
  for (i in 1:max(as.integer(Dat$group))){
    
    # reset parameter and moment
    beta_old_sgd  <- beta_new_sgd
    tau_old_sgd   <- tau_new_sgd
    phi_old_sgd   <- phi_new_sgd
    
    X <- as.matrix(Dat[group==i, 3:12])
    Y <- as.matrix(Dat[group==i, 1])
    
    # ------------------------ ADAM SGD ---------------------------------
    # Updating 1st and 2nd moment
    
    score_beta <- g_h(X, Y, gamma_temp[i,], beta_old_sgd)$g_beta - 1/2*g_d_h_2(X, Y, gamma_temp[i,], beta_old_sgd, tau_old_sgd)$g_beta
    m_beta  <- b_1*m_beta  + (1-b_1)*score_beta
    nu_beta <- b_2*nu_beta + (1-b_2)*score_beta^2
    steps_1      <- rate/sqrt(nu_beta/(1-b_1))*m_beta/(1-b_1)
    beta_new_sgd <- beta_old_sgd + steps_1
    
    score_phi <- g_h(X, Y, gamma_temp[i,], beta_old_sgd)$g_phi - 1/2*g_d_h_2(X, Y, gamma_temp[i,], beta_old_sgd, tau_old_sgd)$g_phi
    m_phi     <- b_1*m_phi  + (1-b_1)*score_phi
    nu_phi    <- b_2*nu_phi + (1-b_2)*score_phi^2
    steps_2   <- rate/sqrt(nu_phi/(1-b_1))*m_phi/(1-b_1)
    phi_new_sgd <- phi_old_sgd + steps_2

    score_tau  <- g_h(X, Y, gamma_temp[i,], beta_old_sgd)$g_tau - 1/2*t(self_solve(d_h_2(X,Y,gamma_temp[i,],beta_old_sgd,tau_old_sgd)))  + 1/2*D_old
    m_D  <- b_1*m_D + (1-b_1)*score_tau
    nu_D <- b_2*nu_D + (1-b_2)*score_tau^2
    steps_3      <- rate*(1/sqrt(nu_D/(1-b_2)))*m_D/(1-b_2)
    tau_new_sgd  <- tau_old_sgd + steps_3
    tau_new_sgd  <- eigen(tau_new_sgd)$vectors %*% diag(pmax(eigen(tau_new_sgd)$values,0)) %*% t(eigen(tau_new_sgd)$vectors)

    # 
    # ----------------------------- Momentum SGD -----------------------------------
    
    # nu_beta_new <- 0.9*nu_beta_old + rate*score_beta
    # nu_D_new    <- 0.9*nu_D_old + rate*score_psi

    # beta_new_sgd <- beta_old_sgd + nu_beta_new
    # D_new_sgd    <- D_old_sgd + nu_D_new
    #   
    if(any(is.nan(beta_new_sgd)==TRUE)) stop("NaN appears again !!!!!!!!!!!")
    
    cond1 <- all(abs(beta_old_sgd - beta_new_sgd) < 0.00001 )
    cond2 <- all(abs(phi_old_sgd - phi_new_sgd) < 0.00001 )
    cond3 <- all(abs(tau_old_sgd - tau_new_sgd) < 0.1 )
    
    if(cond1 & cond2 & cond3) {cat("SGD converge","\n"); break}
    
    if(verbose){
      cat(i,"th group descent (beta|phi)", paste0(round(beta_new_sgd,digit=4), collapse=" "), " | ", phi_new_sgd, "\n")
      # cat("biggest tau_ij", max(tau_new_sgd),"\n")
    }
    
  }
  
  # Give up iteratively reweighted algorithm.....
  # I need some time to really understand IRWLS
  # Use Newtown Rapson instead
  marker_inner <- 1
  converge <- FALSE
  
  gamma_old     <- matrix(rep(gamma_temp,each=100), ncol=5, byrow=FALSE)
  gamma_new     <- aggregate(gamma_old, by=list(Dat$group), unique, simplify=TRUE)
  gamma_new     <- as.matrix(gamma_new[,-1])
  
  while(!converge){
    
    gamma_old     <- gamma_new
    working_value <- working_vec(Dat[,3:12], Dat[,1], gamma_old, beta_new_sgd, tau_new_sgd)
    
    # if(marker_inner > 100) {
    #   factor <- marker_inner
    # }else{
    #   factor <- 101
    # }
    # working_value <- working_vec(Dat[,3:12], Dat[,1], gamma_old, beta_new_sgd, tau_new_sgd)
    
    # All these loops seems stupid
    for (i in 1:length(unique(Dat$group))){
      
      gamma_new[i,] <- gamma_old[i,] + self_solve(working_value$H[[i]]) %*% working_value$d_h[[i]]
      
    }
    
    if(all(abs(gamma_new - gamma_old) < 2)) {converge <- TRUE}
    
    if(any(is.nan(gamma_new)==TRUE)) {stop("NaN appears again !!!!!!!!!!")}
    
    marker_inner <- marker_inner + 1
    cat(paste0(marker_inner, "th iteration: Newtown's Way ", max(abs(gamma_new - gamma_old))), "\n")
    cat("Difference between D and D_hat(frobenious): ", sum((D_old-var(gamma_new))^2), "\n")
    
    # if(marker_inner > 19) {cat("Fail to converge","\n"); break}
    
  }
  
  # check convergence
  
  beta_new <- beta_new_sgd
  phi_new  <- phi_new_sgd
  D_new    <- self_solve(tau_new_sgd)
  
  cond_1 <- (all(abs(beta_new - beta_old) < tolerance_1))
  cond_2 <- (all(abs(phi_new - phi_old) < tolerance_2))
  cond_3 <- (sum((D_new-D_old)^2) < tolerance_3)
  
  if ( cond_1 & cond_2 & cond_3) {
    
    terminate <- TRUE 
    cat("Estimation Finished","\n")
    
    (beta_hat <- beta_new)
    (D_hat    <- D_new)
    (phi_hat  <- phi_new)
    
  } 
  
  marker <- marker + 1
  
  if(marker > 30) {

    warning("Fail to converge, estimation is not optimal")
    terminate <- TRUE
    
    (beta_hat <- as.numeric(beta_new))
    (D_hat    <- D_new)
    (phi_hat  <- phi_new)

  }
  
  cat("--------",paste0(marker, "th iteration"),"--------","\n")
}

# 
gamma_hat  <- matrix(rep(gamma_new,100), ncol=5, byrow=FALSE)
eta_hat    <- as.matrix(Dat[,3:7]) %*% beta_hat + apply(as.matrix(Dat[,8:12])*gamma_hat, 1, sum)
mu_hat   <- 1/(1+exp(eta_hat)^(-1))

y[y > 1-.Machine$double.eps] <- 1-.Machine$double.eps
y[y < .Machine$double.eps]    <- .Machine$double.eps

mu_real <- y

RMSE <- sqrt(sum((mu_hat - mu_real)^2))

dev.cur()
png("trial .png", width=600, height=400, res=100, units="px")
plot(density(y), type="l", lwd=2)
lines(density(mu_hat, bw=density(y)$bw), lwd=2, lty=2, main="Line: Real vs Fitted", col="red")
legend("topright", legend=c("Fitted","Real"), lty=c(2,1), lwd=2, col=c(2,1))
dev.off()


# Simulate some data
library(MASS)
set.seed(234)
fixed  <- mvrnorm(1,m=rep(0,5), Sigma=diag(rep(1,5)))
D <- outer(rep(sqrt(3),5), rep(3,5))
diag(D) <- 9
random <- as.matrix(mvrnorm(100, m=rep(3,5), Sigma=D))
random <- matrix(rep(random, each=100),ncol=5, byrow=FALSE)

X <- matrix(0,10000,5)
for(i in 1:100){
  X[((i-1)*100+1):(i*100),] <- mvrnorm(100, m=rep(i/10,5), Sigma=diag(rep(1,5)))
}
Z <- X

group <- rep(seq_len(100), each=100)
Dat <- as.data.frame(cbind(X,Z))
Dat <- cbind(group, Dat)

#  Let X=Z at present time, but of course we could sim different Z
eta <- X %*% fixed + apply(Z*random, 1, sum)
mu  <- pmax(exp(eta)/(1+exp(eta)), .Machine$double.eps)
phi <- 10
y <- rbeta(10000, shape1 = mu*phi, shape2 = (1-mu)*phi)
names(Dat) <- c("group",paste0("x",seq_len(5)),paste0("z",seq_len(5)))
Dat <- cbind(y=y, Dat)

# Estimation with Laplace Approximation 

# initialization 

gamma_0 <- matrix(10, 100, 5)
beta_0  <- rep(3,5)
D_0     <- matrix(100,5,5)
diag(D_0) <- 100+50

# gamma_0 <- as.matrix(mvrnorm(100, m=rep(3,5), Sigma=D))
# beta_0  <- fixed
# D_0     <- D

# Use mini batch (possiblly find a extra package do SGD to accelerate or rewrite with C)
terminate   <- FALSE 
tolerance_1 <- 0.0001
tolerance_2 <- 0.01 
gamma_new   <- gamma_0 
beta_new_sgd  <- beta_0
D_new_sgd     <- D_0 
rate      <- 0.05
cond_1 <- FALSE
cond_2 <- FALSE

marker <- 1
verbose <- TRUE
while(!terminate){
  
  beta_old <- beta_new_sgd
  tau_old  <- tau_new
  
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
  
  # m_beta_old  <- 0
  # nu_beta_old <- 0
  # m_D_old  <- 0
  # nu_D_old <- 0
  
  # Give gamma_i (gamma_old), use SGD to optimize beta and D
  for (i in 1:max(as.integer(Dat$group))){
    
    # reset parameter and moment
    beta_old_sgd  <- beta_new_sgd
    tau_old_sgd   <- self_solve(D_old_sgd)
    
    X <- Dat[group==i, 3:12]
    Y <- Dat[group==i, 1]
    
    score_beta <- apply(g_h(X, Y, gamma_temp[i,], beta_old_sgd),2,sum) + 1/2 * g_d_h_2(X, Y, gamma_temp[i,], beta_old_sgd, tau_old_sgd)
    score_psi  <- g_psi(X, Y, gamma_temp[i,], beta_old_sgd, tau_old_sgd)
    
    # ------------------------ ADAM SGD ---------------------------------
    
    # Updating 1st and 2nd moment
    m_beta  <- b_1*m_beta  + (1-b_1)*score_beta
    nu_beta <- b_2*nu_beta + (1-b_2)*score_beta^2
    steps_1      <- rate/sqrt(nu_beta/(1-b_1))*m_beta/(1-b_1)
    beta_new_sgd <- beta_old_sgd + steps_1

    m_D  <- b_1*m_D + (1-b_1)*score_psi
    nu_D <- b_2*nu_D + (1-b_2)*score_psi^2
    steps_2  <- rate*self_solve(sqrt(nu_D/(1-b_2)))*m_D/(1-b_2)
    tau_new_sgd  <- tau_old_sgd + steps_2
    
    # 
    # ----------------------------- Momentum SGD -----------------------------------
    
    # nu_beta_new <- 0.9*nu_beta_old + rate*score_beta
    # nu_D_new    <- 0.9*nu_D_old + rate*score_psi
    # 
    # beta_new_sgd <- beta_old_sgd + nu_beta_new
    # D_new_sgd    <- D_old_sgd + nu_D_new
  #   
    if(any(is.nan(beta_new_sgd)==TRUE)) stop("NaN appears again !!!!!!!!!!!")

    cond_1 <- all(abs(beta_old_sgd - beta_new_sgd) < 0.001 )
    cond_2 <- all(abs(tau_old_sgd - tau_new_sgd) < 0.1 )

    if(cond_1 & cond_2) {cat("SGD converge","\n"); break}

    if(verbose){
      cat(paste0(i,"th group descent ",paste0(round(beta_new_sgd,digit=4),collapse=" "), "\n"))
      cat("biggest Dii", max(tau_new_sgd),"\n")
      Sys.sleep(0.1)
    }
  }
  
  # Give up iteratively reweighted algorithm.....
  # I need some time to really understand IRWLS
  # Use Newtown Rapson instead
  marker_inner <- 1
  converge <- FALSE
  
  gamma_new     <- aggregate(gamma_old, by=list(Dat$group), unique, simplify=TRUE)
  gamma_new     <- as.matrix(gamma_new[,-1])
  
  while(!converge){

    gamma_old     <- gamma_new
    working_value <- working_vec(Dat[,3:12], Dat[,1], gamma_old, beta_new_sgd, tau_new_sgd)
    
    # All these loops seems stupid
    for (i in 1:length(unique(Dat$group))){
      
      temp_inverse  <- eigen(working_value$H[[i]])$vector %*% diag(1/eigen(working_value$H[[i]])$values) %*% t(eigen(working_value$H[[i]])$vector)
      gamma_new[i,] <- gamma_old[i,] - temp_inverse %*% working_value$d_h[[i]]
      
    }
    
    if(all(abs(gamma_new - gamma_old) < 0.1)) {converge <- TRUE}
    
    if(any(is.nan(gamma_new)==TRUE)) {stop("NaN appears again !!!!!!!!!!")}

    marker_inner <- marker_inner + 1
    cat(paste0(marker_inner, "th iteration: Newtown's Way ", max(abs(gamma_new - gamma_old))), "\n")
    cat("Difference between D and D_hat(frobenious): ", sum((tau_new_sgd-var(gamma_new))^2), "\n")
    
    # if(marker_inner > 19) {cat("Fail to converge","\n"); break}

  }

  # check convergence
  beta_new <- beta_new_sgd
  tau_new  <- tau_new_sgd
  
  cond_1 <- (all(abs(beta_new - beta_old) < tolerance_1))
  cond_2 <- (sum((tau_new-tau_old)^2) < tolerance_2)
  
  if ( cond_1 & cond_2 ) {
    
    terminate <- TRUE 
    cat("Estimation Finished","\n")
    
    (beta_hat <- beta_new)
    (D_hat    <- self_solve(tau_new))
    
  } 
  
  marker <- marker + 1
  
  # if(marker > 29) {
  #   
  #   warning("Fail to converge, esti")
  #   terminate <- TRUE
  #   
  # }
  
  cat("--------",paste0(marker, "th iteration"),"--------","\n")
}





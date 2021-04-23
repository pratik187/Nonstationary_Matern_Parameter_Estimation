library(geoR)
library(LaplacesDemon)

NS.cov <- function(n = 100,locations = NULL,lam = NULL,sigma = NULL, nu = NULL){
  kernel.local <- array(0, dim = c(2, 2, n))
  for(i in 1 : n){
    kernel.local[, , i] <-  lam[i]^2*diag(c(1,1))
  }
  
  ##Calculate Matern form Nonstationary Covariance function 
  Sigma.mat <- matrix(rep(NA, n^2), nrow = n)
  Q.mat <- nu.mat <- cov <- matrix(rep(NA, n^2), nrow = n)
  Inv_ij <- matrix(rep(NA,4),2,2)
  for (i in 1:n) {
    Sigma.mat[i, i] <- 1
    Q.mat[i, i] <- 0
    Kernel_i <- kernel.local[, , i]
    det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
    if (i < n) {
      for (j in (i + 1):n) {
        Kernel_j <- kernel.local[, , j]
        det_j <- Kernel_j[1,1] * Kernel_j[2,2] - 
          Kernel_j[1,2] * Kernel_j[2,1]
        Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
        Inv_ij[1,1] <- Kernel_ij[2,2] 
        Inv_ij[2,2] <- Kernel_ij[1,1] 
        Inv_ij[2,1] <- - Kernel_ij[2,1] 
        Inv_ij[1,2] <- - Kernel_ij[1,2] 
        det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - 
          Kernel_ij[1,2] * Kernel_ij[2,1]
        x <- locations[i, ] - locations[j, ]
        Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij)
        Q.mat[i, j] <- sqrt(t(x) %*% Inv_ij %*% x/det_ij)
        Sigma.mat[j, i] <- Sigma.mat[i, j]
        Q.mat[j, i] <- Q.mat[i, j]
      }
    }
  }
  for(i in 1:n){
    for(j in 1:n){
      nu.mat[i,j] <- (nu[i]+nu[j])/2
      a = nu.mat[i,j]
      b = Q.mat[i,j]
      cov[i,j] <- geoR::matern(u = Q.mat[i,j], phi = 1, 
                         kappa = nu.mat[i,j])
      }
  }
  
  NS.cov <- Sigma.mat * cov
  
  return(NS.cov)
}



n <- 400
N.mc <- 9
x <- y <-seq(0,1,length.out = sqrt(n))
locations <- as.matrix(expand.grid(x,y))
lam_true <- 0.04*exp((sin(2*pi*locations[,1]/4)
                      +sin(2*pi*locations[,2]/4)))
sig_true <- exp(-locations[,1]-locations[,2])/3+0.8
nu_true <- exp((locations[,1]+locations[,2])/2)/10+0.3

simulation <- function(n = 100, m = 1, locations = NULL,lam = NULL,sigma = NULL, nu = NULL){
  
  kernel.local <- array(0, dim = c(2, 2, n))
  for(i in 1 : n){
    kernel.local[, , i] <-  lam[i]*diag(c(1,1))
  }
  
  ##Calculate Matern form Nonstationary Covariance function 
  Sigma.mat <- matrix(rep(NA, n^2), nrow = n)
  Q.mat <- nu.mat <- cov <-matrix(rep(NA, n^2), nrow = n)
  Inv_ij <- matrix(rep(NA,4),2,2)
  print("first iteration")
  for (i in 1:n) {
    Sigma.mat[i, i] <- sigma[i]^2
    Q.mat[i, i] <- 0
    Kernel_i <- kernel.local[, , i]
    det_i <- Kernel_i[1,1] * Kernel_i[2,2] - Kernel_i[1,2] * Kernel_i[2,1]
    if (i < n) {
      for (j in (i + 1):n) {
        Kernel_j <- kernel.local[, , j]
        det_j <- Kernel_j[1,1] * Kernel_j[2,2] - 
          Kernel_j[1,2] * Kernel_j[2,1]
        Kernel_ij <- 0.5 * (Kernel_i + Kernel_j)
        Inv_ij[1,1] <- Kernel_ij[2,2] 
        Inv_ij[2,2] <- Kernel_ij[1,1] 
        Inv_ij[2,1] <- - Kernel_ij[2,1] 
        Inv_ij[1,2] <- - Kernel_ij[1,2] 
        det_ij <- Kernel_ij[1,1] * Kernel_ij[2,2] - 
          Kernel_ij[1,2] * Kernel_ij[2,1]
        x <- locations[i, ] - locations[j, ]
        Sigma.mat[i, j] <- sqrt(sqrt(det_i * det_j)/det_ij) * sigma[i] * sigma[j]
        Q.mat[i, j] <- (t(x) %*% Inv_ij %*% x/det_ij)
        Sigma.mat[j, i] <- Sigma.mat[i, j]
        Q.mat[j, i] <- Q.mat[i, j]
        
        if( j %% 10 == 0){
          
          print(j)
          print("############################")
          
        }
      }
    }
    print("#####################################")
    
    print(paste0("iteration",i,"complete"))
    print("#####################################")
  }
  print("2nd iteration")
  for(i in 1:n){
    for(j in 1:n){
      nu.mat[i,j] <- (nu[i]+nu[j])/2
      val = 2*sqrt(Q.mat[i,j] * nu.mat[i,j])
      cov[i,j] <- matern(u = val, phi = 1, 
                               kappa = nu.mat[i,j])
      if( j %% 10 == 0)print(j)
    }
  }
  
  NS.cov <- Sigma.mat * cov
  set.seed(1)
  sim.data <- LaplacesDemon::rmvnc(m, rep(0,n), NS.cov)
  simdata <- list(locations = locations, kernel.local = kernel.local, 
                  NS.cov = NS.cov, sim.data = sim.data)
  
}


sims <- simulation(n = 400, m = 10, locations = locations,lam = lam_true, sigma = sig_true, nu = nu_true)
#NS_cov = NS.cov(n = n, locations = locations,lam = lam_true, sigma = sig_true, nu = nu_true)

data1 <- sims$sim.data[1,]
the_data = list(x = locations[,1], y = locations[,2], z = data)
df = data.frame(the_data)
df$z <- df$z * sig_true
min(data1)

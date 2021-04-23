N = 10000
x <- seq(0,1, length.out = 100)
y <- seq(0,1, length.out = 100)
d1 <- expand.grid(x = x, y = y)
X = d1$x              # X, Y co-ordinates getting generated here
Y = d1$y

co_ords = data.frame(X = X, Y = Y)    # data frame of the co-ordinates

# function for spatial range
lambda = function(u){
  return(0.04*exp(sin(0.5*pi*u[1])+sin(0.5*pi*u[2])))
}             

# function for partial sill
sigma = function(u){
  return((0.33 * exp(-(u[1]+u[2]))) + 0.8)
}             

# Function for smoothness 
neu = function(u){
  return(0.1 * exp(-0.5*(u[1]+u[2])) + 0.3)
}                


# non stationary matern cov function

NS_matern = function(u1, u2){
  cap_sigma_i = lambda(u1)*diag(2)
  cap_sigma_j = lambda(u2)*diag(2)
  term1 = sigma(u1)*sigma(u2)*(sqrt(lambda(u1)))*(sqrt(lambda(u2)))
  
  term2 = 2/(lambda(u1)+lambda(u2))
  
  # neuij = (neu(u1) + neu(u2))/2
  neuij = 0.98
  Qij = term2* (((u1[1]-u2[1])^2) + ((u1[2]-u2[2])^2))
  prod1 = 2*sqrt(neuij * Qij)
  term3 = ((prod1)^neuij) * besselK(prod1,neuij)
  
  return(term1*term2*term3)
}            

# generating the matrix of order N X N 

NS_Cov_matrix = matrix(, nrow = N, ncol = N)

for(i in 1:N){
  for(j in i:N){
    u1 <- unlist(co_ords[i,],use.names = FALSE)
    u2 <- unlist(co_ords[j,],use.names = FALSE)
    if(j == i){
      NS_Cov_matrix[i,j] <- sigma(u1)^2
    }else{
      value = NS_matern(u1,u2)
      NS_Cov_matrix[i,j] <- value
      NS_Cov_matrix[j,i] <- value
      
    }
    
  }
  if(i %% 10 == 0){
    print("######################################")
    print(paste0("row ", i," is complt"))
    print("######################################")
  }
  
}                          

library(MASS)
Z = mvrnorm(1,rep(0,N),NS_Cov_matrix)

co_ords$Z <- Z
# R = chol(NS_Cov_matrix)

library('exageostatr')
dmetric = "euclidean"

exageostat_init(hardware = list(ncores = 20, ngpus=0,
                                ts = 320, lts = 600, pgrid = 1, qgrid = 1))

# data = simulate_data_exact(sigma_sq, beta, nu, dmetric, n, seed)
result = exact_mle(list(x = co_ords[, 1],
                        y = co_ords[, 2],
                        z = co_ords[, 3]),
                   dmetric,
                   optimization = list(clb = c(0.001,0.001, 0.001),
                                       cub = c(5,5,5),
                                       tol = 1e-5,
                                       max_iters = 2000))

exageostat_finalize()
sink("final_result.txt")
print(list(result))
sink()

sigma = function(u){
  return((0.33 * exp(-(u[1]+u[2]))) + 0.8)
}

sigma_vec = rep(NA, N)

lambda_vec = rep(NA, N)
lambda = function(u){
  return(0.04*exp(sin(0.5*pi*u[1])+sin(0.5*pi*u[2])))
} 

for(i in 1:N){
  u = u1 <- unlist(co_ords[i,],use.names = FALSE)
  sigma_vec[i] = sigma(u)
}
for(i in 1:N){
  u = u1 <- unlist(co_ords[i,],use.names = FALSE)
  lambda_vec[i] = lambda(u)
}

quilt.plot(d1$x,d1$y, lambda_vec)

library(geoR)
data = read.csv('simulated_data_sigma.csv')
the_data = as.geodata(data,coords.col = 2:3, data.col = 4)
ml.n <- likfit(the_data , ini = c(0.81,0.7), nugget = 0, fix.nugget = TRUE)




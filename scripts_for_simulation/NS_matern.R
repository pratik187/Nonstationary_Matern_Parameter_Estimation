library(geoR)

N = 100
x <- seq(0, 1, length.out = 10)
y <- seq(0, 1, length.out = 10)
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
  
  neuij = (neu(u1) + neu(u2))/2
  # neuij = 0.95
  Qij = term2* (((u1[1]-u2[1])^2) + ((u1[2]-u2[2])^2))
  prod1 = 2*sqrt(neuij * Qij)
  term3 = matern(prod1, 1, neuij)
  
  return(term1*term2*term3)
}            

# generating the matrix of order N X N 

NS_Cov_matrix = matrix(, nrow = N, ncol = N)

for(i in 1:N){
  for(j in i:N){
    u1 <- unlist(co_ords[i,],use.names = FALSE)
    u2 <- unlist(co_ords[j,],use.names = FALSE)
    
      value = NS_matern(u1,u2)
      NS_Cov_matrix[i,j] <- value
      NS_Cov_matrix[j,i] <- value
    
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

setwd("C:/Users/Monalisha/Desktop/spatial \ stats/project work")
# write.table(x = NS_Cov_matrix ,file = "cov_matrix.csv", sep = ',', 
#             row.names = FALSE, col.names = FALSE)
write.csv(x = co_ords,file = "simulated_data_NS_matern.csv")




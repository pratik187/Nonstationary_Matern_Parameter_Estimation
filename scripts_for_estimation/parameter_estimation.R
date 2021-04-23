library(fields)

N = 10000
x <- seq(0,1, length.out = 100)
y <- seq(0,1, length.out = 100)
d1 <- expand.grid(x = x, y = y)
X = d1$x              # X, Y co-ordinates getting generated here
Y = d1$y

co_ords = data.frame(X = X, Y = Y) 

k = function(u1,u2,h){
  distance = ((u1[1] - u2[1])^2) + ((u1[2] - u2[2])^2)
  return(exp(-distance/(2*h)))
}


nod_points = list(n1 = c(0.25,0.25), n2 = c(0.25,0.75), n3 = c(0.75,0.25),
                  n4 = c(0.75,0.75))

param_est = function(u,nod_points,sigma,h){
  k_part = rep(NA,4)
  for(i in 1:4){
    
    k_part[i] = k(u,unlist(nod_points[i],use.names = FALSE),h)
  }
  sum_k = sum(k_part)
  param_est1 = 0
  for(i in 1:4){
    param_est1 = param_est1 + (sigma[i]*(k_part[i]/sum_k))
  }
  return(param_est1)
}

estimation_vec = function(param,nod_points, sigma_l, h ){
  count = 0
  original_param_vec = rep(NA,N)
  param_est_vec = rep(NA,N)
  for(i in 1:N){
    u1<- unlist(co_ords[i,],use.names = FALSE)
    original_param_vec[i] = param(u1)
    param_est_vec[i] = param_est(u1,nod_points,sigma_l,h)
    count = count +1
    if(count %% 50 == 0)print(count)
  }
  return(list(original = original_param_vec, estimated = param_est_vec))
}

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

##########################################################################
####### exact MLE
## Sigma estimate
sigma_node = c(1.0340602, 0.4044752, 0.6029332, 0.5432597)
sigma_node = sqrt(sigma_node)

exact_sigma_estimation_vec = estimation_vec(sigma, nod_points, sigma_node, h = 0.09)


#quilt.plot(d1$x,d1$y, abs(sigma_est_vec - sigma_vec))

##########################################################################
####### exact MLE
## lambda estimate
lambda_node = c(0.13705527, 0.09809074, 0.29796933, 0.41240071)


ex_lambda_estimation_vec = estimation_vec(lambda, nod_points, lambda_node, h = 0.09)


##########################################################################
####### exact MLE
## nu estimate
nu_node = c(0.3991226,0.3498865, 0.3632282, 0.3030800)


ex_nu_estimation_vec = estimation_vec(neu, nod_points, nu_node, h = 0.09)

##########################################################################
####### dst MLE
## Sigma estimate
sigma_node = c(0.9641589, 0.4134835, 0.7929231, 0.5432651)
sigma_node = sqrt(sigma_node)

dst_sigma_estimation_vec = estimation_vec(sigma, nod_points, sigma_node, h = 0.09)


#quilt.plot(d1$x,d1$y, abs(sigma_est_vec - sigma_vec))

##########################################################################
####### dst MLE
## lambda estimate
lambda_node = c(0.13707335, 0.09809084, 0.32796309, 0.41240731)


dst_lambda_estimation_vec = estimation_vec(lambda, nod_points, lambda_node, h = 0.09)



##########################################################################
####### dst MLE
## nu estimate
nu_node = c(0.4098864,0.3201226, 0.3832280, 0.3038182)


dst_nu_estimation_vec = estimation_vec(neu, nod_points, nu_node, h = 0.09)

##########################################################################
####### geoR package
## Sigma estimate
sigma_node = c(0.5742, 0.3408, 0.5903, 0.3331)
sigma_node = sqrt(sigma_node)

geor_sigma_estimation_vec = estimation_vec(sigma, nod_points, sigma_node, h = 0.09)


##########################################################################
####### geoR package
## lambda estimate
lambda_node = c(0.4550, -0.0920, 0.6854, 0.3996)


geor_lambda_estimation_vec = estimation_vec(lambda, nod_points, lambda_node, h = 0.09)


##########################################################################





par(mfrow=c(3,3))

# original parameter values 
quilt.plot(d1$x,d1$y, exact_sigma_estimation_vec$original, main = 
             "Parameter values", xlab = "sigma")
quilt.plot(d1$x,d1$y, ex_lambda_estimation_vec$original,xlab = "lambda")
quilt.plot(d1$x,d1$y, ex_nu_estimation_vec$original,xlab = "neu")

# MLE exact estimates using ExaGeoStatR
quilt.plot(d1$x,d1$y, exact_sigma_estimation_vec$estimated, main = 
             "exact MLE estimates", xlab = "sigma")
quilt.plot(d1$x,d1$y, ex_lambda_estimation_vec$estimated, xlab = "lambda")
quilt.plot(d1$x,d1$y, ex_nu_estimation_vec$estimated, xlab = "neu")

# MLE DST estimates using ExaGeoStatR
quilt.plot(d1$x,d1$y, dst_sigma_estimation_vec$estimated,main = 
             "DST MLE estimtaes", xlab = "sigma")
quilt.plot(d1$x,d1$y, dst_lambda_estimation_vec$estimated,xlab = "lambda")
quilt.plot(d1$x,d1$y, dst_nu_estimation_vec$estimated,xlab = "neu")

# MLE estimates using geoR
par(mfrow=c(1,2))
quilt.plot(d1$x,d1$y, geor_sigma_estimation_vec$estimated,main = 
             "geoR estimate", xlab = "sigma")
quilt.plot(d1$x,d1$y, geor_lambda_estimation_vec$estimated, xlab = "lambda")

  
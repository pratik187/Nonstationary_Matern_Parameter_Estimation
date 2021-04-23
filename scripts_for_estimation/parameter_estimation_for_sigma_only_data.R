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


nod_points = list(n1 = c(0.165,0.165), n2 = c(0.165,0.495), n3 = c(0.165,0.825),
                  n4 = c(0.495,0.165),n5 = c(0.495,0.495),n6 = c(0.495,0.825),
                  n7 = c(0.825,0.165),n8 = c(0.825,0.495),n9 = c(0.825,0.825))

param_est = function(u,nod_points,sigma,h){
  k_part = rep(NA,4)
  for(i in 1:9){
    
    k_part[i] = k(u,unlist(nod_points[i],use.names = FALSE),h)
  }
  sum_k = sum(k_part)
  param_est1 = 0
  for(i in 1:9){
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

# function for partial sill
sigma = function(u){
  return((0.33 * exp(-(u[1]+u[2]))) + 0.8)
}  

##########################################################################
####### exact MLE
## Sigma estimate
sigma_node = c(1.3582995, 
               0.4790440 ,
               0.4355629,
               0.7239371 ,
               0.3554572 ,
               0.5093547 ,
               0.9539438 ,
               0.6456552 ,
               0.4004263 )
sigma_node = sqrt(sigma_node)

exact_sigma_estimation_vec = estimation_vec(sigma, nod_points, sigma_node, h = 0.09)


par(mfrow=c(1,2))

# original parameter values 
quilt.plot(d1$x,d1$y, exact_sigma_estimation_vec$original, main = 
             "Parameter values", xlab = "sigma")

# MLE exact estimates using ExaGeoStatR
quilt.plot(d1$x,d1$y, exact_sigma_estimation_vec$estimated, main = 
             "DST MLE estimates", xlab = "sigma")



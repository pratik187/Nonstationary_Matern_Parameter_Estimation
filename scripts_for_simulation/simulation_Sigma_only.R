library('exageostatr')

dmetric = "euclidean"
n = 6400

sigma_sq = 1
beta = 0.1
nu = 0.6
seed = 0

exageostat_init(hardware = list(ncores = 20, ngpus=0,
                                ts = 320, lts = 600, pgrid = 1, qgrid = 1))

data = simulate_data_exact(sigma_sq, beta, nu, dmetric, n, seed)

exageostat_finalize()

data = data.frame(data)

sigma = function(u){
  return((0.33 * exp(-(u[1]+u[2]))) + 0.8)
}

sigma_vec = rep(NA, n)

for(i in 1:n){
  u1<- unlist(data[i,1:2],use.names = FALSE)
  sigma_vec[i] = sigma(u1)
}

data$z = data$z * sigma_vec

setwd("/home/nagp/project/")
write.csv(x = data,file = "simulated_data_sigma.csv")





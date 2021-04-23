library('exageostatr')
print("exageostatr loaded")
data = read.csv('simulated_data_sigma.csv')
#data1 = data.frame(data)
data1 = data[,2:4]
colnames(data1) <- c("x","y","z")

df_xy11 = data1[data1$x >= 0 & data1$x <= 0.33 & data1$y >=0 & data1$y < 0.33,]
df_xy12 = data1[data1$x >= 0 & data1$x <= 0.33 & data1$y > 0.33 & data1$y <= 0.66,]
df_xy13 = data1[data1$x >= 0 & data1$x <= 0.33 & data1$y > 0.66 & data1$y <= 1,]
df_xy21 = data1[data1$x > 0.33 & data1$x <= 0.66 & data1$y >=0 & data1$y < 0.33,]
df_xy22 = data1[data1$x > 0.33 & data1$x <= 0.66 & data1$y > 0.33 & data1$y <= 0.66,]
df_xy23 = data1[data1$x > 0.33 & data1$x <= 0.66 & data1$y > 0.66 & data1$y <= 1,]
df_xy31 = data1[data1$x > 0.66 & data1$x <= 1 & data1$y >=0 & data1$y < 0.33,]
df_xy32 = data1[data1$x > 0.66 & data1$x <= 1 & data1$y > 0.33 & data1$y <= 0.66,]
df_xy33 = data1[data1$x > 0.66 & data1$x <= 1 & data1$y > 0.66 & data1$y <= 1,]
#write.csv(data_xy_11,file = "sample.csv")
print("third split done")
df_list = list(df_xy11,df_xy12,df_xy13,df_xy21,df_xy22,df_xy23,df_xy31,df_xy32,df_xy33)

for(i in 1:length(df_list)){
  # print(length(df_list[i]$X))
  print("i")
  df = data.frame(df_list[i])
  
  dmetric = "euclidean"

  exageostat_init(hardware = list(ncores = 20, ngpus=0,
                                  ts = 320, lts = 600, pgrid = 1, qgrid = 1))

  # data = simulate_data_exact(sigma_sq, beta, nu, dmetric, n, seed)
  result = exact_mle(list(x = df[, 1],
                          y = df[, 2],
                          z = df[, 3]),
                     dmetric,
                     optimization = list(clb = c(0.001,0.001, 0.001),
                                         cub = c(5,5,5),
                                         tol = 1e-5,
                                         max_iters = 2000))

  exageostat_finalize()
  sink(paste0("result", i, ".txt"))
  print(list(result))
  sink()
  
}









data = read.csv('simulated_data_NS_matern.csv')
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
#write.csv(df_xy11,file = "sample.csv")



data = df_xy33

library('exageostatr')

dmetric = "euclidean"
tlr_acc = 2
tlr_maxrank = 450
dst_thick = 6

exageostat_init(hardware = list(ncores = 20, ngpus=0,
                                ts = 320, lts = 0, pgrid = 1, qgrid = 1))

#data = simulate_data_exact(sigma_sq, beta, nu, dmetric, n, seed)
# result = dst_mle(list(x = data[, 2],
#                           y = data[, 3],
#                           z = data[, 4]),
# 			dst_thick,
#                         dmetric,
#                      optimization = list(clb = c(0.001,0.001, 0.001),
#                                          cub = c(5,5,5),
#                                          tol = 1e-3,
#                                          max_iters = 2000))
result = exact_mle(list(x = data[, 2],
                          y = data[, 3],
                          z = data[, 4]),
                     dmetric,
                     optimization = list(clb = c(0.001,0.001, 0.001),
                                         cub = c(5,5,5),
                                         tol = 1e-4,
                                         max_iters = 2000))

exageostat_finalize()

result
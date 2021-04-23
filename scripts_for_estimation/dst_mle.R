library('exageostatr')
print("exageostatr loaded")
data = read.csv('sample.csv')
#data1 = data.frame(data)
data1 = data[,2:4]
#colnames(data1) <- c("x","y","z")

# df_xy11 = data1[data1$x >= 0 & data1$x <= 0.33 & data1$y >=0 & data1$y < 0.33,]
# df_xy12 = data1[data1$x >= 0 & data1$x <= 0.33 & data1$y > 0.33 & data1$y <= 0.66,]
# df_xy13 = data1[data1$x >= 0 & data1$x <= 0.33 & data1$y > 0.66 & data1$y <= 1,]
# df_xy21 = data1[data1$x > 0.33 & data1$x <= 0.66 & data1$y >=0 & data1$y < 0.33,]
# df_xy22 = data1[data1$x > 0.33 & data1$x <= 0.66 & data1$y > 0.33 & data1$y <= 0.66,]
# df_xy23 = data1[data1$x > 0.33 & data1$x <= 0.66 & data1$y > 0.66 & data1$y <= 1,]
# df_xy31 = data1[data1$x > 0.66 & data1$x <= 1 & data1$y >=0 & data1$y < 0.33,]
# df_xy32 = data1[data1$x > 0.66 & data1$x <= 1 & data1$y > 0.33 & data1$y <= 0.66,]
# df_xy33 = data1[data1$x > 0.66 & data1$x <= 1 & data1$y > 0.66 & data1$y <= 1,]
#write.csv(df_xy11,file = "sample.csv")


df_xy11 = data1[data1$x >= 0 & data1$x <= 0.50 & data1$y >=0 & data1$y < 0.50,]
df_xy12 = data1[data1$x >= 0 & data1$x <= 0.50 & data1$y > 0.50 & data1$y <= 1,]
df_xy21 = data1[data1$x > 0.50 & data1$x <= 1 & data1$y >=0 & data1$y < 0.50,]
df_xy22 = data1[data1$x > 0.50 & data1$x <= 1 & data1$y > 0.50 & data1$y <= 1,]
print("third split done")
# df_list = list(df_xy11,df_xy12,df_xy13,df_xy21,df_xy22,df_xy23,df_xy31,df_xy32,df_xy33)
df_list = list(df_xy11,df_xy12,df_xy21,df_xy22)

exageostat_init(hardware = list(ncores = 20, ngpus = 0, ts = 320, lts = 0, 
                                pgrid = 1, qgrid = 1))#Initiate exageostat instance
Exa_local_mle <- matrix(0,4,3)
Exa_local_time<- matrix(0,4,3) 
dst_thick = 3
dmetric = "euclidean"
for (k in 1:4) {
  df = data.frame(df_list[k])
  result = dst_mle(list(x = df[, 1],
                        y = df[, 2],
                        z = df[, 3]),
                   dst_thick,
                   dmetric,
                   optimization = list(clb = c(0.001,0.001, 0.001),
                                       cub = c(5,5,5),
                                       tol = 1e-5,
                                       max_iters = 2000))
  
  Exa_local_mle[k,] = c(result$sigma_sq, result$beta, result$nu)
  Exa_local_time[k,] = c(result$time_per_iter, result$total_time, 
                         result$no_iters)
}

exageostat_finalize()
Exa_local_mle
Exa_local_time


sink(paste0("result_dst_part4.txt"))
print(Exa_local_mle)
print(Exa_local_time)
sink()

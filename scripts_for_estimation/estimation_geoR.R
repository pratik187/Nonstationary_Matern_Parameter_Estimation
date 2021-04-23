library(geoR)
library('exageostatr')
print("exageostatr loaded")
data = read.csv('simulated_data_NS_matern.csv')
#data1 = data.frame(data)
data1 = data[,2:4]
colnames(data1) <- c("x","y","z")


df_xy11 = data1[data1$x >= 0 & data1$x <= 0.50 & data1$y >=0 & data1$y < 0.50,]
df_xy12 = data1[data1$x >= 0 & data1$x <= 0.50 & data1$y > 0.50 & data1$y <= 1,]
df_xy21 = data1[data1$x > 0.50 & data1$x <= 1 & data1$y >=0 & data1$y < 0.50,]
df_xy22 = data1[data1$x > 0.50 & data1$x <= 1 & data1$y > 0.50 & data1$y <= 1,]
print("third split done")
# df_list = list(df_xy11,df_xy12,df_xy13,df_xy21,df_xy22,df_xy23,df_xy31,df_xy32,df_xy33)
df_list = list(df_xy11,df_xy12,df_xy21,df_xy22)


Exa_local_mle <- matrix(0,4,3)
Exa_local_time<- matrix(0,4,3) 
dst_thick = 2
dmetric = "euclidean"
start_time <- Sys.time()
sink(paste0("result_geoR_part4.txt"))

for (k in 1:4) {
  df = data.frame(df_list[k])
  the_data = as.geodata(df,coords.col = 1:2, data.col = 3)
  ml.n <- likfit(the_data , ini = c(0.01,0.01),fix.nugget = TRUE, nugget = 0)
  print("#############################################")
  print(paste0("part ",k," is complete"))
  print(ml.n)
}
print("#############################################")
end_time <- Sys.time()
print("total time taken :")
print(start_time - end_time)
sink()

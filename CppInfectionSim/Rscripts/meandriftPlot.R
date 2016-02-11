library(data.table)
tmp <- fread("~/tmp/outputs/voutput1_1.csv",data.table=FALSE)

means <- NULL
nullWindow <- c(0,10)
for(i in 1:3000){
  print(i)
  window <- nullWindow + i
  means[i] <- mean(tmp[tmp$birth >= window[1] & tmp$birth < window[2],"distRoot"])
}

tmp1 <- fread("~/tmp/outputs/voutput1_3.csv",data.table=FALSE)

means1 <- NULL

for(i in 1:3000){
  print(i)
  window <- nullWindow + i
  means1[i] <- mean(tmp1[tmp1$birth >= window[1] & tmp1$birth < window[2],"distRoot"])
}

tmp2 <- fread("~/tmp/outputs/voutput1_4.csv",data.table=FALSE)

means2 <- NULL
for(i in 1:3000){
  print(i)
  window <- nullWindow + i
  means2[i] <- mean(tmp2[tmp2$birth >= window[1] & tmp2$birth < window[2],"distRoot"])
}

plot(means,type='l',col="blue", ylim=c(0,50))
lines(means1,col="red")
lines(means2,col="green")

fromV <- NULL
fromR <- NULL
for(i in 1:3000){
  print(i)
  fromV[i] <- max(tmp[tmp$birth==i,"changeFromV"])
  fromR[i] <- max(tmp[tmp$birth==i,"changeFromR"])
}
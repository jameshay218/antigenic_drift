slebs <- integer(101)

for(i in 1:length(grebs)){
  slebs[grebs[i]+1] <- slebs[grebs[i]+1] + 1
}

voutput1_2 <- read.csv("~/tmp/outputs/voutput1_2.csv")
tmp <- voutput1_2[sample(nrow(voutput1_2),10000),]
plot(tmp$bindingAvidityFinal~tmp$death)

hostKs_1 <- read.csv("hostKs_2.csv",header=TRUE)
hostK <- as.numeric(hostKs_1[nrow(hostKs_1),1:25])
hostK1 <- NULL
hostK1[1] <- hostK[1] + hostK[2] + hostK[3] + hostK[4]
for(i in 4:27){
  hostK1[i-2] <- hostK[i+1] 
}
plot(hostK,ylim=c(0,600000))
plot(hostK1,ylim=c(0,600000))


voutput1_2 <- read.csv("~/tmp/outputs/voutput1_2.csv")
grebs <- NULL
grebs1 <- NULL
median <- NULL
upper <- NULL
lower <- NULL
for(i in 1:1000){
  grebs[i] <- mean(voutput1_2[voutput1_2$birth==i,"bindingAvidityFinal"])
  median[i] <- median(voutput1_2[voutput1_2$birth==i,"bindingAvidityFinal"])
  tmp <- quantile(voutput1_2[voutput1_2$birth==i,"bindingAvidityFinal"], seq(0,1,0.25))
  lower[i] <- tmp[2]
  upper[i] <- tmp[4]
  grebs1[i] <- mean(voutput1_2[voutput1_2$birth==i,"bindingAvidityIni"])
}
plot(median,type='l',col="red",ylim=c(0,1.5))
lines(grebs,col="blue")
lines(lower)
lines(upper)
abline(h=0.8)
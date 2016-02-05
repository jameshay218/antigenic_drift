tmp <- read.csv("~/tmp/outputs/voutput1_1.csv")

x <- NULL #Time
y <- NULL #Binding avidity
z <- NULL #Mutation size

samps <- tmp[(tmp$death == -1 & tmp$birth > 0),]
i <- sample(nrow(samps),1)

wow <- startVid <- tmp[tmp$vid==samps[i,"vid"],"vid"]

index <- 1
nextVid <- 1
while(nextVid != 0){
  nextVid <- tmp[tmp$vid==startVid,"parentid"]
  x[index] <- tmp[tmp$vid==startVid,"birth"]
  y[index] <- tmp[tmp$vid==startVid,"bindingAvidityFinal"]
  z[index] <- tmp[tmp$vid==startVid,"distRoot"]
  index <- index + 1
  startVid <- nextVid
  print(startVid)
}
dat <- data.frame(time = x, V=y ,delta=z )
dat1 <- dat[with(dat,order(time,delta)),]
dat1 <- dat1[dat1$time > 0,]
plot(dat1$delta~dat1$time,type='l',ylim=c(0,15))
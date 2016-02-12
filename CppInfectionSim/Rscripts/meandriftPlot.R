library(data.table)
tmp <- fread("~/tmp/outputs/firstDriftRun/voutput1_1.csv",data.table=TRUE)
dt <- data.table(tmp)
dt[dt$immJ<0,"immJ"] <- 0
means <- data.frame(dt[,list(mean=mean(distRoot)),by=birth])
lower <- data.frame(dt[,list(lower=quantile(distRoot,c(0.025))),by=birth])
upper <- data.frame(dt[,list(upper=quantile(distRoot,c(0.975))),by=birth])

tmp1 <- fread("~/tmp/outputs/firstDriftRun/voutput1_3.csv",data.table=TRUE)
dt1 <- data.table(tmp1)
dt1[dt1$immJ<0,"immJ"] <- 0
means1 <- data.frame(dt1[,list(mean=mean(distRoot)),by=birth])
lower1 <- data.frame(dt1[,list(lower=quantile(distRoot,c(0.025))),by=birth])
upper1 <- data.frame(dt1[,list(upper=quantile(distRoot,c(0.975))),by=birth])

tmp2 <- fread("~/tmp/outputs/firstDriftRun/voutput1_4.csv",data.table=FALSE)
dt2 <- data.table(tmp2)
dt2[dt1$immJ<0,"immJ"] <- 0
means2 <- data.frame(dt2[,list(mean=mean(distRoot)),by=birth])
lower2 <- data.frame(dt2[,list(lower=quantile(distRoot,c(0.025))),by=birth])
upper2 <- data.frame(dt2[,list(upper=quantile(distRoot,c(0.975))),by=birth])

plot(means,type='l',lwd=1,col="blue", ylim=c(0,40),xlab="Time (days)",ylab="Average Antigenic Distance")
#lines(lower,col="blue",lwd=0.5)
#lines(upper,col="blue",lwd=0.5)
#lines(means1,col="red")
#lines(lower1,col="red",lwd=0.5)
#lines(upper1,col="red",lwd=0.5)
lines(means2,col="green")
#lines(lower2,col="green",lwd=0.5)
#lines(upper2,col="green",lwd=0.5)
#legend(1500,45,c("R only","R + BP","R + no BP"),lty=c(1,1),lwd=c(2,2,2),col=c("blue","red","green"))
legend(2000,40,c("Random Drift Only","Random Drift + V"),lty=c(1,1),lwd=c(2,2),col=c("blue","green"))
#plot(data.frame(dt[,list(mean=mean(distRoot)),by=birth]),type='l')
tmp1 <- fread("~/tmp/outputs/firstDriftRun/voutput1_3.csv",data.table=TRUE)
dt1 <- data.table(tmp1)
dt1[dt1$immJ<0,"immJ"] <- 0
means2 <- data.frame(dt1[,list(mean=mean(bindingAvidityFinal)),by=birth])
tmp <- dt1[dt1$distToParent > 0.5,]
plot(tmp$distToParent~tmp$birth,ylim=c(0,4))
plot(means2,type='l')
# plot host immunity dynamics
# Simulated immunodynamice for single peak with adaptive binding avidity

# Heatmap plot
library(ggplot2)
library(reshape2)
library(reldist)
time <- 1001
maximm <- 12 # the maximum value at Y 

myhome <- "F:/Documents/GitHub/antigenic_drift/bindingavid"
# file for single epidemic fixed binding with low initial
#inputAa <- paste(myhome, "/outputs/single_adaptive_low/.001/hostKs_",1,".csv",sep="")
inputAa <- paste(myhome, "/outputs/single_adaptive_high/.004/hostKs_",1,".csv",sep="")

dat <- read.csv(inputAa,nrows=time+1)
#dat <- read.csv("~/tmp/outputs/hostKs_1.csv", header=FALSE)
dat <- unname(as.matrix(dat))
dat1 <- t(apply(dat,1, function(x) x/sum(x)))
dat1 <- dat1[1:((time/5)+0),1:maximm]

means <- NULL
for(i in 1:nrow(dat1)){
  x1 <- seq(0,ncol(dat1)-1,by=1)
  weights <- dat1[i,]
  #quantiles <- wtd.quantile(x1,probs=c(0.025,0.975),weight=weights)
  tmp <- NULL
  for(j in 1:ncol(dat1)){
    tmp[j] <- dat1[i,j] * j
  }
  means[i] <- sum(tmp)/sum(dat1[i,1:ncol(dat1)])
}
meanDat <- data.frame(y=means,x=seq(0,time-5,by=5))

x <- seq(0,time,by=5)
y <- seq(0,ncol(dat1),by=1)
#image(log(dat+1))
colr <- rainbow(1000,start=4/6,end=0)
#image(x=x,y=y,z=dat1,col=colr)

tmp <- melt(dat1)
colnames(tmp) <- c("Var1","Var2","Proportion")
tmp$Var1 <- (tmp$Var1-1)*5

plot <- ggplot(tmp) + geom_raster(aes(x=Var1,y=Var2,fill=Proportion),interpolate=FALSE) +
  scale_fill_gradientn(colours=c("darkblue","dodgerblue","orange","red")) +
  scale_y_continuous(expand=c(0,0),breaks=seq(0,max(tmp$Var2),by=1),labels=seq(-1,max(tmp$Var2)-1,by=1))+
  #scale_y_continuous(expand=c(0,0),breaks=seq(0,10,by=1),labels=seq(-1,10-1,by=1))+
  
  scale_x_continuous(expand=c(0,0),breaks=seq(0,time,by=100))+
  xlab("Time (days)")+
  ylab("Host Immunity (J)")+
  geom_line(data=meanDat,aes(y=y,x=x),linetype=2,colour="white",size=1) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor=element_blank(),
      panel.background = element_blank(),
      text=element_text(size=16,colour="gray20"),
      axis.line=element_line(colour="gray20"),
      axis.line.x = element_line(colour = "gray20"),
      axis.line.y=element_line(colour="gray20")
      )
print(plot)
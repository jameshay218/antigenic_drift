# plot infecteds dynamics
# Simulated epidemic for single peak
# Can we also calculate the maximum peak?

library(graphics)
library(ggplot2)

runs <- 5
duration <- 750
max_xlim <- 700
duration <- max_xlim

allDatAa <- matrix(nrow=runs,ncol=duration+1)
allDatBa <- matrix(nrow=runs,ncol=duration+1)

allDatAb<- matrix(nrow=runs,ncol=duration+1)
allDatBb <- matrix(nrow=runs,ncol=duration+1)
upperAa <- NULL
lowerAa <- NULL
meanAa <- NULL

upperBa <- NULL
lowerBa<- NULL
meanBa <- NULL

upperAb <- NULL
lowerAb <- NULL
meanAb <- NULL

upperBb <- NULL
lowerBb <- NULL
meanBb <- NULL


for(i in 1:runs){
    #
    print(i)
    IBb <- numeric()
    myhome <- "E:/James/Documents/antigenic_drift/bindingavid"
    # file for single epidemic adaptive binding with high initial binding
    inputBb <- paste(myhome, "/outputs/single_adaptive_high/.001/voutput1_",i,".csv",sep="")
    
    #datBb <- read.csv(inputBb,header=TRUE,nrows=duration+1)
    datBb <- data.table::fread(inputBb,header=TRUE,data.table=FALSE)
    datBbBirth <- datBb[,c("birth","bindingAvidityIni","bindingAvidityFinal")]
    #transform into date and value
    for (day in 1:max_xlim+1) {
      #apply(dat,1,function(x) mean(na.omit(x)))
      if (day < max(datBbBirth)) {
      #meandb <- apply(datBbBirth[datBbBirth$birth == day, 'bindingAvidityIni'],1,function(x) mean(na.omit(x)))
      
      meanbd <- mean(datBbBirth[datBbBirth$birth == day, 'bindingAvidityIni'])
      if (is.nan(meanbd)) {
        print(meanbd)
        meanbd <- IBb[day-1]
        print(meanbd)
      }
      #max_xlim
      IBb[day] <- meanbd
      }
      
      if (day >= max(datBbBirth)) {
        #IBb[day] <- mean(IBb[max(datBbBirth)-2])
        IBb[day] <- 0
      }
      #print(meanbd)
    }
    
    #
    allDatBb[i,1:length(IBb)] <- IBb 
}

allDataBb <- allDatBb[,1:duration+1]

for(i in 1:ncol(allDatBb)){
    # Adaptive high
    tmp <- quantile(allDatBb[,i],probs=c(0.025,0.5,0.975), na.rm=TRUE)
    upperBb[i] <- tmp[3]
    lowerBb[i] <- tmp[1]
    meanBb[i] <- tmp[2]
}

CI.x.topBb <- seq(0,duration,by=1)
CI.x.botBb <- rev(seq(0,duration,by=1))
CI.x.Bb <- c(CI.x.topBb,CI.x.botBb)
CI.y.topBb <- upperBb
CI.y.botBb <- rev(lowerBb)
CI.y.Bb <- c(CI.y.topBb,CI.y.botBb)
CI.colBb <- adjustcolor("red",alpha.f=0.25)




#png("low.png")
#plot(meanAa,ylim=c(0,50000),type='l',col="blue")
#lines(meanBa,col="red")
#polygon(CI.x.Aa,CI.y.Aa,col=CI.colAa,border=NA)
#polygon(CI.x.Ba,CI.y.Ba,col=CI.colBa,border=NA)
#dev.off()

#png("high.png")
#plot(meanAb,ylim=c(0,50000),type='l',col="blue")
#lines(meanBb,col="red")

#polygon(CI.x.Ab,CI.y.Ab,col=CI.colAb,border=NA)
#polygon(CI.x.Bb,CI.y.Bb,col=CI.colBb,border=NA)
#dev.off()

meanB <- data.frame(x=seq(1,length(meanBb),by=1),y=meanBb,Scenario="Binding Avidity")

allDat <- rbind(meanB)

boundB <- data.frame(x=CI.x.Bb,y=CI.y.Bb,Scenario="Binding Avidity")

allBounds <- rbind(boundB)

#allBounds[,1] <- as.numeric(allBounds[,1])
                                        #allBounds[,2] <- as.numeric(allBounds[,2])
#xaxislim <- 1200
ggplot() +
  geom_polygon(data=allBounds,aes(x=x,y=y),alpha=0.3) +  
  geom_line(data=allDat,aes(x=x,y=y,fill=Scenario,colour=Scenario),size=0.75) +
    
    xlab("Time (days)") +
    scale_y_continuous(limits=c(0,1),expand=c(0,0))+
    scale_x_continuous(limits=c(0,duration),expand=c(0,0))+
    ylab("Binding Avidity") +
    #xlim(0, xaxislim) +
    theme(
        #legend.justification=c(1,1),
        #legend.position=c(1,1),
        text=element_text(size=12,colour="gray20"),
        legend.text=element_text(size=10,colour="gray20"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line=element_line(colour="gray20"),
        axis.line.x = element_line(colour = "gray20"),
        axis.line.y=element_line(colour="gray20"),
        plot.margin=unit(c(0.5,1,0.5,0.5),"cm"),
        panel.background=element_blank(),
        legend.title=element_blank(),
        axis.text.y=element_text(colour="gray20",size=10))
#print(plot1)
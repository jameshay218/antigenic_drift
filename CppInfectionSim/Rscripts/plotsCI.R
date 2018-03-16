# plot infecteds dynamics
# Simulated epidemic for single peak
# Can we also calculate the maximum peak?

library(graphics)
library(ggplot2)

runs <- 20
duration <- 2000
max_xlim <- 1000
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
    #myhome <- "F:/Documents/GitHub/antigenic_drift/bindingavid"
    myhome <- "D:/Dropbox (NHRI Taiwan)/Workspace/GitHub/antigenic_drift/bindingavid"  
    # file for single epidemic fixed binding with low initial
    inputAa <- paste(myhome, "/outputs/single_fixed_low/.005/scenario_SIR_",i,".csv",sep="")
    # file for single epidemic fixed binding with high initial
    inputAb <- paste(myhome, "/outputs/single_fixed_high/.007/scenario_SIR_",i,".csv",sep="")
    # file for single epidemic adaptive binding with high initial binding
    inputBa <- paste(myhome, "/outputs/single_adaptive_high/.004/scenario_SIR_",i,".csv",sep="")
    # file for single epidemic adaptive binding with low initial binding
    #inputBa <- paste(myhome, "/outputs/single_adaptive_high/.001/scenario_SIR_",i,".csv",sep="")
    

    datAa <- read.csv(inputAa,header=FALSE,nrows=duration+1)
    datBa <- read.csv(inputBa,header=FALSE,nrows=duration+1)
    datAb <- read.csv(inputAb,header=FALSE,nrows=duration+1)
    #datBb <- read.csv(inputBa,header=FALSE)

    
    IAa <- datAa[,2]
    IBa <- datBa[,2]
    
    IAb <- datAb[,2]
    #IBb <- datBb[,2]
    
    allDatAa[i,] <- IAa
    allDatBa[i,] <- IBa
    allDatAb[i,] <- IAb
    #allDatBb[i,] <- IBb    
}


allDataAa <- allDatAa[,1:duration+1]
allDataBa <- allDatBa[,1:duration+1]
allDataAb <- allDatAb[,1:duration+1]

for(i in 1:ncol(allDatAa)){
    # Fixed low
    tmp <- quantile(allDatAa[,i],probs=c(0.025,0.5,0.975))
    upperAa[i] <- tmp[3]
    lowerAa[i] <- tmp[1]
    meanAa[i] <- tmp[2]

    # Adaptive high
    tmp <- quantile(allDatBa[,i],probs=c(0.025,0.5,0.975))
    upperBa[i] <- tmp[3]
    lowerBa[i] <- tmp[1]
    meanBa[i] <- tmp[2]

    #Fixed high
    tmp <- quantile(allDatAb[,i],probs=c(0.025,0.5,0.975))
    upperAb[i] <- tmp[3]
    lowerAb[i] <- tmp[1]
    meanAb[i] <- tmp[2]

    # Adaptive high
    #tmp <- quantile(allDatBb[,i],probs=c(0.025,0.5,0.975))
    #upperBb[i] <- tmp[3]
    #lowerBb[i] <- tmp[1]
    #meanBb[i] <- tmp[2]

    
}


CI.x.topAa <- seq(0,duration,by=1)
CI.x.botAa <- rev(seq(0,duration,by=1))
CI.x.Aa <- c(CI.x.topAa,CI.x.botAa)
CI.y.topAa <- upperAa
CI.y.botAa <- rev(lowerAa)
CI.y.Aa <- c(CI.y.topAa,CI.y.botAa)
CI.colAa <- adjustcolor("blue",alpha.f=0.25)


CI.x.topBa <- seq(0,duration,by=1)
CI.x.botBa <- rev(seq(0,duration,by=1))
CI.x.Ba <- c(CI.x.topBa,CI.x.botBa)
CI.y.topBa <- upperBa
CI.y.botBa <- rev(lowerBa)
CI.y.Ba <- c(CI.y.topBa,CI.y.botBa)
CI.colBa <- adjustcolor("red",alpha.f=0.25)


CI.x.topAb <- seq(0,duration,by=1)
CI.x.botAb <- rev(seq(0,duration,by=1))
CI.x.Ab <- c(CI.x.topAb,CI.x.botAb)
CI.y.topAb <- upperAb
CI.y.botAb <- rev(lowerAb)
CI.y.Ab <- c(CI.y.topAb,CI.y.botAb)
CI.colAb <- adjustcolor("blue",alpha.f=0.25)


#CI.x.topBb <- seq(0,1000,by=1)
#CI.x.botBb <- rev(seq(0,1000,by=1))
#CI.x.Bb <- c(CI.x.topBb,CI.x.botBb)
#CI.y.topBb <- upperBb
#CI.y.botBb <- rev(lowerBb)
#CI.y.Bb <- c(CI.y.topBb,CI.y.botBb)
#CI.colBb <- adjustcolor("red",alpha.f=0.25)



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
polygon(CI.x.Bb,CI.y.Bb,col=CI.colBb,border=NA)
#dev.off()

meanA <- data.frame(x=seq(1,length(meanBa),by=1),y=meanBa,Scenario="Adaptive Binding Avidity (V0=0.7)")
meanFa <- data.frame(x=seq(1,length(meanAa),by=1),y=meanAa,Scenario="Low Fixed V=0.45")
meanFb <- data.frame(x=seq(1,length(meanAb),by=1),y=meanAb,Scenario="High Fixed V=0.7")

allDat <- rbind(meanA,meanFa,meanFb)


boundA <- data.frame(x=CI.x.Ba,y=CI.y.Ba,Scenario="Adaptive Binding Avidity (V0=0.7)")
boundFa <- data.frame(x=CI.x.Aa,y=CI.y.Aa,Scenario="Low Fixed V=0.45")
boundFb <- data.frame(x=CI.x.Ab,y=CI.y.Ab,Scenario="High Fixed V=0.7")

allBounds <- rbind(boundA,boundFa,boundFb)
#allBounds[,1] <- as.numeric(allBounds[,1])
                                        #allBounds[,2] <- as.numeric(allBounds[,2])
#xaxislim <- 1200
plot1 <- ggplot() +
    geom_line(data=allDat,aes(x=x,y=y,fill=Scenario,colour=Scenario),size=0.75) +
    geom_polygon(data=allBounds,aes(x=x,y=y,fill=Scenario),alpha=0.3) +
    xlab("Time (days)") +
    scale_y_continuous(limits=c(0,25000),expand=c(0,0))+
    scale_x_continuous(limits=c(0,duration),expand=c(0,0))+
    ylab("Incidence (individuals)") +
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
print(plot1)
# plot infecteds dynamics
# Simulated epidemic for single peak
# Can we also calculate the maximum peak?
library("plyr")
library(ggplot2)
library(reshape2)
library(graphics)
library(ggplot2)
library(ggthemes)
library(ggpubr)

runs <- 200
duration <- 1000
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
    myhome <- "E:/James/Google Drive/Binding.Yuan.Shared/figures/fig5/fig5data"
  #myhome <- "E:/James/Documents/tmp1/out"
    #myhome <- "H:/all_results_backup/binding_avidity/may2018_analysis/fig5data"
    #myhome <- "D:/Dropbox (NHRI Taiwan)/Workspace/GitHub/antigenic_drift/bindingavid"  
    # file for single epidemic fixed binding with low initial
  #inputAa <- paste(myhome, "/SIRfixed_low_",i,".csv",sep="")
  #inputAb <- paste(myhome, "/SIRfixed_high_",i,".csv",sep="")
  #inputBa <- paste(myhome, "/SIRadaptive_high_",i,".csv",sep="")
    inputAa <- paste(myhome, "/sfl/.005/scenario_SIR_",i,".csv",sep="")
    # file for single epidemic fixed binding with high initial
    inputAb <- paste(myhome, "/sfh/.007/scenario_SIR_",i,".csv",sep="")
    # file for single epidemic adaptive binding with high initial binding
    inputBa <- paste(myhome, "/sah/.004/scenario_SIR_",i,".csv",sep="")
    # file for single epidemic adaptive binding with low initial binding
    #inputBa <- paste(myhome, "/outputs/single_adaptive_high/.001/scenario_SIR_",i,".csv",sep="")
    

    datAa <- read.csv(inputAa,header=FALSE,nrows=duration+1)
    datBa <- read.csv(inputBa,header=FALSE,nrows=duration+1)
    datAb <- read.csv(inputAb,header=FALSE,nrows=duration+1)
    #datBb <- read.csv(inputBa,header=FALSE)

    NAa <- rowSums(datAa)
    NBa <- rowSums(datBa)
    NAb <- rowSums(datAb)
    
    IAa <- cumsum(datAa[,2])/NAa
    IBa <- cumsum(datBa[,2])/NBa
    IAb <- cumsum(datAb[,2])/NAb
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


CI.x.topBb <- seq(0,1000,by=1)
CI.x.botBb <- rev(seq(0,1000,by=1))
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

meanA <- data.frame(x=seq(1,length(meanBa),by=1),y=meanBa,Scenario="Adaptive binding\n avidity (V0=0.7)")
meanFa <- data.frame(x=seq(1,length(meanAa),by=1),y=meanAa,Scenario="Low fixed V0=0.45")
meanFb <- data.frame(x=seq(1,length(meanAb),by=1),y=meanAb,Scenario="High fixed V0=0.7")

allDat <- rbind(meanA,meanFa,meanFb)


boundA <- data.frame(x=CI.x.Ba,y=CI.y.Ba,Scenario="Adaptive binding\n avidity (V0=0.7)")
boundFa <- data.frame(x=CI.x.Aa,y=CI.y.Aa,Scenario="Low fixed V0=0.45")
boundFb <- data.frame(x=CI.x.Ab,y=CI.y.Ab,Scenario="High fixed V0=0.7")

allBounds <- rbind(boundA,boundFa,boundFb)
#allBounds[,1] <- as.numeric(allBounds[,1])
                                        #allBounds[,2] <- as.numeric(allBounds[,2])
#xaxislim <- 1200
mycols <- c("#0072B2","#009E73","#D55E00")
options(scipen=999)
plot1 <- ggplot() +
    geom_line(data=allDat,aes(x=x,y=y,fill=Scenario,colour=Scenario),size=0.75) +
    geom_polygon(data=allBounds,aes(x=x,y=y,fill=Scenario),alpha=0.3) +
    xlab("Time (days)") +
  #coord_cartesian(ylim=c(1,100000)) +
  #scale_y_log10(breaks=c(1,10,100,1000,10000,100000),labels=c(1,10,100,1000,10000,100000)) +
    #scale_y_continuous(limits=c(0,25000),expand=c(0,0))+
    #scale_x_continuous(limits=c(0,duration),expand=c(0,0))+
    ylab("Cumulative incidence per capita") +
  scale_y_continuous(breaks=seq(0,3,by=0.5),limits=c(0,3),expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0,1000,by=100),limits=c(0,1000),expand=c(0,0)) +
  scale_fill_manual(values=mycols) +
  scale_color_manual(values=mycols) +
    theme_pubr() +
    #xlim(0, xaxislim) +
  theme(axis.title=element_text(size=8,face="bold"),
        #legend.position=c(0.9,0.7),
        axis.title.x=element_blank(),
        legend.position="none",
        legend.title.align = 0.5,
        legend.justification = "center",
        legend.text = element_text(size=6),
        legend.background = element_rect(colour="grey70"),
        text = element_text(size=8))#,
        #panel.grid.major=element_line(colour="black",size=0.1))

hostKs <- read.csv("E:/James/Documents/tmp1/hostKs_1.csv",header=FALSE)
hostKs_propn <- hostKs/rowSums(hostKs)
titres <- seq_len(ncol(hostKs))-1
mean_titre <- apply(hostKs_propn,1, function(x) sum(x*titres))
#hostKs_propn <- hostKs_propn[,1:12]
colnames(hostKs_propn) <- titres
hostKs_propn$time <- seq(0,1000,by=5)
hostKs_melted <- reshape2::melt(hostKs_propn,id.vars="time")
hostKs_melted$variable <- as.numeric(as.character(hostKs_melted$variable))


high_immunity <- hostKs_melted[hostKs_melted$variable >= 11,]
high_immunity <- ddply(high_immunity,~time, function(x) sum(x$value))
high_immunity$variable <- ">10"
colnames(high_immunity) <- c("time","value","variable")

low_immunity <- hostKs_melted[hostKs_melted$variable < 11,]
hostKs_melted <- rbind(low_immunity, high_immunity[,c("time","variable","value")])
hostKs_melted$variable <- factor(hostKs_melted$variable, 
                                    levels=c("0","1","2","3","4","5","6","7","8","9","10",">10"))
mean_titres <- data.frame(time=seq(0,1000,by=5),mean=mean_titre)

p2 <- ggplot(hostKs_melted) + 
  geom_raster(aes(x=time,y=variable,fill=value)) +
  geom_line(data=mean_titres,aes(x=time,y=mean),linetype=2,colour="white",size=1) +
  scale_y_discrete(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0,1000,by=100),limits=c(0,1000),expand=c(0,0)) +
  scale_fill_gradientn(colours=c("darkblue","dodgerblue","orange","red")) +
  theme_pubr() +
  xlab("Time (days)") +
  ylab("Host immunity (J)") +
  theme(axis.title=element_text(size=8,face="bold"),
        #legend.position=c(0.9,0.7),
        legend.position="none",
        legend.title.align = 0.5,
        legend.justification = "center",
        legend.text = element_text(size=6),
        legend.background = element_rect(colour="grey70"),
        text = element_text(size=8))

overall_p <- plot1 + p2 + plot_layout(ncol=1,heights=c(1,1))

pdf("fig5.pdf",height=5,width=5)
print(overall_p)
dev.off()


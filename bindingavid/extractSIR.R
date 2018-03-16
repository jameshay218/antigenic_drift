

# read all files in a folder
#Kc = seq(0, 0.5, by=0.05) # check job setting
#path = "outputs/multiple_adaptive_test2/.013/03"


#Figures in Manuscript

#One figure used as manuscript for single season
path = "outputs/single_adaptive_high/.002"  
maxdate <- 1500
ylimmax <- 1500
SIR_multiple_flag <- FALSE #' Flag to plot SIR dynamics for single peak

#One figure used as manuscript for multiple seasons
#path = "outputs/multiple_adaptive_test2/.017" 
#maxdate <- 10000 # check job setting
#ylimmax <- 7000
#SIR_multiple_flag <- TRUE #' Flag to plot SIR dynamics for multiple peaks

Kc = seq(0, 0.6, by=0.05) # check job setting

filelist <- list.files(path,pattern="scenario.*.csv")

# read the index
string <- gsub(".csv", "", filelist)
file_id <- gsub("s.*SIR_", "", string)


extinctdate <- 0
kcval <- 0 # Kc value
Inc_max <- 0 # max Incidence

interval = 200 # check job setting
j <- 0
# read the csv file
for(i in 1:length(file_id)) {
  print(filelist[i])
  filename <- paste("scenario_SIR_", i, ".csv", sep="")
  if (i %% interval == 1) {
    j<-j+1
  }
  
  #read sir numbers from the csv file
  sir <- read.csv(file=paste(path,"/",filename,sep=""))
  
  names(sir) <- cbind("s","i","r")
  Inc_ext <- sir$i==0
  Inc_t <- sir$i
  
  # calculate the extinction time 
  stopdate <- min(which(Inc_ext == TRUE))
  if (is.infinite(stopdate)) {
    stopdate <- maxdate
  }
  if (stopdate > maxdate) {
    stopdate <- maxdate
  }
  extinctdate[i] = stopdate
  kcval[i] = Kc[j]
  
  Inc_max[i] <- max(Inc_t)
}

#create dataframe with columns extincdate and kc
x <- as.data.frame(list(extinctdate=extinctdate, kc=kcval))
### S3 method for class 'formula'
#aggregate(formula, data, FUN, ...,
#          subset, na.action = na.omit)
ext_kc <- aggregate(extinctdate ~ kc,data=x,sum)

#plot(ext_kc$kc, ext_kc$extinctdate)

survtable <- as.data.frame(list(stop=extinctdate, kc = kcval, inc_max=Inc_max)) #stop=extinctdate

# Boxplot of extinction time 
p <- ggplot(survtable, aes(x = kc*0.5, y = stop, group=kc)) + 
  geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=16, size=2,show_guide = FALSE) +
  #geom_bar() +
  coord_flip() + 
  xlab("Rates of binding avidity changes (Kc)") +
  ylab("Extinction time (Days)") +
  ylim(0, ylimmax)
p

# Boxplot of max incidence
# Used for single peak analysis
if (SIR_multiple_flag) {
p2 <- ggplot(survtable, aes(x = kc, y = inc_max/5000, group=kc)) + 
  geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=16, size=2,show_guide = FALSE) +
  xlab("Rates of binding avidity changes") +
  ylab("Max Incidence(%)")
p2
}


#boxplot(stop~kc,data=survtable, main="", 
#        xlab="kc ( The rate of binding avidity adatation)", ylab="Extinction date")





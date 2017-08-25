

# read all files in a folder
#path = "outputs/multiple_adaptive_test2/.013/02"
#path = "outputs/multiple_adaptive_test1/test"
#path = "outputs/multiple_adaptive/test1_highdrift/20170801"
#path = "outputs/multiple_adaptive_test2/.013/02"
#Kc = seq(0, 0.5, by=0.05) # check job setting
#path = "outputs/multiple_adaptive_test2/.013/03"
path = "outputs/multiple_adaptive_test2/.017"  # Used as manuscript for multiple seasons
Kc = seq(0, 0.6, by=0.05) # check job setting

filelist <- list.files(path,pattern="scenario.*.csv")

# read the index
string <- gsub(".csv", "", filelist)
file_id <- gsub("s.*SIR_", "", string)

maxdate <- 10000 # check job setting
extinctdate <- 0
kcval <- 0 # Kc value


interval = 200 # check job setting
j <- 0
# read the csv file
for(i in 1:length(file_id)) {
  print(filelist[i])
  filename <- paste("scenario_SIR_", i, ".csv", sep="")
  if (i %% interval == 1) {
    j<-j+1
  }
  
  sir <- read.csv(file=paste(path,"/",filename,sep=""))
  
  names(sir) <- cbind("s","i","r")
  lv <- sir$i==0
  stopdate <- min(which(lv == TRUE))

  
  if (is.infinite(stopdate)) {
    stopdate <- maxdate
  }
  
  if (stopdate > maxdate) {
    stopdate <- maxdate
  }
# calculate the extinction time 
  extinctdate[i] = stopdate
  kcval[i] = Kc[j]
}

x <- as.data.frame(list(extinctdate=extinctdate, kc=kcval))
ext_kc <- aggregate(extinctdate ~ kc,data=x,sum)
#plot(ext_kc$kc, ext_kc$extinctdate)

survtable <- as.data.frame(list(stop=extinctdate, kc = kcval)) #stop=extinctdate

# Boxplot of MPG by Car Cylinders 
p <- ggplot(survtable, aes(x = kc, y = stop, group=kc)) + 
  geom_boxplot() +
  #geom_bar() +
  coord_flip() + 
  xlab("Rates of binding avidity changes") +
  ylab("Extinction time (Days)") +
  ylim(0, 7000)
p
#boxplot(stop~kc,data=survtable, main="", 
#        xlab="kc ( The rate of binding avidity adatation)", ylab="Extinction date")





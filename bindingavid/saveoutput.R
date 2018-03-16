## Produces output files
#' outputdir: the rtarget folder of output
#' runName: string for the name of the job (outputs files will be saved with this name)
#' input_pars: a named vector of all the parameters needed for the simulation, as in input_params.csv"

saveoutput <- function(outputdir, runName, input_pars){
####################
    ## Setup some pars
####################

  print("create subfolders")
  outdir <- paste("outputs/",runName,sep="")
  subdir = list.dirs(path = outdir, full.names = TRUE, recursive = TRUE)
  subfolders <- sapply(strsplit(subdir, "/"), "[", 3)
  subfolders[is.na(subfolders)]<-0
  v <- as.numeric(subfolders)
  maxv <- max(v*1000)  ## to do 20170620
  nextv = maxv + 1
  newfolder <- ""
  if (nextv >=10) {
    if (nextv >=100) {
      newfolder = paste(".", toString(nextv),sep = "")
    } else {
      newfolder = paste(".0", toString(nextv),sep = "")  
    }
  } else {
    newfolder = paste(".00", toString(nextv),sep = "")  
  }
  newpath <- paste(outdir,"/",toString(newfolder),sep="")
  dir.create(newpath)
  print(newpath)
  if (file.exists(newpath)) {
    filelist <- dir(outdir)
    for (fileID in 1:length(filelist)) {
      #print(filelist[fileID])
      file.rename(from=paste(outdir,"/",filelist[fileID],sep=""),to=paste(newpath,"/",filelist[fileID],sep=""))
    }
  }
  write.csv(input_pars, file = paste(newpath,"/","parameters_test.csv",sep=""), row.names=FALSE)
  return(TRUE)
}


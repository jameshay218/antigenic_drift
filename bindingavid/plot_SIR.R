plot_SIR <- function(filename,N){
            #N <- isolate(inputs$s0) + isolate(inputs$i0) +  isolate(inputs$r0) + 500
            dat <- read.csv(filename,header=0)
            dat <- cbind(seq(1,nrow(dat),by=1),dat)
            #dur <- length(dat[,1]);
            dur <- 400
            colnames(dat) <- c("t","S","I","R")
            dat <- melt(dat,id="t")
            colnames(dat) <- c("t","Population","value")
            SIR_plot <- ggplot() +
                geom_line(data=dat,aes(x=t,y=value,colour=Population,group=Population)) +
                xlab("Time (days)") +
                ylab("Number Individuals") +
                scale_y_continuous(limits=c(0,N),expand=c(0,0))+
                #scale_x_continuous(limits=c(0,inputs$dur+1),expand=c(0,0))+
                scale_x_continuous(limits=c(0,dur+1),expand=c(0,0))+
                theme(
                    text=element_text(colour="gray20",size=14),
                    plot.title=element_text(size=28),
                    legend.text=element_text(size=20,colour="gray20"),
                    legend.title=element_text(size=20,colour="gray20"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x = element_line(colour = "gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    axis.text.x=element_text(colour="gray20"),
                    panel.background=element_blank(),
                    axis.text.y=element_text(colour="gray20"))
            return(SIR_plot)
        }

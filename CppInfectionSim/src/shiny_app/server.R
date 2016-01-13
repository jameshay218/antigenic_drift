library(shiny)
library(ggplot2)
library(reshape2)
library(Rcpp)
library(gridExtra)
source("~/Documents/antigenic_drift/antigenic_drift/CppInfectionSim/src/antigenic_drift.R")

setup_env()

shinyServer(
    function(inputs, output, session){
        
        run_sim <- observeEvent(inputs$run, {
            print("Running simulation")
            #' Make sure Rcpp will compile and run
               
            SIR_flag <- 1 %in% inputs$flags
            voutput1_flag <- 2 %in% inputs$flags
            voutput2_flag <- 3 %in% inputs$flags
            time_flag <- 4 %in% inputs$flags
            VERBOSE <- 5 %in% inputs$flags

            flags <- c(SIR_flag, voutput1_flag, voutput2_flag, time_flag)
            flags <- as.numeric(flags)
            print(flags)
            
            S0 <- as.numeric(inputs$s0)
            I0 <- as.numeric(inputs$i0)
            R0 <- as.numeric(inputs$r0)

            contactRate <- as.numeric(inputs$contact)
            mu <- 1/(as.numeric(inputs$mu)*365)
            wane <- 1/as.numeric(inputs$wane)
            gamma <- 1/as.numeric(inputs$gamma)
            iniBind <- as.numeric(inputs$iniBinding)

            hostpars <- c(S0,I0, R0,contactRate,mu,wane,gamma,iniBind)

            p = as.numeric(inputs$p) 
            r = as.numeric(inputs$r)
            q = as.numeric(inputs$q)
            a =as.numeric(inputs$a)
            b = as.numeric(inputs$b)
            n = as.numeric(inputs$n)
            v = as.numeric(inputs$v)

            probMut = inputs$probMut
            expDist = inputs$expDist
            kc = inputs$kc
            VtoD = inputs$VtoD

            viruspars <- c(p,r,q,a,b,n,v,probMut,expDist,kc,VtoD)
            print("Host pars:")
            print(hostpars)
            print("Virus pars:")
            print(viruspars)
            if(length(inputs$scenarios) > 0){
                for(i in inputs$scenarios){
                    filename1 <- paste("scenario_",i,"_SIR.csv",sep="")
                    filename2 <- paste("voutput1_",i,".csv",sep="")
                    filename3 <- paste("voutput2_",i,".csv",sep="")
                    filenames <- c(filename1, filename2, filename3)
                    run_simulation(flags,hostpars,viruspars,0,inputs$dur,filenames,VERBOSE, i)                }
            }
            
            
        })

        output$sim_main_1<- renderPlot({
            inputs$run
            if(1 %in% inputs$scenarios){
                N <- isolate(inputs$s0) + isolate(inputs$i0) +  isolate(inputs$r0)
                dat <- read.csv("scenario_1_SIR.csv",header=0)
                dat <- cbind(seq(1,nrow(dat),by=1),dat)
                colnames(dat) <- c("t","S","I","R")
                dat <- melt(dat,id="t")
                colnames(dat) <- c("t","Population","value")
                SIR_plot <- ggplot() +
                    geom_line(data=dat,aes(x=t,y=value,colour=Population,group=Population)) +
                    xlab("Time (days)") +
                    ylab("Number Individuals") +
                    scale_y_continuous(limits=c(0,N),expand=c(0,0))+
                    scale_x_continuous(limits=c(0,inputs$dur+1),expand=c(0,0))+
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
                SIR_plot
            }
            else{
                NULL
            }
        })

         output$sim_main_3<- renderPlot({
                         inputs$run
            if(3 %in% inputs$scenarios){
                N <- isolate(inputs$s0) +isolate(inputs$i0) +  isolate(inputs$r0)
                dat <- read.csv("scenario_3_SIR.csv",header=0)
                dat <- cbind(seq(1,nrow(dat),by=1),dat)
                colnames(dat) <- c("t","S","I","R")
                dat <- melt(dat,id="t")
                colnames(dat) <- c("t","Population","value")
                SIR_plot <- ggplot() +
                    geom_line(data=dat,aes(x=t,y=value,colour=Population,group=Population)) +
                    xlab("Time (days)") +
                    ylab("Number Individuals") +
                    scale_y_continuous(limits=c(0,N),expand=c(0,0))+
                    scale_x_continuous(limits=c(0,inputs$dur+1),expand=c(0,0))+
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
                SIR_plot
            }
            else{
                NULL
            }
        })
        output$sim_main_4<- renderPlot({
                        inputs$run
            if(4 %in% inputs$scenarios){
                N <- isolate(inputs$s0) +isolate(inputs$i0) +  isolate(inputs$r0)
                dat <- read.csv("scenario_4_SIR.csv",header=0)
                dat <- cbind(seq(1,nrow(dat),by=1),dat)
                colnames(dat) <- c("t","S","I","R")
                dat <- melt(dat,id="t")
                colnames(dat) <- c("t","Population","value")
                SIR_plot <- ggplot() +
                    geom_line(data=dat,aes(x=t,y=value,colour=Population,group=Population)) +
                    xlab("Time (days)") +
                    ylab("Number Individuals") +
                    scale_y_continuous(limits=c(0,N),expand=c(0,0))+
                    scale_x_continuous(limits=c(0,inputs$dur+1),expand=c(0,0))+
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
                SIR_plot
            }
            else{
                NULL
            }
        })
        output$sim_main_2<- renderPlot({
                        inputs$run
            if(2 %in% inputs$scenarios){
                N <- isolate(inputs$s0) +isolate(inputs$i0) +  isolate(inputs$r0)
                dat <- read.csv("scenario_2_SIR.csv",header=0)
                dat <- cbind(seq(1,nrow(dat),by=1),dat)
                colnames(dat) <- c("t","S","I","R")
                dat <- melt(dat,id="t")
                colnames(dat) <- c("t","Population","value")
                SIR_plot <- ggplot() +
                    geom_line(data=dat,aes(x=t,y=value,colour=Population,group=Population)) +
                    xlab("Time (days)") +
                    ylab("Number Individuals") +
                    scale_y_continuous(limits=c(0,N),expand=c(0,0))+
                    scale_x_continuous(limits=c(0,inputs$dur+1),expand=c(0,0))+
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
                SIR_plot
            }
            else{
                NULL
            }
        })
        
        output$Main <- renderPlot({
            p = as.numeric(inputs$p) #' parameter to control degree by which changes in binding avidity affect probability of escape from immune response
            r = as.numeric(inputs$r) #' parameter to control degree by which previous exposure reduce probability of immune escape
            b = as.numeric(inputs$b) #' parameter to control the shape of the relationship between probability of successful replication and changes in binding avidity
            a =as.numeric(inputs$a) #' controls rate of changes of relationship between probability of successful replication and change in binding avidity.
            c = as.numeric(inputs$c) #' per day contact rate
            n = as.numeric(inputs$n) #' average number of virus copies: n is number of offspring per virus replication
            v = as.numeric(inputs$v) #' v is number of virions initialyl transmitted
            nv = n*v
            q = as.numeric(inputs$q) #' parameter to control the shape of the relationship between binding avidity and immune escape (shift on the x-axis)
            V = seq(0, 2, by = 0.005) #' Bidning avidity
            N_reinfect = as.numeric(inputs$N_reinfect)
            max_reinfect = N_reinfect
         
            delta = as.numeric(inputs$delta)
       

            #' Get a colour scale gradient blue to red

            Trans_array <- NULL
            Trans_Pr_array <- NULL
            for(k in 0:(N_reinfect-1)){
                P_Ab = exp(-p*(V+q)) #' probability of being targetted by immune system. As binding avidity increases, this probability decreases
                P_s0 = 1-P_Ab #' probability of escape
                
                immK = r*k-delta #' Strength of host immune respone. As k increases, virus must escape more antibodies. As delta increases, this effect diminishes as infecting virus is further from host immunity.
                if(immK < 0) immK = 0
                
                P_Trans = (P_s0)^(immK) #' probability of escape from immunity
                
                if(k >= 1) P_Trans_Pr = immK*p*((1-P_Ab)^(immK-1))*(P_Ab) #' derivative of this. ie. rate of change of relationship between binding avidity and probability of immune escape
                else P_Trans_Pr = rep(0, length(V))
                Trans_array[[k+1]] <- P_Trans  
                Trans_Pr_array[[k+1]] <- P_Trans_Pr
            }

            
            P_Rep <- exp(-a*(V^b)) #' probability of successful infection within a host. as binding avidity increases in a naive host, chance of successfully replicating decreases
            P_Rep_Pr = -a*b*(V^(b-1))*exp(-a*(V^b)) #' rate of change of this relationship

            Rho_Trans_array <- NULL
            B_Pr_array <- NULL
            B_Trans_array <- NULL
            for(i in 1:length(Trans_Pr_array)){
                R0_Trans = Trans_array[[i]]*P_Rep*nv
                Rho_Trans = 1 - R0_Trans^-1
                Rho_Trans[Rho_Trans < 0] = 0
                Rho_Trans_array[[i]] <- Rho_Trans
                B_Pr = Trans_Pr_array[[i]]*P_Rep + Trans_array[[i]]*P_Rep_Pr
                B_Pr_array[[i]] <- B_Pr
                
                B_Trans = Trans_array[[i]]*P_Rep
                B_Trans_array[[i]] = B_Trans
            }

            rho0 = max(Rho_Trans_array[[1]])
            rho1 = max(Rho_Trans_array[[2]])
            rho2 = max(Rho_Trans_array[[3]])

            all_data <- NULL
            for(i in 1:length(Trans_array)){
                data <- data.frame(x=V,y=B_Trans_array[[i]],z=as.character(i))
                all_data <- rbind(all_data,data)
            }

            colours <- NULL
            blue <- rev(seq(0,1,by=1/N_reinfect))
            red <- seq(0,1,by=1/N_reinfect)
            for(i in 1:N_reinfect){
                colours[i] <- rgb((i-1)/(max_reinfect-1),0,(N_reinfect-i)/(max_reinfect-1))
                #'colours[i] <- rgb(red[i],0,blue[i])
            }

            #' Replication
            A <- ggplot() + geom_line(data=all_data,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
                theme(
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    axis.text.x=element_text(colour="gray20",size=12),
                    axis.text.y = element_text(colour="gray20",size=12),
                    text = element_text(colour="gray20",size=14),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x =element_line(colour="gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    legend.position = "none",
                    axis.title=element_text(size=12)
                ) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
                ylab("Probability of Successful Replication Within a Host") + 
                xlab("Binding Avidity")


            all_data2 <- NULL
            for(i in 1:length(Trans_array)){
                data <- data.frame(x=V,y=Rho_Trans_array[[i]],z=as.character(i))
                all_data2 <- rbind(all_data2,data)
            }

            #' Infection
            B <- ggplot() + geom_line(data=all_data2,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
                theme(
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    axis.text.x=element_text(colour="gray20",size=12),
                    axis.text.y = element_text(colour="gray20",size=12),
                    text = element_text(colour="gray20",size=14),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x =element_line(colour="gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    legend.position = "none",
                    axis.title=element_text(size=12)
                ) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
                ylab("Probability of Infection Between Hosts") + 
                xlab("Binding Avidity") +
                geom_hline(yintercept = rho0,linetype="longdash",colour="dodgerblue4")+
                geom_hline(yintercept = rho1,linetype="longdash",colour="dodgerblue4")+
                geom_hline(yintercept = rho2,linetype="longdash",colour="dodgerblue4")


            all_data3 <- NULL
            for(i in 1:length(Trans_array)){
                data <- data.frame(x=V,y=Trans_Pr_array[[i]],z=as.character(i))
                all_data3 <- rbind(all_data3,data)
            }

            #' Infection
            C <- ggplot() + geom_line(data=all_data3,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
                theme(
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    axis.text.x=element_text(colour="gray20",size=12),
                    axis.text.y = element_text(colour="gray20",size=12),
                    text = element_text(colour="gray20",size=14),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x =element_line(colour="gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    legend.position = "none",
                    axis.title=element_text(size=12)
                ) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0)) +
                ylab(expression(d*f/d*V)) + 
                xlab("Binding Avidity")


            all_data4 <- NULL
            for(i in 1:length(Trans_array)){
                data <- data.frame(x=V,y=B_Pr_array[[i]],z=as.character(i))
                all_data4 <- rbind(all_data4,data)
            }

            #' Infection
            D <- ggplot() + geom_line(data=all_data4,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
                theme(
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    axis.text.x=element_text(colour="gray20",size=12),
                    axis.text.y = element_text(colour="gray20",size=12),
                    text = element_text(colour="gray20",size=12),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x =element_line(colour="gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    legend.position = "none",
                    axis.title=element_text(size=12)
                ) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0)) +
                ylab(expression(d*beta/d*V)) + 
                xlab("Binding Avidity") +
                geom_hline(yintercept=0,linetype='longdash',colour="gray20")


            all_data5 <- NULL
            for(i in 1:length(Trans_array)){
                data <- data.frame(x=V,y=Trans_array[[i]],z=as.character(i))
                all_data5 <- rbind(all_data5,data)
            }

            E <- ggplot() + geom_line(data=all_data5,aes(x=x,y=y,colour=z)) + scale_color_manual(values=colours) +
                theme(
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    axis.text.x=element_text(colour="gray20",size=12),
                    axis.text.y = element_text(colour="gray20",size=12),
                    text = element_text(colour="gray20",size=14),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x =element_line(colour="gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    legend.position = "none",
                    axis.title=element_text(size=12)
                ) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0)) +
                ylab("Probability of Evading Immune System, f(k,V)") + 
                xlab("Binding Avidity")

            P_rep_dat <- data.frame(x=V,y=P_Rep)
            F <- ggplot() + geom_line(data=P_rep_dat,aes(x=x,y=y,colour="red")) +
                theme(
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    axis.text.x=element_text(colour="gray20",size=12),
                    axis.text.y = element_text(colour="gray20",size=12),
                    text = element_text(colour="gray20",size=14),
                    axis.line=element_line(colour="gray20"),
                    axis.line.x =element_line(colour="gray20"),
                    axis.line.y=element_line(colour="gray20"),
                    legend.position = "none",
                    axis.title=element_text(size=12)

                ) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    ylab("Probability of Successful Replication, g(V)") + 
    xlab("Binding Avidity")

            
            g <- grid.arrange(A,B,E,F,C,D,ncol=2)
        },height=1000,width=1000)
    })

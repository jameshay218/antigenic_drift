p = 1
r = 1
b = 3
a = 0.7
c = 1
nv = 2
q = 1
V = seq(0, 2, by = 0.005)
lifespan= 70
infectious_period = 5
N_reinfect = 20
max_reinfect = 20
D = 0.5
delta = 0.2
mu = 1/(lifespan*365.25)
gamma = 1/infectious_period

#' Get a colour scale gradient blue to red

Trans_array <- NULL
Trans_Pr_array <- NULL
for(k in 0:(N_reinfect-1)){
  P_Ab = exp(-p*(V+q))
  P_s0 = 1-P_Ab
  
  immK = r*k-delta
  if(immK < 0) immK = 0
  
  P_Trans = (P_s0)^(immK)
  
  if(k >= 1) P_Trans_Pr = immK*p*((1-P_Ab)^(immK-1))*(P_Ab)
  else P_Trans_Pr = rep(0, length(V))
  Trans_array[[k+1]] <- P_Trans  
  Trans_Pr_array[[k+1]] <- P_Trans_Pr
}

P_Rep <- exp(-a*(V^b))
P_Rep_Pr = -a*b*(V^(b-1))*exp(-a*(V^b))

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
                               legend.position = "none"
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
    legend.position = "none"
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
    legend.position = "none"
  ) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("f'(k,V)") + 
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
    text = element_text(colour="gray20",size=14),
    axis.line=element_line(colour="gray20"),
    axis.line.x =element_line(colour="gray20"),
    axis.line.y=element_line(colour="gray20"),
    legend.position = "none"
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
    legend.position = "none"
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
    legend.position = "none"
  ) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Probability of Successful Replication, g(V)") + 
  xlab("Binding Avidity")

  
g <- grid.arrange(A,B,E,F,C,D,ncol=2)
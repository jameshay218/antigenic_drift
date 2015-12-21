library(deSolve)

testing <- function(V, params){
  tmp <- NULL
  max <- 200
  for(i in 1:max){
      params[5] <- i
      tmp[i] <- y(V,params)
  }
  plot(tmp)
}


difeq <- function(V, t, params){
  p <- params[1]
  q <- params[2]
  r <- params[3]
  k <- params[4]
  delta <- params[5]
  a <- params[6]
  b <- params[7]
  c <- params[8]
    
  immK <- r*k - delta
  if(immK <0) immK <- 0
  
  f_x <- (1-exp(-p*(V+q)))^(immK)
  g_x <- exp(-a*(V^b))
  f_dx <- p*(r*k-delta)*((1-exp(-p*(V+q)))^(immK-1))*(exp(-p*(V+q)))
  g_dx <- -a*b*V^(b-1)*exp(-a*(V^b))
    
  dV <- f_x*g_dx + g_x*f_dx
  print(dV)
  dt <- dV*t*c
  return(dt)
}

y <- function(V, params){
  p <- params[1]
  q <- params[2]
  r <- params[3]
  k <- params[4]
  delta <- params[5]
  a <- params[6]
  b <- params[7]
  
  immK <- r*k - delta
  if(immK <0) immK <- 0
  
  f_x <- (1-exp(-p*(V+q)))^(immK)
  g_x <- exp(-a*(V^b))
  return(f_x*g_x)
}

sir_ode <- function (t, y, pars, C) 
{
  S <- y[1]
  I <- y[2]
  R <- y[3]
  N <- S + I + R
  beta <- pars[1]
  gamma <- pars[2]
  mu <- pars[3]
  dS <- mu*N -beta*((S*I)/N) - mu*S
  dI <- beta*((S*I)/N) - I*gamma - mu*I
  dR <- I*gamma - mu*R
  return(list(c(dS, dI, dR, beta*(S*I)/N)))
}

calculate_EE <- function(gamma, mu, beta){
  R0 <- (beta*mu)/(mu*(mu+gamma))
  S <- (gamma+mu)/beta
  I <- (mu/beta)*(R0-1)
  R <- (gamma/beta)*(R0-1)
  return(list(S,I,R))
}




N <- 1000000
pars <- c(0.5,0.333,1/100)
y0 <- c(900000,100,N-900000-100,0)
y1 <- ode(y0, seq(0,1000,by=1),sir_ode,pars)
sim <- read.csv("output.csv")
plot(sim[,2],ylim=c(0,200000),type='l',col="blue",lwd=3)

#lines(sim[,2],col="red")
#lines(y1[,3],col="green",lwd=5)

greb <- read.csv("voutput.csv")
#greb <- greb[sample(nrow(greb),10000,FALSE),]
#plot(greb$distRoot~greb$birth)
#plot(greb$immK~greb$birth)
# 
# 
# grebbons <- NULL
# for(i in 1:10){
#   grebbons[[i]] <- 1 - 1/(4*y(seq(0,2,by=0.01),c(3,1,45,(i-1)*10,1000,0.7,3)))
# }
# 
# plot(grebbons[[1]],type='l',ylim=c(0,1),col=1,lwd=2)
# for(i in 2:length(grebbons)) lines(grebbons[[i]],lwd=2,col=i)
# 
grebbons <- read.csv("voutput2.csv")
grebbons <- grebbons[,2:ncol(grebbons)]
test <- unlist(grebbons)
test <- na.omit(test)
hist(test)

greb1 <- greb[sample(nrow(greb),100,FALSE),]
for(i in 1:nrow(greb1)){
  target <- 0
  next_id <- greb1[i,"vid"]
  dist <- 0
  print(paste("Stored distance to root: ",greb1[i,"distRoot"],sep=""))
while(greb[greb$vid == next_id,"vid"] != target){
  
  dist <- dist + greb[greb$vid==next_id,"distance_to_parent"]
  next_id <- greb[greb$vid==next_id,"parentid"]
}
print(dist)
}

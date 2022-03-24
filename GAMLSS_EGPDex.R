library(gamlss)

options(error = rlang::entrace)
setwd("Documents/GAMLSS_EGPD/")

# 1. Create EGPD Family object
source ("GenericEGPDg.R")
EGPD1Family <- MakeEGPD (function (z,nu) z^nu, Gname = "Model1")

EGPD1Family(mu.link="identity",mu.init=1,
            sigma.link="log",sigma.init=3,
            nu.link="log",nu.init=0.5)

# From the EGPD1Family function, we create the object EGPDModel1 (name concatenating EGPD1 + Gname):
EGPDModel1 <- EGPD1Family()

# From there, we can use the functions:
# dEGPDModel1(x,mu=?,sigma=?,nu=?)  PDF
# pEGPDModel1(q,mu=?,sigma=?,nu=?)  CDF
# qEGPDModel1(q,mu=?,sigma=?,nu=?)  Inverse CDF (quantile function)
# rEGPDModel1(n,mu=?,sigma=?,nu=?)  Random number generation function

# show.link(EGPD1Family) #not yet fully implemented at the moment

# Create your own link function, here a shifted log
own.linkfun <- function (mu) {
  .shift = 0.0001
  log (mu - .shift) }
own.linkinv <- function (eta) {
  shift = 0.0001 
  thresh <- - log (.Machine$double.eps)
  eta <- pmin (thresh, pmax(eta, -thresh))
  exp (eta) + shift
} 
own.mu.eta <- function (eta) {
  shift = 0.0001  
  pmax (exp (eta), .Machine$double.eps)
}
own.valideta <- function (eta) TRUE

con <- gamlss.control (n.cyc = 150,
                       mu.step = 0.01, sigma.step = 0.01, nu.step = 0.01,tau.step = 0.1,autostep=TRUE)
# mu(and others).step = 0.1 is enough however for very low shape parameter distributions, 
# i.e with thin tails, 0.01 can help convergence without ending up with very high vu 
# (multiple local optima for log-likelihood)
# I imagine the model is less interesting that others in these situations...
con.i=glim.control(glm.trace = TRUE)

# Draw from EGPD1:
N <- 1000
p <- runif(N)
pr.Sim <- numeric(N)
for (i in 1:N) {
  pr.Sim[i] <- qEGPDModel1(p[i],mu=-0.1,sigma=1,nu=2)
}
Fn <- ecdf(pr.Sim)
plot(Fn,xlim=c(0,300))
par(mfrow=c(1,1))
hist(pr.Sim,xlab="x",ylab="Frequency",main="",cex.lab=1.5)

## Simulation study to assess we fit/estimate the same parameters that we would get 
# from a perfect and controlled EGPD simulation
mu0 <- -0.5
sigma0 <- 1
nu0 <- 2

Ny <- 20
N <- Ny*365
Nstat <- 10
db <- data.frame(matrix(ncol = 4, nrow = N*Nstat))
colnames(db) <- c("longitude","latitude","x","cyc")
Lat=runif(10,min=0,max=10)
Long=runif(10,min=0,max=10)
time=1:365
for (j in 1:Nstat){
  p <- runif(N)
  pr.Sim <- numeric(N)
  pr.Sim <- qEGPDModel1(p,mu=mu0,sigma=sigma0,nu=nu0)
  db$latitude[((j-1)*N+1):(j*N)] <- Lat[j] # repmat(Lat[j],N,1)
  db$longitude[((j-1)*N+1):(j*N)] <- Long[j] #repmat(Long[j],N,1)
  db$cyc[((j-1)*N+1):(j*N)] <- rep(time,Ny) #repmat(time,1,Ny)
  db$x[((j-1)*N+1):(j*N)] <-pr.Sim
}

# Empirical CDF
Fn <- ecdf(db$x)
# Compute theoretical CDF and PDF
if(mu0<0){
  q <- 0.1*(1:((-1/mu0)*10)) # support defined for 0 <= z <= -1/mu for mu < 0
}else{
  q <- 0.1*(1:(max(db$x)*10)) }
simCDF <- numeric(length(q))
simPDF <- numeric(length(q))
simCDF <- pEGPDModel1(q,mu=mu0,nu=nu0,sigma=sigma0)
simPDF <- dEGPDModel1(q,mu=mu0,nu=nu0,sigma=sigma0,log=FALSE)

# Sample directly from random number generation function
indTS <- rEGPDModel1(10000,mu=mu0,sigma=sigma0,nu=nu0)
Fn <- ecdf(indTS)

# Compare empirical and theoretical distributions
par(mfrow=c(1,3))
hist(db$x,xlab="x",ylab="Frequency",main="",cex.lab=1.5)

plot(q,simCDF,type ="b",cex=1.5,col="red",xlim=c(0,max(db$x)),ylim=c(0,1),xlab="x",ylab="CDF",cex.lab=1.5)
lines(Fn,lty = 1)
#lines(Fn,lty = 1)
legend("bottomright", legend=c("EGPD theorical", "EGPD simulation"),
       col=c("red", "black"), lty = c(NA,1),pch=c(1,NA), cex=1.2)

plot(q,simPDF,type ="b",cex=1.5,col="blue",xlim=c(0,max(db$x)),ylim=c(0,1),xlab="x",ylab="PDF",cex.lab=1.5)
lines(density(db$x),lty = 1)
legend("topright", legend=c("EGPD theorical", "EGPD simulation"),
       col=c("blue", "black"), lty = c(NA,1),pch=c(1,NA), cex=1.2)

## QQ plot
p<-(1:length(db$x))/(length(db$x)+1)
plot(qEGPDModel1(p,mu=mu0,nu=nu0,sigma=sigma0),sort(db$x), 
     xlab="Theoretical Quantiles",ylab="Sample Quantiles",
     main="QQ-plot")
abline(a=0,b=1)

## Fit EGPD (constant) parameters:
con <- gamlss.control (n.cyc = 100,
                       mu.step = 0.1, sigma.step = 0.1, nu.step = 0.1,tau.step = 0.1,autostep=TRUE)
con.i=glim.control(glm.trace = TRUE)
mod.egpd.0 <- gamlss(x ~ 1, 
                        data = db,
                        family = EGPD1Family(mu.link = "identity"),
                        # mu.link must be "identity" if we want to allow mu <0
                        control = con,
                        mu.start=1,
                        sigma.start=1,
                        nu.start=1,
                        #i.control=con.i,
                        method=CG())
# Compare fitted parameters to original ones (0.1,1,2)
muFit <- fitted(mod.egpd.0,"mu")[1]
sigmaFit <- predict(mod.egpd.0,what="sigma", type="response")[[1]]
nuFit <- predict(mod.egpd.0,what="nu", type="response")[[1]]
print(muFit)
print(sigmaFit)
print(nuFit)

## Simulation study but for a non-stationary field with cyclic dependence in day of the year.
# Namely, we still assume the linear trend in latitude, that is: mu(lat) = a lat + b
# with mu(lat=0)=mu0 and mu(lat=10)=mu10 for a given cyc
# but we also add a dependence in day of the year cyc, that is:
# mu(lat,cyc)= mu(lat)*(1+0.4*sin(2*cyc*pi/200)), for a given lat
# mu(lat,cyc)= ( ((mu10-mu0)/10)*lat+mu0 )*(1+0.4*sin(2*cyc*pi/200))

mu0 <- 0.1
mu10 <- 0.5
sigma0 <- 1
nu0 <- 2

Ny <- 20
time=1:365
N <- Ny*length(time)
Nstat <- 10
db <- data.frame(matrix(ncol = 4, nrow = N*Nstat))
colnames(db) <- c("longitude","latitude","x","cyc")
Lat=runif(10,min=0,max=10)
Long=runif(10,min=0,max=10)
for (j in 1:Nstat){
  mu_lat <- ((mu10-mu0)/10)*Lat[j]+mu0
  mu_lat_year <- mu_lat * (1+ 0.4* sin(2*time*pi/200))
  mu_lat_allyears <- rep(mu_lat_year,Ny)
  p <- runif(N)
  pr.Sim <- numeric(N)
  for (i in 1:N) {
    pr.Sim[i] <- qEGPDModel1(p[i],mu=mu_lat_allyears[i],sigma=sigma0,nu=nu0)
  }
  db$latitude[((j-1)*N+1):(j*N)] <- Lat[j] 
  db$longitude[((j-1)*N+1):(j*N)] <- Long[j] 
  db$cyc[((j-1)*N+1):(j*N)] <- rep(time,Ny)
  db$x[((j-1)*N+1):(j*N)] <-pr.Sim
}


# Create db2 for later cross-validation
Nstat <- 5
db2 <- data.frame(matrix(ncol = 4, nrow = N*Nstat))
colnames(db2) <- c("longitude","latitude","x","cyc")
Lat=runif(10,min=0,max=10)
Long=runif(10,min=0,max=10)
for (j in 1:Nstat){
  mu_lat <- ((mu10-mu0)/10)*Lat[j]+mu0
  mu_lat_year <- mu_lat * (1+ 0.4* sin(2*time*pi/200))
  mu_lat_allyears <- rep(mu_lat_year,Ny)
  p <- runif(N)
  pr.Sim <- numeric(N)
  for (i in 1:N) {
    pr.Sim[i] <- qEGPDModel1(p[i],mu=mu_lat_allyears[i],sigma=sigma0,nu=nu0)
  }
  db2$latitude[((j-1)*N+1):(j*N)] <- Lat[j] 
  db2$longitude[((j-1)*N+1):(j*N)] <- Long[j] 
  db2$cyc[((j-1)*N+1):(j*N)] <- rep(time,Ny)
  db2$x[((j-1)*N+1):(j*N)] <-pr.Sim
}

# Plot station maps

plot(latitude~longitude,data=db,pch=20,cex=2,col="blue",xlim=c(0,10),ylim=c(0,10),
     cex.lab=1.5,cex.axis=1.2)
points(latitude~longitude,data=db2,pch=1,cex=2,col="red")
legend("bottomleft", legend=c("db (fitting)", "db2 (validation)"),
       col=c("blue", "red"),pch=c(20,1), cex=1.5)

# Fit different models with EGPD1 or GA marginal distributions
con.i=glim.control(glm.trace = FALSE)

mod.egpd.0 <- gamlss(x~1, 
                     data=db,
                     family = EGPD1Family,
                     control = con,
                     mu.start=1,
                     sigma.start=1,
                     nu.start=1,
                     i.control=con.i,
                     #start.from = mod.egpd.mu.lat,
                     method=CG())

mod.egpd.1 <- gamlss(x~latitude+cyc, 
                     data=db,
                     family = EGPD1Family,
                     control = con,
                     mu.start=1,
                     sigma.start=1,
                     nu.start=1,
                     i.control=con.i,
                     start.from = mod.egpd.0,
                     method=CG())

mod.egpd.2 <- gamlss(x~latitude+pbc(cyc), 
                     data=db,
                     family = EGPD1Family,
                     control = con,
                     mu.start=1,
                     sigma.start=1,
                     nu.start=1,
                     i.control=con.i,
                     start.from = mod.egpd.0,
                     method=CG())

mod.egpd.3 <- gamlss(x~cs(latitude,longitude)+pbc(cyc), 
                     nu.fo=~cs(latitude,longitude)+pbc(cyc), 
                     data=db,
                     family = EGPD1Family,
                     control = con,
                     mu.start=1,
                     sigma.start=1,
                     nu.start=1,
                     i.control=con.i,
                     start.from = mod.egpd.0,
                     method=CG())

mod.ga.0 <- gamlss(x~1, 
                   data=db,
                   family = GA)

mod.ga.1 <- gamlss(x~latitude+cyc, 
                   data=db,
                   family = GA)

mod.ga.2 <- gamlss(x~latitude+pbc(cyc), 
                   data=db,
                   family = GA)

mod.ga.3 <- gamlss(x~cs(latitude,longitude)+pbc(cyc), 
                   nu.fo=~cs(latitude,longitude)+pbc(cyc), 
                   data=db,
                   family = GA)

# Visualise results on a model
mod <- mod.egpd.2

fittedPlot(mod, x=db$latitude, line.type = TRUE,xlab="latitude")
fittedPlot(mod, x=db$cyc, line.type = TRUE,xlab="Cyclic time index")

# Fit parameters for each value of covariates in the training set
muFit <- predict(mod.egpd.1,what="mu", type="response")
sigmaFit <- predict(mod.egpd.1,what="sigma", type="response")
nuFit <- predict(mod.egpd.1,what="nu", type="response")
# also: 
# muFit <- fitted(mod,"mu")

# Compare empirical and fitted CDF
indTS <- rEGPDModel1(length(muFit),mu=muFit,sigma=sigmaFit,nu=nuFit)

par(mfrow=c(1,3))
nb=5000
xl=c(0,30)
hist(db$x,breaks=2000,xlim=xl,freq=FALSE,xlab="x",ylab="Density",main="",cex.lab=1.5,cex.axis=1.5,col='skyblue',border=F)
hist(indTS,breaks=5000,xlim=xl,freq=FALSE,add=T,col=scales::alpha('red',.4),border=F,
     cex.axis=1.5,cex.lab=1,5)
title(main="",cex.main=2)
legend('topright',c("Empirical",'Simulated after fit'),
       fill = c('skyblue',scales::alpha('red',.5)), bty = 'n',
       border = NA,cex=1.2)

Fn <- ecdf(db$x)
plot(Fn,main=" ",cex=1.5,col='skyblue',xlim=xl,
     ylim=c(0,1),xlab="x",ylab="CDF",cex.lab=1.5,lwd=5)
# If constant parameters, we can compute it theoreticaly:
# curve(pEGPDModel1(x,mu=mu0,nu=nu0,sigma=sigma0),
#      from=min(db$x),to=max(db$x),add=TRUE,col="red",type="p")
Fn2 <- ecdf(indTS)
lines(Fn2,lty=2,col="red",lwd=2)
legend("bottomright", 
       legend=c("empirical","simulated"), col=c('skyblue',"red"), lty = c(1,1))

## QQ plot
par(mfrow=c(1,1))
p<-(1:length(db$x))/(length(db$x)+1)
plot(qEGPDModel1(p,mu=muFit,nu=nuFit,sigma=sigmaFit),sort(db$x), 
     xlab="Theoretical Quantiles",ylab="Sample Quantiles",
     main="QQ-plot")
abline(a=0,b=1)

par(mfrow=c(1,1))
lim=c(0,max(indTS,db$x))
plot(sort(indTS),sort(db$x), 
     xlab="Ordered simulated x",ylab="Ordered sampled x",
     main="",cex.lab=1.5,cex.axis=1.2,xlim=lim,ylim=lim)
abline(a=0,b=1)

#............................ Repeat but for out-of-sample stations - validation

# For new values of covariates:
#muFit2 <- predict(mod.egpd.1,what="mu",newdata=db2,data=db) 
#sigmaFit2 <- predict(mod.egpd.1,what="sigma",newdata=db2,data=db)
#nuFit2 <- predict(mod.egpd.1,what="nu",newdata=db2,data=db)
# or 
all <- predictAll(mod.egpd.1,newdata=db2,data=db)
muFit2 <- all$mu
sigmaFit2 <- all$sigma
nuFit2 <- all$nu

# Compare empirical and fitted CDF in a given station
is = 4
cond <- db2$latitude==Lat[is] & db2$longitude==Long[is]
indTS <- rEGPDModel1(length(muFit2[cond]),mu=muFit2[cond],sigma=sigmaFit2[cond],nu=nuFit2[cond])

par(mfrow=c(1,3))
xl=c(0,50)
hist(db2$x[cond],breaks=100,xlim=xl,freq=FALSE,xlab="x",ylab="Density",main="",cex.lab=1.5,cex.axis=1.5,col='skyblue',border=F)
hist(indTS,breaks=100,xlim=xl,freq=FALSE,add=T,col=scales::alpha('red',.4),border=F,
     cex.axis=1.5,cex.lab=1,5)
title(main="",cex.main=2)
legend('topright',c("Empirical",'Simulated after fit'),
       fill = c('skyblue',scales::alpha('red',.5)), bty = 'n',
       border = NA,cex=1.5)

Fn <- ecdf(db2$x[cond])
plot(Fn,main=" ",cex=1.5,cex.lab=1.5,cex.axis=1.5,col='skyblue',xlim=xl,
     ylim=c(0,1),xlab="x",ylab="CDF",cex.lab=1.5,lwd=5)
# If constant parameters, we can compute it theoreticaly:
# curve(pEGPDModel1(x,mu=mu0,nu=nu0,sigma=sigma0),
#      from=min(db$x),to=max(db$x),add=TRUE,col="red",type="p")
Fn2 <- ecdf(indTS)
lines(Fn2,lty=2,col="red",lwd=2)
legend("bottomright", cex=2,
       legend=c("empirical","simulated"), col=c('skyblue',"red"), lty = c(1,1))

#par(mfrow=c(1,1))
lim=c(0,max(indTS,db2$x[cond]))
plot(sort(indTS),sort(db2$x[cond]), 
     xlab="Ordered simulated x",ylab="Ordered sampled x",
     main="",cex.lab=1.5,cex.axis=1.5,xlim=lim,ylim=lim)
abline(a=0,b=1)

#............................... Compare empirical and fitted CDF - all stations
indTS <- rEGPDModel1(length(muFit2),mu=muFit2,sigma=sigmaFit2,nu=nuFit2)

par(mfrow=c(1,3))
xl=c(0,50)
hist(db2$x,breaks=200,xlim=xl,freq=FALSE,xlab="x",ylab="Density",main="",cex.lab=1.5,cex.axis=1.5,col='skyblue',border=F)
hist(indTS,breaks=200,xlim=xl,freq=FALSE,add=T,col=scales::alpha('red',.4),border=F,
     cex.axis=1.5,cex.lab=1,5)
title(main="",cex.main=2)
legend('topright',c("Empirical",'Simulated after fit'),
       fill = c('skyblue',scales::alpha('red',.5)), bty = 'n',
       border = NA,cex=1.5)

Fn <- ecdf(db2$x)
plot(Fn,main=" ",cex=1.5,cex.lab=1.5,cex.axis=1.5,col='skyblue',xlim=xl,
     ylim=c(0,1),xlab="x",ylab="CDF",cex.lab=1.5,lwd=5)
# If constant parameters, we can compute it theoreticaly:
# curve(pEGPDModel1(x,mu=mu0,nu=nu0,sigma=sigma0),
#      from=min(db$x),to=max(db$x),add=TRUE,col="red",type="p")
Fn2 <- ecdf(indTS)
lines(Fn2,lty=2,col="red",lwd=2)
legend("bottomright", cex=2,
       legend=c("empirical","simulated"), col=c('skyblue',"red"), lty = c(1,1))

#par(mfrow=c(1,1))
lim=c(0,max(indTS,db2$x))
plot(sort(indTS),sort(db2$x), 
     xlab="Ordered simulated x",ylab="Ordered sampled x",
     main="",cex.lab=1.5,cex.axis=1.5,xlim=lim,ylim=lim)
abline(a=0,b=1)


#......... Model quality
# Plot normalized quantile residuals against fitted values and against index
# Plot density estimate vs normalized quantile residuals
# Plot normal QQ plots
plot(mod)
# Worm plot
par(mfrow=c(1,1))
wp(mod)

# x as a function of latitude
par(mfrow=c(1,2))
plot(x~latitude, col="lightblue", data=db,cex.lab=1.5,ylab="x")

# Plot fitted parameter wrt covariates
muFit <- predict(mod,what="mu", type="response")
sigmaFit <- predict(mod,what="sigma", type="response")
nuFit <- predict(mod,what="nu", type="response")
par(mfrow=c(1,2))
plot(muFit~db$latitude,col="hotpink",cex=0.8,data=db,cex.lab=1.5,xlab="Latitude",ylab="mu")
lines(c(0,10),c(mu0,mu10),col="black", lty=2, lwd=2)
legend("bottomright",legend=c("Fitted","Theoretical"),
       pch = c(1,NA), lty = c(NA,2),col=c("hotpink","black"),cex=1.2)

plot(muFit~db$cyc,col="hotpink", data=db,cex=0.8,cex.lab=1.5,xlab="Cyclic time index",ylab="mu")
lines(c(0,10),c(mu0,mu10),col="black", lty=2, lwd=2)

# Plot predicted parameter wrt space/month
par(mfrow=c(1,2))
plot(muFit2~db2$latitude,col="hotpink",cex=0.8,data=db,cex.lab=1.5,xlab="Latitude",ylab="mu")
lines(c(0,10),c(mu0,mu10),col="black", lty=2, lwd=2)
legend("bottomright",legend=c("Fitted","Theoretical"),
       pch = c(1,NA), lty = c(NA,2),col=c("hotpink","black"),cex=1.2)

plot(muFit2~db2$cyc,col="hotpink", data=db,cex=0.8,cex.lab=1.5,xlab="Cyclic time index",ylab="mu")
lines(c(0,10),c(mu0,mu10),col="black", lty=2, lwd=2)


# Plot parameter maps
Lat=seq(0,10,by=0.2)
Long=seq(0,10,by=0.2)
db3 <- data.frame(matrix(ncol = 4, nrow = length(Lat)*length(Long)*2))
colnames(db3) <- c("longitude","latitude","x","cyc")
ii <- 1
for (j in 1:length(Lat)){
  for (i in 1:length(Long)){
    for (tt in c(1,90)){
  db3$latitude[ii] <- Lat[j] 
  db3$longitude[ii] <- Long[i] 
  db3$cyc[ii] <- tt
  db3$x[ii] <-NA
  ii <- ii+1
    }
  }
}

all <- predictAll(mod.egpd.2,newdata=db3[,c("longitude","latitude","cyc")],data=db)
muFit3 <- all$mu
sigmaFit3 <- all$sigma
nuFit3 <- all$nu

db3$mu <- muFit3
db3$sigma <- sigmaFit3
db3$nu <- nuFit3

tm <- min(db3$mu,na.rm = TRUE)
tM <- max(db3$mu,na.rm = TRUE)
mList <- c(1,150)
par(mfrow=c(1,2))
for(ii in c(1,90)){
  print(ii)
  cond = db3$cyc == ii
  xyz <- db3[cond,c("longitude","latitude","mu")]
  bubblePlot (as.matrix(xyz[,c("longitude","latitude")]),as.matrix(xyz[,c("mu")]), 
              xlab ="lon",ylab ="lat",
              col=tim.colors, highlight=FALSE,
              zlim=c(tm,tM),size=1.5,legend.cex = 1.2)
  title(paste(ii))
}

# Plot predicted parameter wrt space/month
muFit <- predict(mod.egpd.2,what="mu", type="response")
sigmaFit <- predict(mod.egpd.2,what="sigma", type="response")
nuFit <- predict(mod.egpd.2,what="nu", type="response")
par(mfrow=c(1,2))
plot(muFit~db$latitude,col="hotpink",cex=0.8,data=db,cex.lab=1.5,xlab="Latitude",ylab="mu")
lines(c(0,10),c(mu0,mu10),col="black", lty=2, lwd=2)
legend("bottomright",legend=c("Fitted","Theoretical"),
       pch = c(1,NA), lty = c(NA,2),col=c("hotpink","black"),cex=1.2)

plot(muFit~db$cyc,col="hotpink", data=db,cex=0.8,cex.lab=1.5,xlab="Cyclic time index",ylab="mu")
lines(c(0,10),c(mu0,mu10),col="black", lty=2, lwd=2)


# Model assessment by means of information criteria
GAIC(mod.egpd.0,mod.egpd.1,mod.egpd.2,mod.egpd.3,mod.ga.0,mod.ga.1,mod.ga.2,mod.ga.3,k=0)
GAIC(mod.egpd.0,mod.egpd.1,mod.egpd.2,mod.egpd.3,mod.ga.0,mod.ga.1,mod.ga.2,mod.ga.3)
GAIC(mod.egpd.0,mod.egpd.1,mod.egpd.2,mod.egpd.3,mod.ga.0,mod.ga.1,mod.ga.2,mod.ga.3,k=log(length(db$x)))

# On validation dataset

gg.egpd.0 <- getTGD(mod.egpd.0,newdata=db2)
gg.ga.0 <- getTGD(mod.ga.0,newdata=db2) 
gg.egpd.1 <- getTGD(mod.egpd.1,newdata=db2)
gg.ga.1 <- getTGD(mod.ga.1,newdata=db2) 
gg.egpd.2 <- getTGD(mod.egpd.2,newdata=db2)
gg.ga.2 <- getTGD(mod.ga.2,newdata=db2) 
gg.egpd.3 <- getTGD(mod.egpd.3,newdata=db2)
gg.ga.3 <- getTGD(mod.ga.3,newdata=db2) 

TGD(gg.egpd.0,gg.ga.0,gg.egpd.1,gg.ga.1,gg.egpd.2,gg.ga.2,gg.egpd.3,gg.ga.3)

# Since, in this example, only one explanatory variable is used in the fit, centile
# estimates for the fitted distribution can be shown using the functions centiles()
# or centiles.fan()
# ERROR
centiles.fan(mod.egpd.1, xvar=db$latitude, cent=c(3,10,25,50,75,90),
             colors="terrain",ylab="x", xlab="latitude")

# Repeat experience to test robustness

remov <- sample.int(length(db$x),length(db$x)/4)
db0 <- db
db$x[remov] <- NA
db <- na.omit(db)
db <- db[keep,]






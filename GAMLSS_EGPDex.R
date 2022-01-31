library(gamlss)

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
  mu_lat_allyears <- repmat(mu_lat_year,1,Ny)
  p <- runif(N)
  pr.Sim <- numeric(N)
  for (i in 1:N) {
    pr.Sim[i] <- qEGPDModel1(p[i],mu=mu_lat_allyears[i],sigma=sigma0,nu=nu0)
  }
  db$latitude[((j-1)*N+1):(j*N)] <- repmat(Lat[j],N,1)
  db$longitude[((j-1)*N+1):(j*N)] <- repmat(Long[j],N,1)
  db$cyc[((j-1)*N+1):(j*N)] <- repmat(time,1,Ny)
  db$x[((j-1)*N+1):(j*N)] <-pr.Sim
}

# Fit parameters (here only for mu, others are considered constant)
mod.0 <- gamlss(x~latitude+pbc(cyc), 
                              data=db,
                              family = EGPD1Family,
                              control = con,
                              mu.start=1,
                              sigma.start=1,
                              nu.start=1,
                              i.control=con.i,
                              #start.from = mod.egpd.mu.lat,
                              method=CG())

## Same with a constant mu and variable sigma
mu0 <- 0.1
sigma0 <- 1
sigma10 <- 2
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
  sigma_lat <- ((sigma10-sigma0)/10)*Lat[j]+sigma0
  for (i in 1:N) {
    pr.Sim[i] <- qEGPDModel1(p[i],mu=mu0,sigma=sigma_lat,nu=nu0)
  }
  db$latitude[((j-1)*N+1):(j*N)] <- repmat(Lat[j],N,1)
  db$longitude[((j-1)*N+1):(j*N)] <- repmat(Long[j],N,1)
  db$cyc[((j-1)*N+1):(j*N)] <- repmat(time,1,Ny)
  db$x[((j-1)*N+1):(j*N)] <-pr.Sim
}

mod.0 <- gamlss(x~1, 
                sigma.fo=~ latitude,
                data=db,
                family = EGPD1Family,
                #start.from = mod.egpd.0,
                control = con,
                mu.start=1,
                sigma.start=1,
                nu.start=1,
                method=CG())

## Simulation study to assess we fit/estimate the same parameters that we would get 
# from a perfect and controlled EGPD simulation - constant parameters
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
  #    for (i in 1:N) {
  #        pr.Sim[i] <- qEGPDModel1(p[i],mu=mu0,sigma=sigma0,nu=nu0)
  #    }
  pr.Sim <- qEGPDModel1(p,mu=mu0,sigma=sigma0,nu=nu0)
  db$latitude[((j-1)*N+1):(j*N)] <- repmat(Lat[j],N,1)
  db$longitude[((j-1)*N+1):(j*N)] <- repmat(Long[j],N,1)
  db$cyc[((j-1)*N+1):(j*N)] <- repmat(time,1,Ny)
  db$x[((j-1)*N+1):(j*N)] <-pr.Sim
}

mod.egpd.0 <- gamlss(x ~ 1, 
                        data = db, #[1:300,], #[1:300,],
                        family = EGPD1Family(mu.link = "identity"),
                        control = con,
                        mu.start=1,
                        sigma.start=1,
                        nu.start=1,
                        #i.control=con.i,
                        method=CG()) #CG() mixed(5,75)

# Fit parameters for each value of covariates in the training set
muFit <- predict(mod.egpd.0,what="mu", type="response")
sigmaFit <- predict(mod.egpd.0,what="sigma", type="response")
nuFit <- predict(mod.egpd.0,what="nu", type="response")
# For new values of covariates:
muFit2 <- predict(mod.egpd.o,what="mu",newdata=db2,type="response") 
sigmaFit2 <- predict(mod.egpd.0,what="sigma",newdata=db2,type="response")

# Compare empirical and fitted CDF
Fn <- ecdf(db$x)
plot(Fn,main=" ",xlim=c(0,quantile(db$x,0.999)))
# If constant parameters, we can cimompute it theoreticaly:
# curve(pEGPDModel1(x,mu=mu0,nu=nu0,sigma=sigma0),
#      from=min(db$x),to=max(db$x),add=TRUE,col="red",type="p")

indTS <- rEGPDModel1(length(muFit),mu=muFit,sigma=sigmaFit,nu=nuFit)
Fn <- ecdf(indTS)
lines(Fn,lty = 2,col="blue")
#legend("bottomright", 
#       legend=c("true", "empirical","simulated"), col=c("black","red","blue"), lty = c(1,NA,1),pch=c(NA,1,NA))
legend("bottomright", 
       legend=c("empirical","simulated"), col=c("black","blue"), lty = c(1,1))

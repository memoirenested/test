rm(list=ls())
par(mfrow = c(2, 2))

Evidence <- function(Y,u0=0,var0=1,var=1){
  S = sum(Y)
  S2 = sum(Y**2)
  n = length(Y)
  vpost = 1/(1/var0+n/var) #looool
  Z = -0.5*(u0**2/var0++S2/var) +0.5*(u0/var0+S/var)**2*vpost
  Z = Z - n/2*log(2*pi) - log(var0) -n*0.5*log(var) + 0.5*log(vpost)
  return(Z)
}

#s = 1
#m = 0
#Y = rnorm(20,mean=m,sd=s)
#Z = exp(Evidence(Y,u0=m,var0=s**2,var=1))

likelihood <- function(mu,y,v){
  L <- exp(sum(log(dnorm(y, mu, sd=v))))
  return(L)
}


nested_sampling <- function(n,j,N,mprior,vprior,vdata, data,ev=Z){
  E = ev
  T = runif(N,-10,10)
  Z = 0
  X = 1
  L = sapply(T,likelihood ,y = data ,v=sqrt(vdata))
  print(length(L))
  for(i in 1:j){
    L0 = min(L)
    X1 = exp(-i/n)
    w = X-X1
    Z = Z + L0*w
    print(Z)
    X = X1
    theta = rnorm(1,mean=mprior,sd=sqrt(vprior))
    Lnew = likelihood(mu=theta,y = data, v=sqrt(vdata))
    while(Lnew<=L0){
      theta = rnorm(1,mean=mprior,sd=sqrt(vprior))
      Lnew = likelihood(mu=theta,y = data, v=sqrt(vdata))
    }
    L[L==L0] = Lnew
  }
  Z = Z + (sum(L)-Lnew)*X/n
  return(c(Z,E,Z/E,log(Z)-log(E),Z-E))
}

nested_sampling2 <- function(n,j,N,mprior,vprior,vdata, data,ev=Z){
  E = ev
  T = rnorm(N,mean=mprior,sd=sqrt(vprior))
  Z = 0
  X = 1
  L = sapply(T,likelihood ,y = data ,v=sqrt(vdata))
  t = runif(N)
  for(i in 1:j){
    L0 = min(L)
    t0 = max(t)
    X1 = t0*X
    t[t==t0] = runif(1)
    w = X-X1
    Z = Z + L0*w
    X = X1
    theta = rnorm(1,mean=mprior,sd=sqrt(vprior))
    Lnew = likelihood(mu=theta,y = data, v=sqrt(vdata))
    while(Lnew<=L0){
      theta = rnorm(1,mean=mprior,sd=sqrt(vprior))
      Lnew = likelihood(mu=theta,y = data, v=sqrt(vdata))
    }
    L[L==L0] = Lnew
  }
  Z = Z + (sum(L)-Lnew)*X/n
  return(c(Z,E,Z/E,log(Z)-log(E),Z-E))
}

#####################################################################


f <- function(x, u, sigma, Lh, data){       # densité sous laquelle on veut simuler à une constante près
  y <- exp(-(x-u)*(x-u)/(2*sigma))*( likelihood(data, x, sqrt(sigma)) >= Lh )
  return(y)
}

q <- function(x, x_j, sigma){       # densité instrumentale, sous laquelle il est facile de simuler, ici N(x_j, sigma)
  y <- exp(-(x - x_j)*(x - x_j)/(2*sigma))
  return(y)
}

MCMC <- function(n_iteration, depart, u, sigma_f, sigma_q ,Lh, data){
  list_mcmc <- c(depart) 
  for (i in 1:n_iteration){
    actual_x <- list_mcmc[length(list_mcmc)]
    x_star <- rnorm(1, actual_x, sqrt(sigma_q))
    if(f(actual_x,  u , sigma_f, Lh, data)*q(x_star, actual_x, sigma_q) == 0){
      proba_acceptation <- 1
    } 
    else {
      proba_acceptation <- min(1, f(x_star, u , sigma_f, Lh, data)*q(actual_x, x_star, sigma_q)/
                                 (f(actual_x, u , sigma_f, Lh, data)*q(x_star, actual_x, sigma_q)))
    }
    new_sample <- sample(c(x_star, actual_x), size = 1, prob = c(proba_acceptation, 1 - proba_acceptation))
    list_mcmc <- append(list_mcmc, new_sample)
  }
  return(list_mcmc[length(list_mcmc)])                   # retourne dernier élément de la chaine
}


nested_samplingMC <- function(nombre_MCMC,pas_MCMC,j,N,mprior,vprior,vdata, data,ev=Z){
  E = ev
  #List_verif <- c()
  T = runif(N,-10,10)
  Z = 0
  X = 1
  L = sapply(T,likelihood,y = data ,v=sqrt(vdata))
  for(i in 1:j){
    L0 = min(L)
    X1 = exp(-i/N)
    w = X-X1
    Z = Z + L0*w
    X = X1
    i = order(L)[1]
    t = T[i]
    a <- MCMC(nombre_MCMC, t, mprior, vprior, pas_MCMC, L0, data)
    Lnew = likelihood(mu=a, y=data, v=sqrt(vdata))
    L[L==L0] = Lnew
    T[i] <- a
    #List_verif <- append(List_verif, Lnew - L0)
    
    '
    theta = rnorm(1,mean=mprior,sd=sqrt(vprior))
    Lnew = likelihood(data,mu=theta,v=sqrt(vdata))
    while(Lnew<=L0){
      theta = rnorm(1,mean=mprior,sd=sqrt(vprior))
      Lnew = likelihood(data,mu=theta,v=sqrt(vdata))
    }
    L[L==L0] = Lnew
    '
  }
  Z = Z + (sum(L)-Lnew)*X/N
  return(c(Z,E,Z/E,log(Z)-log(E),Z-E)) #List_verif
}

nested_samplingMC2 <- function(nombre_MCMC,n,j,N,mprior,vprior,vdata, data,ev=Z){
  E = ev
  T = rnorm(N,mean=mprior,sd=sqrt(vprior))
  Z = 0
  X = 1
  L = sapply(data,likelihood,mu=T,v=sqrt(vdata))
  t = runif(j)
  for(i in 1:j){
    L0 = min(L)
    t0 = max(t)
    X1 = t0*X
    t[t==t0]=runif(1)
    w = X-X1
    Z = Z + L0*w
    X = X1
    i = order(L)[1]
    ti = T[i]
    a <- MCMC(nombre_MCMC, 0, mprior, vprior, 1/j, ti, data)
    Lnew = likelihood(data,mu=a,v=sqrt(vdata))
    L[L==L0] = Lnew
    '
    theta = rnorm(1,mean=mprior,sd=sqrt(vprior))
    Lnew = likelihood(data,mu=theta,v=sqrt(vdata))
    while(Lnew<=L0){
      theta = rnorm(1,mean=mprior,sd=sqrt(vprior))
      Lnew = likelihood(data,mu=theta,v=sqrt(vdata))
    }
    L[L==L0] = Lnew
    '
  }
  Z = Z + (sum(L)-Lnew)*X/N
  return(c(Z,E,Z/E,log(Z)-log(E),Z-E))
}


#####variables locales 

md=0 #meandata
vd=1 #vardata
mp = #meanprior
vp = 1 #varprior

n = 20
j = 100
N <- 5000
j_mcmc <- 1000


D = rnorm(n,mean=md,sd=sqrt(vd)) #data
Z = exp(Evidence(D,u0=mp,var0=vp,var=vd))


##### estimateurs 

nest_mc = nested_samplingMC(50, 1/5, j_mcmc, N, mp, vp, vd, data=D)
print(nest_mc)

nest_mc2 = nested_samplingMC2(1, n, j_mcmc, N, mp, vp, vd, data=D)
nest = nested_sampling(n, j, N, mp, vp, vd, data=D)
nest2 = nested_sampling2(n, j, N, mp, vp, vd, data=D)

######### test smc

library(BayesianTools)
likelihood1 <- function(param){ #N(0,1) parametre moyenne
  singlelikelihoods = dnorm(D, mean = param[1], sd = 1,log= TRUE)
  return(sum(singlelikelihoods))  
}

bayesianSetup = createBayesianSetup(likelihood = likelihood1, lower = c(-10,0.01), upper = c(10,10))
iter = 1000

settings <- list(initialParticles = iter, iterations= 10)

out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "SMC", settings = settings)
a=marginalLikelihood(out)
smc_results=exp(a[["ln.ML"]])

print(smc_results-nest[1])
##############################
##### histogramme MCMC pas exponentiel

list_nest_mcmc_1_log <- c()
list_nest_mcmc_1_error <- c()
for (i in 1:200){
  nest_mc <- nested_samplingMC(50, 1/5, j_mcmc, N, mp, vp, vd, data=D)
  print(nest_mc)
  list_nest_mcmc_1_log <- append(list_nest_mcmc_1_log, log(nest_mc[3]))
  list_nest_mcmc_1_error <- append(list_nest_mcmc_1_error, nest_mc[5])
}

hist(list_nest_mcmc_1_log ,breaks=100)
hist(list_nest_mcmc_1_error, breaks=100)


##### histogramme MCMC pas uniforme

list_nest_mcmc_2_log <- c()
list_nest_mcmc_2_error <- c()
for (i in 1:200){
  nest_mc2 <- nested_samplingMC2(50, 1/5, j_mcmc, N, mp, vp, vd, data=D)
  list_nest_mcmc_2_log <- append(list_nest_mcmc_2_log, log(nest_mc2[3]))
  list_nest_mcmc_2_error <- append(list_nest_mcmc_2_error, nest_mc2[5])
}

hist(list_nest_mcmc_2_log ,breaks=100)
hist(list_nest_mcmc_2_error, breaks=100)


##### histogramme classique pas exponentiel

list_nest_1_log <- c()
list_nest_1_error <- c()
for (i in 1:200){
  nest <- nested_sampling(n, j, N, mp, vp, vd, data=D)
  list_nest_1_log <- append(list_nest_1_log, log(nest[3]))
  list_nest_1_error <- append(list_nest_1_error, nest[5])
}

hist(list_nest_1_log ,breaks=100)
hist(list_nest_1_error, breaks=100)


##### histogramme classique pas uniforme

list_nest_2_log <- c()
list_nest_2_error <- c()
for (i in 1:200){
  nest2 <- nested_sampling2(n, j, N, mp, vp, vd, data=D)
  list_nest_2_log <- append(list_nest_2_log, log(nest2[3]))
  list_nest_2_error <- append(list_nest_2_error, nest2[5])
}

hist(list_nest_2_log ,breaks=100)
hist(list_nest_2_error, breaks=100)

##### test smc 

smc <- function(N, echantillon_depart){   # echantillon depart represente une simulation sous le prior
  t <- 0
  L1 <- 0
  Z <- 1
  for(k in 1:N){
    
  }
}

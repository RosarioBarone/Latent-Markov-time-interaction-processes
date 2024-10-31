#########################################################################################
#########################################################################################
##########   Parametric MCMC inference for Hidden Time-Interaction Processes  ###########
#########################################################################################
#########################################################################################

lambda=function(psi,gamma,eta,X,j,individual.history){
  mu<-psi[1]
  theta<-psi[2]
  alpha<-psi[3]
  beta<-psi[4]
  #N<-j
  tk<-individual.history[j];
  # if(j>1) ts<-individual.history[1:(j-1)] else ts<-0
  # out=log(alpha*sum(exp(-beta*(tk-ts)))+eta*(tk^(eta-1))*exp(X*gamma+mu)*exp(-theta*(sum(individual.history<tk)-tk^eta)))
  #### Se la prima componente e' nulla nel caso della prima osservazione, allora ok.
  out=log(alpha*sum(exp(-beta*(individual.history[j]-individual.history[individual.history<tk])))+eta*(tk^(eta-1))*exp(X*gamma+mu)*exp(-theta*(sum(individual.history<tk)-tk^eta)))
  return(out)
}


lambda.integral.se=function(psi,gamma,eta,TT,individual.history){
  alpha<-psi[3]
  beta<-psi[4]
  p1=-(alpha/beta)*sum(exp(-beta*(TT-individual.history))-1) 
  return(p1)
}


library(statmod)
dd=gauss.quad(39,"chebyshev1")
fint=function(x,eta,theta,lower,upper){
  cv<-(upper-lower)*(x/2)+(lower+upper)/2
  out<-eta*cv^(eta-1)*exp(theta*cv^(eta))*sqrt(1-x^2)*(upper-lower)/2
  return(out)
}

lambda.integral.sc.0=function(psi,eta,gamma,X,individual.history){
  mu<-psi[1]
  theta<-psi[2]
  p2<-(exp(mu+gamma*X))*(exp(theta*individual.history[1]^eta)-1)/theta
  return(p2)
}


lambda.integral.sc=function(psi,gamma,eta,X,j,individual.history){
  mu<-psi[1]
  theta<-psi[2]
  N<-sum(individual.history<individual.history[j+1])
  if(j==1){
    p2=(exp(mu+gamma*X)*(exp(theta*individual.history[j]^eta)-1))/theta
    p2<-p2+exp(mu+gamma*X-theta*N)*(exp(theta*individual.history[j+1]^eta)-exp(theta*individual.history[j]^eta))/theta
  }else{
    p2=(exp(mu+gamma*X-theta*N)*(exp(theta*individual.history[j+1]^eta)-exp(theta*individual.history[j]^eta)))/theta
  }
  #p2<-(exp(mu+gamma*X-theta*N))*rowSums(fint(dd$nodes,eta,theta,individual.history[j],individual.history[j+1])*dd$weights)
  return(p2)
}

lambda.integral.last=function(psi,gamma,eta,X,N,Tk,TT){
  mu<-psi[1]
  theta<-psi[2]
  p2=(exp(mu+gamma*X-theta*N)*(exp(theta*TT^eta)-exp(theta*Tk^eta)))/theta
  return(p2)
}

lambda.integral.zerocaptures=function(psi,gamma,eta,X,TT){
  mu<-psi[1]
  theta<-psi[2]
  out=((exp(mu+gamma*X))*(exp(theta*TT^eta)-1))/theta
  return(out)
}


max.na<-function(x) max(x,na.rm = TRUE)


dTIP<-function(data,eta,gamma,psi,log){
  individual.history<-unlist(data$history)
  X<-data$X
  TT<-data$TT
  ### SETTING PER IL CASO CON PROCESSO MARKOVIANO LATENTE: QUI ASSUMIAMO CHE SIAMO SEMPRE NEL PRIMO STATO
  if(length(individual.history)==0){
    obslik<- -lambda.integral.zerocaptures(psi,gamma,eta,X[ceiling(TT)],TT)
    out<-obslik #-capture.conditioning(psi,gamma,eta,X,TT)
  }else{
    l0<-NULL
    sc.part<-0
    ##Part1
    j=1
    while(j<=length(individual.history)){
      l0<-c(l0,lambda(psi,gamma,eta,X[ceiling(individual.history[j])],j,individual.history))
      sc.part<-c(sc.part,lambda.integral.sc(psi,gamma,eta,X[ceiling(individual.history[j])],j,individual.history))
      j<-j+1 
    }
    sc.part<-c(na.omit(sc.part),lambda.integral.last(psi,gamma,eta,X[ceiling(TT)],N=length(individual.history),individual.history[length(individual.history)],TT))
    l0<-sum(l0)
    int1<-lambda.integral.se(psi,gamma,eta,TT,individual.history)
    int2<-sum(sc.part)
    obslik<-l0-(int1+int2)
    out<-obslik 
  }
  if(isFALSE(log)) out<-exp(out)
  return(out)
}

# sumlog<-function (x, lower = -745, upper = 709){
#   if (missing(x)) 
#     stop("'x' missing")
#   s <- tryCatch(.Call("fast_sumlog", x, lower, upper, length(x)), 
#                 `std::range_error` = function(e) {
#                   conditionMessage(e)
#                 })
#   
#   while(!is.finite(s)){
#     x<-x[-which.min(x)]
#     if(length(x)>1){
#       s <- tryCatch(.Call("fast_sumlog", x, lower, upper, length(x)), 
#                     `std::range_error` = function(e) {
#                       conditionMessage(e)
#                     }) 
#     }else{
#       s=x
#     }
#   }
#   #if(!is.finite(s)) stop("Not Summable!")
#   return(s)
# }

sumlog<-function(x) log(sum(exp(x)))

#MCMC.strategies: c("Joint", "Sequential")

MCMCsamplingTIPs<-function(data,MCMCiter,MCMC.strategy="Sequential",sigma.prop=0.05,burn.in=0.5){
  library(mvtnorm)
 # library(extraDistr)
  library(truncnorm)
  #### Set the algorithm ####
  gamma<-rep(0,MCMCiter)

  psi<-array(0,dim = c(MCMCiter,4))
  eta<-rep(0,nrow=MCMCiter)
  liks <- matrix(NA,length(data),MCMCiter)
  #Starting points
  gamma[1]=0.4
  psi[1,1]=0 #mu
  psi[1,2]=0.1 #theta
  psi[1,3]=0.4 #alpha
  psi[1,4]=0.45 #beta
  eta[1]=1

  #Updating strategies for the marginal posterior distributions
  Updating<-MCMC.strategy
  if(Updating=="Joint") s<-diag(sigma.prop,4)
  if(Updating=="Sequential") s<-sigma.prop
  for(iter in 2:MCMCiter){    
    if(iter%%10==0) {print(iter)}
    liks[,iter-1] <- unlist(lapply(data, function(data) dTIP(data,eta[iter-1],gamma[iter-1],psi[iter-1,],log=TRUE)))
    psi[iter,]<-psi[iter-1,]
    ##Update Psi|U # Joint Updating
    if(Updating=="Joint"){
        last.psi<-psi[iter,]
          if(iter>1000) s<-((2.38^2)/dim(psi)[2])*cov(psi[(iter-200):(iter-1),]) 
          psi.prop<-rmvnorm(1,c(last.psi[1],log(last.psi[2:4])),s)
          if(psi.prop[3]<psi.prop[4]){ 
            psi.prop[2:4]<-exp(psi.prop[2:4])
            lnum<-lapply(data, function(data) dTIP(data,eta[iter-1],gamma[iter-1],psi.prop,log=TRUE))
            lden<-lapply(data, function(data) dTIP(data,eta[iter-1],gamma[iter-1],last.psi,log=TRUE))
            lprior.num<-dnorm(psi.prop[1],mean = 0,sd=1000,log = TRUE)+dgamma(psi.prop[2],10,100,log = TRUE) +dunif(psi.prop[3],0,psi.prop[4],log = TRUE) +dunif(psi.prop[4],0,1000,log = TRUE)
            lprior.den<-dnorm(last.psi[1],mean = 0,sd=1000,log=TRUE)+dgamma(last.psi[2],10,100,log = TRUE) +dunif(last.psi[3],0,last.psi[4],log = TRUE) + dunif(last.psi[4],0,1000,log = TRUE)
            ACCEPT<-exp(sum(unlist(lnum))+lprior.num-sum(unlist(lden))-lprior.den)
            if(ACCEPT>runif(1)) psi[iter,]=psi.prop  
          }
      
    }
    
    
    if(Updating=="Sequential"){
      psi[iter,]<-psi[iter-1,]
          for (j in 1:2){
            last.psi<-psi[iter,]
            if(j==1) proposal<-rnorm(1,last.psi[j],s)
            if(j==2)  proposal<-exp(rnorm(1,log(last.psi[j]),s))
            psi.prop<-psi[iter,]
            psi.prop[j]<-proposal
            lnum<-lapply(data, function(data) dTIP(data,eta[iter-1],gamma[iter-1],psi.prop,log=TRUE))
            lden<-lapply(data, function(data) dTIP(data,eta[iter-1],gamma[iter-1],last.psi,log=TRUE))
            if(j==1){ 
              lprior.num<-dnorm(psi.prop[j],mean = 0,sd=10,log = TRUE)
              lprior.den<-dnorm(last.psi[j],mean = 0,sd=10,log = TRUE)
            }
            if(j==2){
              lprior.num<-dgamma(psi.prop[j],1,1,log = TRUE)
              lprior.den<-dgamma(last.psi[j],1,1,log = TRUE)
            }
            ACCEPT<-exp(sum(unlist(lnum))+lprior.num-sum(unlist(lden))-lprior.den)
            if(ACCEPT>runif(1)) psi[iter,j]=proposal 
          }
          #### Update Hawkes part ####
          j=3:4
          last.psi<-psi[iter,]
          proposal<-exp(rnorm(1,log(last.psi[j[1]]),s))
          proposal<-c(proposal,rtruncnorm(1,a=proposal,b=Inf,mean = last.psi[j[2]], sd=s))
          psi.prop<-psi[iter,]
          psi.prop[j]<-proposal
          lprior.num<-dunif(psi.prop[j[1]],0,psi.prop[j[2]],log = TRUE) + dunif(psi.prop[j[2]],0,15,log = TRUE)
          lprior.den<-dunif(last.psi[j[1]],0,last.psi[j[2]],log = TRUE) + dunif(last.psi[j[2]],0,15,log = TRUE)
          lnum<-lapply(data, function(data) dTIP(data,eta[iter-1],gamma[iter-1],psi.prop,log=TRUE))
          lden<-lapply(data, function(data) dTIP(data,eta[iter-1],gamma[iter-1],last.psi,log=TRUE))
          ACCEPT<-exp(sum(unlist(lnum))+lprior.num-sum(unlist(lden))-lprior.den)
          if(ACCEPT>runif(1)) psi[iter,j]=proposal 
      
    }
    #Update gamma
    last.gamma<-gamma[iter-1]
    gamma.prop<-rnorm(1,last.gamma,sigma.prop)
    lnum<-lapply(data, function(data) dTIP(data,eta[iter-1],gamma.prop,psi[iter,],log=TRUE))
    lden<-lapply(data, function(data) dTIP(data,eta[iter-1],last.gamma,psi[iter,],log=TRUE))
    ACCEPT<-exp(sum(unlist(lnum))+dnorm(gamma.prop,0,10,log = TRUE)-sum(unlist(lden))-dnorm(last.gamma,0,10,log=TRUE)); 
    if(ACCEPT>runif(1)) gamma[iter]=gamma.prop else gamma[iter]=last.gamma
    
    # Update the parametric baseline eta 
      last.eta<-eta[iter-1]
      eta.prop<-exp(rnorm(1,log(last.eta),sigma.prop))
      lnum<-lapply(data, function(data) dTIP(data,eta.prop,gamma[iter],psi[iter,],log=TRUE))
      lden<-lapply(data, function(data) dTIP(data,last.eta,gamma[iter],psi[iter,],log=TRUE))
      ACCEPT<-exp(sum(unlist(lnum))+dgamma(eta.prop,1,1,log = TRUE)-sum(unlist(lden))-dgamma(last.eta,1,1,log=TRUE)); 
      if(ACCEPT>runif(1)) eta[iter]=eta.prop else eta[iter]=last.eta
    }

  liks[,iter] <- unlist(lapply(data, function(data) dTIP(data,eta[iter],gamma[iter],psi[iter,],log=TRUE)))

  # Posterior summary statistics #
  burnin<-round(MCMCiter*burn.in)
  post.mean<-c(
    apply(psi[-(1:burnin),], 2, mean),mean(gamma[-(1:burnin)]),mean(eta[-(1:burnin)])
  )
  post.sd<-c(
    apply(psi[-(1:burnin),], 2, sd),sd(gamma),sd(eta)
  )
  quant <- function(x) c(quantile(x,0.025),quantile(x,0.975))
  cred.int<-cbind(apply(psi[-(1:burnin),], 2, quant),quant(gamma[-(1:burnin)]), quant(eta[-(1:burnin)]))
  tab<-as.table(rbind(post.mean,post.sd, cred.int))
  colnames(tab)<-c("mu","theta","alpha","beta","gamma","eta")
  #####
 out<-list() 
 liks<-liks[,-(1:burnin)]
 out$psi<-psi
 out$gamma<-gamma
 out$eta<-eta
 out$summary<-tab
 out$liks <- liks
 out$WAIC <- -2*(sum(apply(liks,1,sumlog)-log(1/iter))-sum(apply(liks,1,var)))
 return(out)
}


#par.mod<-MCMCsamplingTIPs(data,MCMCiter=500)
  

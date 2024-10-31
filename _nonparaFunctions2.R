#########################################################################################
#########################################################################################
######  semiparametric MCMC inference for semi-parametric Time-Interaction Processes   ######
#########################################################################################
#########################################################################################

library(snipEM) 


lambda=function(psi,gamma,lk,X,j,individual.history){
  mu<-psi[1]
  theta<-psi[2]
  alpha<-psi[3]
  beta<-psi[4]
  #N<-j
  tk<-individual.history[j];
  i_tj<-which(abs(lk[,1]-individual.history[j])==min(abs(lk[,1]-individual.history[j])))
  h<-(lk[i_tj,3]-lk[i_tj-1,3])/(lk[i_tj,1]-lk[i_tj-1,1])
  out=log(alpha*sum(exp(-beta*(individual.history[j]-individual.history[individual.history<tk])))+h*exp(X*gamma+mu)*exp(-theta*(sum(individual.history<tk)-lk[i_tj,3]) ))
  return(out)
}

lambda.integral.se=function(psi,gamma,TT,individual.history){
  alpha<-psi[3]
  beta<-psi[4]
  p1=-(alpha/beta)*sum(exp(-beta*(TT-individual.history))-1) 
  return(p1)
}

lambda.integral.sc=function(psi,gamma,lk,X,j,individual.history){
  mu<-psi[1]
  theta<-psi[2]
  N<-sum(individual.history<individual.history[j+1])
  #p2<-(exp(mu+gamma*X-theta*N))*lk[j+1]*(individual.history[j+1]-individual.history[j])*exp(theta*sum(lk[1:(j+1)]))
  i_tj1<-which(abs(lk[,1]-individual.history[j+1])==min(abs(lk[,1]-individual.history[j+1])))
  i_tj0<-which(abs(lk[,1]-individual.history[j])==min(abs(lk[,1]-individual.history[j])))
  p2<-(exp(mu+gamma*X-theta*N))*(exp(theta*lk[i_tj1,3])-exp(theta*lk[i_tj0,3]))/theta
  return(p2)
}

lambda.integral.sc.0=function(psi,gamma,lk,X,individual.history){
  mu<-psi[1]
  theta<-psi[2]
  i_tj<-which(abs(lk[,1]-individual.history[1])==min(abs(lk[,1]-individual.history[1])))
  p2<-(exp(mu+gamma*X))*(exp(theta*lk[i_tj,3])-1)/theta
  return(p2)
}


lambda.integral.sc.last=function(psi,gamma,lk,X,N,Tk,TT,individual.history){
  mu<-psi[1]
  theta<-psi[2]
  j<-length(individual.history)
  i_tj<-which(abs(lk[,1]-individual.history[j])==min(abs(lk[,1]-individual.history[j])))
  i_TT<-which(abs(lk[,1]-TT)==min(abs(lk[,1]-TT)))
  p2<-exp(mu+gamma*X-theta*N)*(exp(theta*lk[i_TT,3])-exp(theta*lk[i_tj,3]))/theta
  return(p2)
}

lambda.integral.zerocaptures=function(psi,gamma,lk,X,TT){
  mu<-psi[1]
  theta<-psi[2]
  out<-(exp(mu+gamma*X))*(exp(theta*lk[dim(lk)[1],3])-1)
  return(out)
}


max.na<-function(x) max(x,na.rm = TRUE)


dTIP<-function(data,gamma,psi,log){
  individual.history<-unlist(data$history)
  X<-data$X
  TT<-data$TT
  #lk<-data$baseline
  lk<-matrix(unlist(data$baseline),ncol=3,byrow=FALSE)
  ### SETTING PER IL CASO CON PROCESSO MARKOVIANO LATENTE: QUI ASSUMIAMO CHE SIAMO SEMPRE NEL PRIMO STATO
  if(length(individual.history)==0){
    obslik<- -lambda.integral.zerocaptures(psi,gamma,lk,X[ceiling(TT)],TT)
    out<-obslik 
  }else{
    l0<-0
    sc.part<-lambda.integral.sc.0(psi,gamma,lk,X[ceiling(individual.history[1])],individual.history)
    ##Part1
    j=1
    while(j<=length(individual.history)){
      l0<-c(l0,lambda(psi,gamma,lk,X[ceiling(individual.history[j])],j,individual.history))
      sc.part<-c(sc.part,lambda.integral.sc(psi,gamma,lk,X[ceiling(individual.history[j])],j,individual.history))
      j<-j+1 
    }
    sc.part<-c(sc.part,lambda.integral.sc.last(psi,gamma,lk,X[ceiling(TT)],N=length(individual.history),Tk=individual.history[length(individual.history)],TT,individual.history))
    l0<-sum(l0)
    int1<-lambda.integral.se(psi,gamma,TT,individual.history)
    int2<-sum(sc.part)
    obslik<-l0-(int1+int2)
    out<-obslik #-log(1-exp(-capture.conditioning(psi,gamma,eta,X,TT)))  
  }
  if(isFALSE(log)) out<-exp(out)
  return(out)
}

rgamma.process.augm<-function(data,from,to,a,eta,time.partition){
  if(is.null(data)) TT=to else TT<-unlist(data$TT)
  x<-seq(from=0,to=TT, length.out=time.partition)
  x1<-x[-1]
  x0<-x[-length(x)]
  proc<-rgamma(length(diff(x)),shape=a*(x1^(eta)-x0^(eta)),rate=a)
  out<-cbind(x,c(0,proc),cumsum(c(0,proc)))
  ### Computational issue ####
  out[out[,2]==0,2]=NA
  out<-na.omit(out)
  out<-rbind(rep(0,3),out)
  colnames(out)<-c("t","h","H")
  return(out)
}

rlgauss.process.augm<-function(data,sigma.proposal,time.partition){
  TT<-unlist(data$TT)
  captures<-c(0,unlist(data$history),TT)
  x<-seq(from=0,to=TT, length.out=time.partition)
  x<-sort(unique(c(captures,x)))
  baseline<-data$baseline
  proc<-exp(rnorm(length(baseline[-1,2]),mean=log(baseline[-1,2]),sd=sigma.proposal))
  out<-cbind(x,c(0,proc),cumsum(c(0,proc)))
  ### Computational issue ####
  out[out[,2]==0,2]=NA
  out<-na.omit(out)
  out<-rbind(rep(0,3),out)
  colnames(out)<-c("t","h","H")
  return(out)
}


dgamma.process<-function(y,a,eta,log){
  x1<-y[-1,1]
  x0<-y[-dim(y)[1],1]
  out<-sum(dgamma(y[-1,2],shape=a*(x1^(eta)-x0^(eta)),rate=a, log=TRUE))
  if(isTRUE(log)){
    out<-out
  }else{
    out<-exp(out)
  }
  return(out)
}



semiparametric_MCMCsamplingTIPs<-function(data,MCMCiter,MCMC.strategy,sigma.prop=0.05,random.gamma=TRUE,time.partition=100,a.fix,eta=1,burn.in=0.5){
  library(mvtnorm)
  library(extraDistr)
  library(truncnorm)
  TT=data[[1]]$TT
  n=length(data)
  #### Set the algorithm ####
  liks <- matrix(NA,length(data),MCMCiter)
  gamma<-rep(0,MCMCiter)
  gamma[1]=0
  psi<-array(0,dim = c(MCMCiter,4))
  # Introducing a random baseline on the data oject
  lambda0<-list(rgamma.process.augm(data=NULL,from=0,to=TT,a=a.fix,eta=1, time.partition=time.partition))
  names(lambda0)<-"baseline"
  data<- lapply(1:n,FUN=function(j) append(data[[j]],lambda0))
  #Starting points
  psi[1,1]=0
  psi[1,2]=0.5
  psi[1,3]=0.5
  psi[1,4]=0.5
  post.baseline<-array(NA,dim=c(MCMCiter,time.partition,3))
  post.baseline[1,,]=matrix(unlist(lambda0),ncol=3,byrow=FALSE)
  #Updating strategies for the marginal posterior distributions
  Updating<-MCMC.strategy
  if(Updating=="Joint") s<-diag(sigma.prop,4)
  if(Updating=="Sequential") s<-sigma.prop
  if(isTRUE(random.gamma)){
    a<-rep(a.fix,MCMCiter) #fixed a
    eta<-rep(0,MCMCiter) ; eta[1]<-1 
    sab.prop<-diag(sigma.prop,2)
  }else{
    a<-rep(a.fix,MCMCiter) #fixed a
    eta<-rep(eta,MCMCiter) #fixed eta 
  }
  
  for(iter in 2:MCMCiter){
    if(iter%%10==0) {print(iter)}
    
    if(iter%%100==0){
      if(iter>1000) burnin=iter-1000 else burnin=1
      par(mar=c(1,1,1,1))
      par(mfrow=c(5,1))
      plot(psi[burnin:(iter-1),1], type = "l", xlab = "", ylab = "", main="mu")
      plot(psi[burnin:(iter-1),2], type = "l", xlab = "", ylab = "", main="theta")
      plot(psi[burnin:(iter-1),3], type = "l", xlab = "", ylab = "", main="alpha")
      plot(psi[burnin:(iter-1),4], type = "l", xlab = "", ylab = "", main="beta")
      # par(mfrow=c(2,1))
      psi.post.mean<-c(t(apply(psi[burnin:(iter-1),],2,mean)), mean(gamma[burnin:(iter-1)]), mean(eta[burnin:(iter-1)]))
      names(psi.post.mean)<-c("mu","theta","alpha","beta", "gamma", "eta")
      print(psi.post.mean)
      waic<--2*(sum(apply(liks[,burnin:(iter-2)],1,sumlog)-log(1/(iter-burnin)))-sum(apply(liks[,burnin:(iter-2)],1,var)))
      print(waic)
    }
    
    
    psi[iter,]<-psi[iter-1,]
    liks[,iter-1]<-unlist(lapply(data, function(data) dTIP(data,gamma[iter-1],psi[iter-1,],log=TRUE)))
    ##Update Psi|U # Joint Updating
    if(Updating=="Joint"){
      last.psi<-psi[iter,]
      if(iter>200) s<-((2.38^2)/dim(psi)[2])*cov(psi[(iter-200):(iter-1),]) 
      if(s[2,2]<0.0001) s[2,2]=0.0001
      psi.prop<-rmvnorm(1,c(last.psi[1],log(last.psi[2:4])),s)
      psi.prop[1]<-0
      lden<-lapply(data, function(data) dTIP(data,gamma[iter-1],last.psi,log=TRUE))
      if(psi.prop[3]<psi.prop[4]){ 
        psi.prop[2:4]<-exp(psi.prop[2:4])
        lnum<-lapply(data, function(data) dTIP(data,gamma[iter-1],psi.prop,log=TRUE))
        lprior.num<-dnorm(psi.prop[1],mean = 0,sd=1000,log = TRUE)+dgamma(psi.prop[2],10,100,log = TRUE) +dunif(psi.prop[3],0,psi.prop[4],log = TRUE) +dunif(psi.prop[4],0,1000,log = TRUE)
        lprior.den<-dnorm(last.psi[1],mean = 0,sd=1000,log=TRUE)+dgamma(last.psi[2],10,100,log = TRUE) +dunif(last.psi[3],0,last.psi[4],log = TRUE) + dunif(last.psi[4],0,1000,log = TRUE)
        ACCEPT<-exp(sum(unlist(lnum))+lprior.num-sum(unlist(lden))-lprior.den)
        if(ACCEPT>runif(1)) psi[iter,]=psi.prop  
      }
      
    }
    
    
    if(Updating=="Sequential"){
      #Update theta
      psi[iter,]<-psi[iter-1,]
      last.psi<-psi[iter,]
      proposal<-exp(rnorm(1,log(last.psi[2]),s))
      psi.prop<-psi[iter,]
      psi.prop[2]<-proposal
      lnum<-lapply(data, function(data) dTIP(data,gamma[iter-1],psi.prop,log=TRUE))
      lden<-lapply(data, function(data) dTIP(data,gamma[iter-1],last.psi,log=TRUE))
      lprior.num<-dgamma(psi.prop[2],1,1,log = TRUE)
      lprior.den<-dgamma(last.psi[2],1,1,log = TRUE)
      ACCEPT<-exp(sum(unlist(lnum))+lprior.num-sum(unlist(lden))-lprior.den)
      if(ACCEPT>runif(1)) psi[iter,2]=proposal 
      #### Update Hawkes part ####
      j=3:4
      last.psi<-psi[iter,]
      proposal<-exp(rnorm(1,log(last.psi[j[1]]),s))
      proposal<-c(proposal,rtruncnorm(1,a=proposal,b=Inf,mean = last.psi[j[2]], sd=s))
      psi.prop<-psi[iter,]
      psi.prop[j]<-proposal
      lprior.num<-dunif(psi.prop[j[1]],0,psi.prop[j[2]],log = TRUE) + dunif(psi.prop[j[2]],0,100,log = TRUE)
      lprior.den<-dunif(last.psi[j[1]],0,last.psi[j[2]],log = TRUE) + dunif(last.psi[j[2]],0,100,log = TRUE)
      lnum<-lapply(data, function(data) dTIP(data,gamma[iter-1],psi.prop,log=TRUE))
      lden<-lapply(data, function(data) dTIP(data,gamma[iter-1],last.psi,log=TRUE))
      ACCEPT<-exp(sum(unlist(lnum))+lprior.num-sum(unlist(lden))-lprior.den)
      if(ACCEPT>runif(1)) psi[iter,j]=proposal 
      
    }
    #Update gamma
    last.gamma<-gamma[iter-1]
    gamma.prop<-rnorm(1,last.gamma,sigma.prop)
    lnum<-lapply(data, function(data) dTIP(data,gamma.prop,psi[iter,],log=TRUE))
    lden<-lapply(data, function(data) dTIP(data,last.gamma,psi[iter,],log=TRUE))
    ACCEPT<-exp(sum(unlist(lnum))+dnorm(gamma.prop,0,10,log = TRUE)-sum(unlist(lden))-dnorm(last.gamma,0,10,log=TRUE)); 
    if(ACCEPT>runif(1)) gamma[iter]=gamma.prop else gamma[iter]=last.gamma
    
    last.lambda0<-post.baseline[iter-1,,]
    prop.lambda0<-list(rgamma.process.augm(data=NULL,from=0,to=TT,a=a.fix,eta=eta[iter-1], time.partition=time.partition))
    data.proposal=data
    for (i in 1:n) data.proposal[[i]]$baseline=prop.lambda0 ### scrivere meglio
    lnum<-lapply(1:n,function(j) dTIP(data.proposal[[j]],gamma[iter],psi[iter,],log=TRUE))
    lden<-lapply(1:n,function(j) dTIP(data[[j]],gamma[iter],psi[iter,],log=TRUE))
    prior.num<- sum(dgamma(matrix(unlist(prop.lambda0),ncol = 3)[,2],shape=1,rate=1,log=TRUE))
    prior.den<- sum(dgamma(last.lambda0[,2],shape=1,rate=1,log=TRUE))
    ACCEPT<-exp(sum(unlist(lnum))+prior.num+ dgamma.process(y=last.lambda0,a=a[iter-1],eta[iter-1],log=TRUE) -sum(unlist(lden))-prior.den- dgamma.process(y=matrix(unlist(prop.lambda0),ncol = 3),a=a[iter-1],eta[iter-1],log=TRUE) ); 
    #ACCEPT<-exp(lnum+prior.num-lden-prior.den);
    if(ACCEPT>runif(1)){
      data=data.proposal
      post.baseline[iter,,]=matrix(unlist(prop.lambda0),ncol = 3)
      }else{
      data=data
      post.baseline[iter,,]=last.lambda0
      } 
    
    if(isTRUE(random.gamma)){
     S.bl<-0.05
     last<-eta[iter-1]
     prop<-exp(rnorm(1,mean = log(last),sd=S.bl))
     ACCEPT<- exp(dgamma.process(post.baseline[iter,,],a=a.fix,eta =prop, log = TRUE) + dgamma(prop,1,1, log = TRUE) 
         - dgamma.process(post.baseline[iter,,],a=a.fix,eta =last, log = TRUE) - dgamma(last,1,1, log = TRUE) )  
     if(ACCEPT>runif(1)){
       eta[iter]=prop
     }else{
       eta[iter]=last
     }
     
     }
   
  }
  
  
  liks[,iter]<-unlist(lapply(data, function(data) dTIP(data,gamma[iter],psi[iter,],log=TRUE)))
  burnin<-round(MCMCiter*burn.in)
  # Posterior summary statistics #
  out<-list() 
  out$psi<-psi
  out$gamma<-gamma
  out$baseline<-post.baseline
  #####
  if(isTRUE(random.gamma)){
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
    out$eta<-eta
  }else{
    post.mean<-c(
      apply(psi[-(1:burnin),], 2, mean),mean(gamma[-(1:burnin)])
    )
    post.sd<-c(
      apply(psi[-(1:burnin),], 2, sd),sd(gamma)
    )
    quant <- function(x) c(quantile(x,0.025),quantile(x,0.975))
    cred.int<-cbind(apply(psi[-(1:burnin),], 2, quant),quant(gamma[-(1:burnin)])) 
    tab<-as.table(rbind(post.mean,post.sd, cred.int))
    colnames(tab)<-c("mu","theta","alpha","beta","gamma")
  }
  out$summary<-tab
  out$WAIC <- -2*(sum(apply(liks[,-(1:burnin)],1,sumlog)-log(1/iter))-sum(apply(liks[,-(1:burnin)],1,var)))
  out$liks=liks
  return(out)
}

 # posterior.baseline<-post.baseline[,,3]
 # mean.baseline<-apply(posterior.baseline, 2, mean)
 # sd.baseline<-apply(posterior.baseline, 2, sd)
 # plot(mean.baseline,type = "l", lwd=2, col="red")
 # 
 # burnin=500
 # post.baseline<-post.baseline[,,-(1:burnin)]
 # 
 # dim(post.baseline)
 # 
 # sum(AR)/length(AR)
 # plot(post.baseline[1,,500])
 # plot(post.baseline[2,,500])

# 

SEQsim<-semiparametric_MCMCsamplingTIPs(data=data, MCMCiter=1500,MCMC.strategy = "Sequential",random.gamma=TRUE,a.fix=10)
#JOINTsim<-semiparametric_MCMCsamplingTIPs(data=data, MCMCiter=2000,MCMC.strategy = "Joint",random.gamma=TRUE,a.fix=1)
# # JOINTsim$summary
# SEQsim$summary
# 
# plot(SEQsim$psi[,1], type = "l")
# plot(SEQsim$psi[,2], type = "l")
# plot(SEQsim$psi[,3], type = "l")
# plot(SEQsim$eta, type = "l")
# plot(SEQsim$a, type = "l")

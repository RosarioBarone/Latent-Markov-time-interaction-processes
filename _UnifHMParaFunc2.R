#########################################################################################
#########################################################################################
######  MCMC inference for parametric Hidden Markov Time-Interaction Processes   ########
#########################################################################################
#########################################################################################

library(compiler)

lambda=function(psi,gamma,eta,X,j,individual.history){
  mu<-psi[1]
  theta<-psi[2]
  alpha<-psi[3]
  beta<-psi[4]
  #N<-j
  tk<-individual.history[j];
  out=log(alpha*sum(exp(-beta*(individual.history[j]-individual.history[individual.history<tk])))+eta*(tk^(eta-1))*exp(X*gamma+mu)*exp(-theta*(sum(individual.history<tk)-tk^eta)))
  return(out)
}


lambda.integral.se=function(psi,gamma,eta,TT,individual.history){
  alpha<-psi[3]
  beta<-psi[4]
  p1=-(alpha/beta)*sum(exp(-beta*(TT-individual.history))-1) 
  return(p1)
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


lambda.integral.zerocaptures.interval=function(psi,gamma,eta,X,t1,t0){
  mu<-psi[1]
  theta<-psi[2]
  out=((exp(mu+gamma*X))*(exp(theta*t1^eta)-exp(theta*t0^eta)))/theta
  return(out)
}

#### FUNZIONE PER CALCOLARE LA VEROSIMIGLIANZA DATE LE TRAIETTORIE INTERAMENTE OSSERVATE

ind.cond.lik<-function(data,st,eta,gamma,psi,log){
  Hidden.Process<-data$U
  states<-Hidden.Process[,1]
  jumps<-Hidden.Process[,2]
  individual.history<-unlist(data$history)
  loglik<-0
  if(length(individual.history)==0&states[1]==st){
    TT<-data$TT
    X<-data$X
    loglik<- -lambda.integral.zerocaptures(psi,gamma,eta,X[ceiling(TT)],TT)
  }else{
    t<-c()
    for (j in individual.history) {
      state<-Hidden.Process[min(which(jumps>j)),1]
      if(state==st) t<-c(t,j)
    }
    if(length(t)==1) loglik<-dTIP.event(data,t=t,eta,gamma,psi[,st])
    if(length(t)>1)  loglik<-sum(sapply(1:length(t),function(i) dTIP.event(data,t=t[i],eta,gamma,psi[,st])))
  }
  return(loglik)
}

#### FUNZIONE PER CALCOLARE LA DENSITA' IN UN PUNTO OSSERVATO (IN CORRISPONDENZA DI UN EVENTO)
dTIP.event<-function(data,t,eta,gamma,psi){
  individual.history<-unlist(data$history)
  TT<-data$TT
  N<-individual.history[individual.history<=t]
  last=FALSE
  if(t==max(individual.history)) last=TRUE 
  X<-data$X
  l0<-lambda(psi,gamma,eta,X[ceiling(t)],which(individual.history==t),individual.history)
  if(isFALSE(last))int2 <-lambda.integral.sc(psi,gamma,eta,X[ceiling(t)],which(individual.history==t),individual.history)
  if(isTRUE(last))int2<-lambda.integral.last(psi,gamma,eta,X[ceiling(TT)],N=length(individual.history),individual.history[length(individual.history)],TT)
  int1<-lambda.integral.se(psi,gamma,eta,TT,t)
  out<-l0-(int1+int2)
  return(out)
}


####################################################################################
####################################################################################
####################################################################################


### FUNZIONI PER UNIFORMIZZAZIONE


UniformSampl <- function(BgnSt,EndSt,RateM,Tm){
  
  ptm1 <- proc.time()[1]
  ## Diagonalization of rate matrix
  Eigen <- eigen(RateM)
  Lam <- Eigen$values
  U <- Eigen$vectors
  InvU <- solve(U)
  PrbBgnEnd <- (U%*%diag(exp(Tm*Lam))%*%InvU)[BgnSt,EndSt]
  if (is.complex(PrbBgnEnd)){
    PrbBgnEnd = Re(PrbBgnEnd)
  }
  ## Determine max diagonal entry and construct transition matrix
  nSt <- nrow(RateM)
  Mx <- max(-diag(RateM))
  TransM <- diag(rep(1,times=nSt))+RateM/Mx
  TransMn <- TransM
  ## TransMn is n'th power of TransM
  ##------------------------------------------------------------
  ## Simulate number of jumps
  rU <- runif(1)
  ifelse(BgnSt==EndSt,cum <- dpois(0,Mx*Tm)/PrbBgnEnd,cum <- 0)
  notExceed <- TRUE
  if (cum>rU) { notExceed <- FALSE }
  nJmp <- 0
  ptm2 <- proc.time()[1]
  while (notExceed){
    nJmp <- nJmp+1 
    prb <- dpois(nJmp,Mx*Tm)*TransMn[BgnSt,EndSt]/PrbBgnEnd
    cum <- cum+prb
    if (cum>rU) notExceed <- FALSE
    ## Update transition matrices
    ## TransArr holds n'th power of TransM
    if (nJmp==1) TransArr <- array(TransM,c(nSt,nSt,1))
    if (nJmp!=1) TransArr <- array(c(TransArr,TransMn),c(nSt,nSt,nJmp))
    TransMn <- TransMn%*%TransM
  }
  #cat("nJmp:",nJmp,"\n")
  ptm3 <- proc.time()[1]
  ##--------------------------------------------------------
  ## if (nJmp==0): done
  if (nJmp==0){
    Path <- list()
    Path$St <- c(BgnSt,EndSt)
    Path$Tm <- c(0,Tm)
    Path$ptm <- c(ptm2-ptm1,ptm3-ptm2)
    return(Path)
  }
  ## if (nJmp==1)
  if (nJmp==1){
    ## if virtual jump: done
    if (BgnSt==EndSt){
      Path <- list()
      Path$St <- c(BgnSt,EndSt)
      Path$Tm <- c(0,Tm)
      Path$ptm <- c(ptm2-ptm1,ptm3-ptm2)
      return(Path)
    }
    ## if true jump: done
    if (BgnSt!=EndSt){
      Path <- list()
      Path$St <- c(BgnSt,EndSt,EndSt)
      Path$Tm <- c(0,Tm*runif(1),Tm)
      Path$ptm <- c(ptm2-ptm1,ptm3-ptm2)
      return(Path)
    }
  }
  ## Case (nJmp >= 2):
  ## Simulate jumping times
  JmpTmV <- Tm*sort(runif(nJmp))
  ## Simulate states (last state always EndSt)
  JmpStV <- rep(0,(nJmp-1))
  Prb1 <- TransM[BgnSt,]
  for (i in 1:(nJmp-1)){
    Prb2Mat <- TransArr[,,(nJmp-i)]
    Prb2 <- Prb2Mat[,EndSt]
    JmpStV[i] <-
      sample(1:nSt,size=1,replace=TRUE,round(Prb1*Prb2,digits=10))
    Prb1 <- TransM[JmpStV[i],]
  }
  ptm3 <- proc.time()[1]
  JmpStV <- c(JmpStV,EndSt)
  ## Remove virtual substitutions
  TrueSub <- c(BgnSt,JmpStV[1:nJmp])!=c(JmpStV[1:nJmp],EndSt)
  State.sequence<-c(BgnSt,JmpStV[TrueSub],EndSt)
  Jump.times <-c(0,JmpTmV[TrueSub],Tm)
  
  Path <- list()
  Path$St <- State.sequence
  Path$Tm <- Jump.times
  return(Path)
  
}

Uniformization<-function(BgnSt,EndS,RateM,Tm){
  out<-UniformSampl(BgnSt,EndS,RateM,Tm)
  EndPath<-rbind(unlist(out$St)[-length(unlist(out$St))],c(diff(unlist(out$Tm))))
  if(is.null(dim(EndPath))) dim(EndPath)<-c(2,1)
  rownames(EndPath)<-c("states","times")
  return(EndPath)
}


####################################################################################
####################################################################################
####################################################################################



PathSampling<-function(data,Q,P,pi,psi,gamma,eta,intervals,grid,S){
  if(length(unlist(data$history))==0){
    #### Assumption: State 1 has lower intensity ####
    u.new<-matrix(c(1,10,10),ncol = 3)
  }else{
  TT<-data$TT
  ### Forward Backward ###
  z<-unlist(ForwardBackward(data,P,pi,psi,gamma,eta,intervals,grid)$s)  #discrete trajectory points 
  u<-matrix(unlist(sapply(2:length(z),function(j) Uniformization(BgnSt=z[j-1],EndS=z[j],RateM=Q,Tm=TT/intervals))),ncol = 2,byrow = TRUE)
  u[,2]<-cumsum(u[,2])
  u<-u[!duplicated(u[,2]),]
  u.new<-matrix(u[1,],nrow = 1)
  for(i in 2:dim(u)[1]){
    if(u[i,1]==u[i-1,1]){
      u.new[dim(u.new)[1],2]<-u[i,2]
    }else{
      u.new<-rbind(u.new,u[i,])
    }
  }
  u.new<-cbind(u.new,diff(c(0,u.new[,2])))
  }
  colnames(u.new)<-c("states","jump times","sojourn times")
  return(u.new)
}

ForwardBackward<-function(data,P,pi,psi,gamma,eta,intervals,grid){
  S<-dim(P)[1]
  individual.history<-unlist(data$history)
  Ti<-length(individual.history)
  Fj <- Fju <- matrix(NA,nrow = intervals+1,ncol = S)
  Fj[1,]<- Fju[1,] <- log(pi)
  for(j in 2:length(grid)){
    y<-individual.history[individual.history>grid[j-1]&individual.history<grid[j]]
    if(length(y)>0) dt<- sapply(1:S, function(s) sum(sapply(1:length(y),function(i) dTIP.event(data,y[i],eta,gamma,psi[,s]))))
    if(length(y)==0) dt<- rep(0,S) #(Zucchini et al.,2016)
  #Fju[j,]<-apply(log(t(P))+Fju[j-1,]+dt,1,sumlog)
  #Fj[j,]<- Fju[j,]-sumlog(Fju[j,])
    Fju[j,]<-log(t(P)%*%exp(Fju[j-1,])*exp(dt))
    Fj[j,]<- Fju[j,]-log(sum(exp(Fju[j,])))
  }
  ### Backward simulations of the Chain ###
  Fj=exp(Fj)
  s<-c()
  s[length(grid)]<-sample(1:S,size = 1,prob = Fj[length(grid),]) 
  for(t in intervals:1) {
    p<-Fj[t,]*P[,s[t+1]]
    s[t]<-sample(1:S,size = 1,prob =p/sum(p))
  }
  out<-list()
  out$s<-s
  return(out)
}


####################################################################################
####################################################################################
####################################################################################

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

sumlog<-function(x) return(log(sum(exp(x))))



permute <- function(items){
  n <- length(items)
  if (n == 0) {
    return(list())
  } else if (n == 1) {
    return(list(items))
  } else {
    result <- list()
    for (i in seq_along(items)) {
      item <- items[i]
      remaining_items <- items[-i]
      sub_permutations <- permute(remaining_items)
      for (j in seq_along(sub_permutations)) {
        sub_permutation <- sub_permutations[[j]]
        permutation <- c(item, sub_permutation)
        result <- append(result, list(permutation))
      }
    }
    return(result)
  }
}

#### Function for Adaptive Metropolis-Hastings 
calibrate<-function(rho,constant){
  if(rho>0.5) constant<-constant+0.5 
  if(rho>0.3&rho<0.5) constant<-constant+0.1
  if(rho<0.2&rho>0.15) constant<-constant-0.1 
  if(rho<0.15) constant<-constant-0.5 
  return(abs(constant))
}


HMM_parametric_MCMCsamplingTIPs<-function(data,MCMCiter,MCMC.strategy,sigma.prop=NULL,time.partition=10,prior.p=NULL,burn.in=0.15,mu.starting=c(-0.4,1.2),theta.starting=c(0.3,0.5),alpha.starting=c(0.5,3),beta.starting=c(1.4,4),gamma.starting=0.4,eta.starting=1,S,initial.distribution="Sampled"){
  
  TT=data[[1]]$TT
  if(S==2){
    if(is.null(sigma.prop)){
      sigma.prop<-list()
      sigma.prop[[1]]<-list()
      sigma.prop[[1]][[1]]<-diag(c(0.1,0.1,0.05,0.05))
      sigma.prop[[1]][[2]]<-diag(c(0.15,0.15,0.1,0.1))
      sigma.prop[[2]]<-0.05
      sigma.prop[[3]]<-0.05
    }
  }
  if(S==3){
    if(is.null(sigma.prop)){
      sigma.prop<-list()
      sigma.prop[[1]]<-list()
      sigma.prop[[1]][[1]]<-diag(c(0.05,0.05,0.05,0.05))
      sigma.prop[[1]][[2]]<-diag(c(0.1,0.1,0.1,0.1))
      sigma.prop[[1]][[3]]<-diag(c(0.2,0.2,0.2,0.2))
      sigma.prop[[2]]<-0.05
      sigma.prop[[3]]<-0.05
    }
  }
  if(is.null(prior.p)) {prior.p<-matrix(1,S,S-1); diag(prior.p)<-1}
  library(mvtnorm)
  library(extraDistr)
  library(truncnorm)
  #library(snipEM)
  library(expm)
  library(castor)
  intervals = time.partition
  n = length(data)  
  #### Set the algorithm ####
  gamma<-rep(0,MCMCiter)
  eta<-rep(0,MCMCiter)
  liks <- matrix(0,length(data),MCMCiter)
  gamma[1]=0
  psi<-array(0,dim = c(MCMCiter,4,S))
  #Starting points
  psi[1,1,]=seq(min(mu.starting),max(mu.starting),length=S)
  psi[1,2,]=seq(min(theta.starting),max(theta.starting),length=S)
  psi[1,3,]=seq(min(alpha.starting),max(alpha.starting),length=S)
  psi[1,4,]=seq(min(beta.starting),max(beta.starting),length=S)
  gamma[1]=gamma.starting
  eta[1]=eta.starting
  #Hidden Markov Transition Probability Matrix
  Q<-array(0,dim = c(MCMCiter,S,S)) ;  Q[1,,]<-0.1; for(s in 1:S) Q[1,s,s]=-sum(Q[1,s,-s])
  Prs<-array(0,dim = c(MCMCiter,S,S-1))
  #Updating strategies for the marginal posterior distributions
  Updating<-MCMC.strategy
  if(Updating=="Joint"){ Sigma.psi<-sigma.prop[[1]] ; sigma.gamma<-sigma.prop[[2]]; sigma.eta<-sigma.prop[[3]] 
  ACCEPTANCE<-list()
  for(s in 1:S) ACCEPTANCE[[s]]<-0
  constant<-rep(2.38,S)
  }
  if(Updating=="Sequential"){ Sigma.psi<-sigma.prop[[1]] ; sigma.gamma<-sigma.prop[[2]]; sigma.eta<-sigma.prop[[3]] }
  P<-expm((TT/intervals)*Q[1,,])
  if(initial.distribution=="Stationary"){
    pi<-get_stationary_distribution(Q[1,,]) 
  }
  if(initial.distribution=="Sampled"){
    pi_0<-matrix(0,nrow = MCMCiter, ncol = S)
    pi_0[1,]<-1/dim(pi_0)[2]
    pi<-pi_0[1,]
  }
  
  grid<-seq(0,TT,length.out=intervals+1)
  permutations<-permute(1:S)
  par(mar=c(1,1,1,1))
  print("starting MCMC") 
  for(iter in 2:MCMCiter){
    if(iter%%10==0) {print(iter)}
    
    if(iter%%100==0) {
      if(iter>500) burnin=500 else burnin=1
      par(mfrow=c(3,S))
      for(s in 1:S){
        plot(psi[burnin:(iter-1),3,s], type = "l", xlab = "", ylab = "", main="alpha")
        plot(psi[burnin:(iter-1),4,s], type = "l", xlab = "", ylab = "", main="beta")
        plot(-Q[burnin:(iter-1),s,s], type = "l", xlab = "", ylab = "", main="q")
      } 
      # par(mfrow=c(2,1))
      psi.post.mean<-t(apply(psi[burnin:(iter-1),,],2:3,mean))
      colnames(psi.post.mean)<-c("mu","theta","alpha","beta")
      print(psi.post.mean)
      print(apply(Q[(iter-100):(iter-1),,], 2:3, mean))
      if(initial.distribution=="Sampled") print(apply(pi_0[(iter-100):(iter-1),],2,mean))
    }
    
    
    #UPDATE PARAMETERS CONDITIONAL TO THE LATENT STRUCTURE
    
     repeat{
       try(Z<-lapply(1:n,function(j) PathSampling(data[[j]],Q[iter-1,,],P,pi,psi[iter-1,,],gamma[iter-1],eta[iter-1],intervals,grid,S)), silent = TRUE)
       if(is.list(Z)) break
      }
    
   # Z<-lapply(1:n,function(j) PathSampling(data[[j]],Q[iter-1,,],P,pi,psi[iter-1,,],gamma[iter-1],eta[iter-1],intervals,grid,S))
    
    for(i in 1:n) data[[i]]$U<-Z[[i]]
    ###Avoid Label switching by using model selection criteria
    modulated.parameters<-psi[iter-1,,]
    lik<-lapply(1:length(permutations), function(i) apply(sapply(1:S, function(j) sapply(data, function(data) ind.cond.lik(data,st=j,eta[iter-1],gamma[iter-1],modulated.parameters[,permutations[[i]]],log=TRUE) )),1,sum) )
    probs<-sapply(1:length(permutations), function(j) sum(lik[[j]])) 
    probs<- exp(probs-sumlog(probs))
    if(any(is.infinite(probs)) | any(is.nan(probs))){
      probs<-sapply(1:length(permutations), function(j) sum(lik[[j]])) 
      probs<-exp(probs)
      }
    if(sum(probs)==0) probs[]<-1
    perm<-sample(x=1:length(permutations),size = 1,prob = probs/sum(probs) )
    if(perm!=1){
      #Relabel states
      for(i in 1:n){
        Z<-data[[i]]$U
        for (j in 1:S){
          Z[which(Z[,1]==j),1]=permutations[[perm]][[j]]  #Forse si può scrivere meglio
        }
        data[[i]]$U<-Z 
      }
    }
    
    liks[,iter-1]<-sapply(1:n,function(i) sum(sapply(1:S, function(j) ind.cond.lik(data[[i]],j,eta[iter-1],gamma[iter-1],psi[iter-1,,],log = TRUE))) )
    
    
    ### Calculate the transitons ###
    N<-matrix(0, nrow=S, ncol=S)
    Sojourn<-NULL
    Visits<-NULL
    last.state<-c()
    starting.state<-NULL
    for(i in 1:n){
      Z<-data[[i]]$U
      states<-Z[,1]
      starting.state<-c(starting.state,states)
      Chain<-matrix(Z[,-2],ncol = 2)
      last.state[i]<-Chain[dim(Chain)[1],1]
      N.star = array(0, dim = c(S,S,length(states)-1))
      Sojourn<-cbind(Sojourn,sapply(1:S,function(j) sum(Chain[which(Chain[,1]==j),2])))
      Visits<-cbind(Visits,sapply(1:S,function(j) length(Chain[which(Chain[,1]==j),1])))
      for (h in 1: length(states)-1){
        p = states[h]
        j = states[h+1]
        N.star[p,j,h] = 1
      }
      rows <- dim(N.star)[1]
      cols <- dim(N.star)[2]
      transitions <- matrix(0, nrow=rows, ncol=cols)
      for (m in seq(rows)) {
        for (g in seq(cols)) {
          transitions[m,g] <- sum(N.star[m,g,])
        }
      }
      N<-N+transitions
    }
    Sojourn<-apply(Sojourn, 1, sum)
    Visits<-apply(Visits, 1, sum)
    #### VERY INFORMATIVE PRIOR ON THE SOJOURN TIMES RATES
    diag(Q[iter,,])<-sapply(1:S, function(j) rgamma(1,shape = 1 + Visits[j] - sum(last.state==j) ,rate= 1 + Sojourn[j] ) )
    ### Update the transition rate matrix ###
    if(S==2) Q[iter,1,2]=Q[iter,1,1]; Q[iter,2,1]=Q[iter,2,2]
    if(S>2){
      Prs[iter,,]<-matrix(unlist(lapply(1:S,function(s){rdirichlet(1,prior.p[s,]+N[s,-s])})),nrow = S, byrow = TRUE)
      for (s in 1:S) {
        Q[iter,s,-s]<-Prs[iter,s,]*Q[iter,s,s]
      }
    }
    diag(Q[iter,,])=-diag(Q[iter,,])
    ### Generate the transition probability matrix for discrete Chain simulation ###
    P<-expm((TT/intervals)*Q[iter,,])
    #Eigen <- eigen(P)
    
      if(initial.distribution=="Stationary"){
        pi<-get_stationary_distribution(Q[iter,,]) #Re(Eigen$vectors[1,]) ### If is complex take the real part
        if(any(pi<0)) pi[pi<0]=0
      }
    
    
      if(initial.distribution=="Sampled"){
        start.points<-sapply(1:S,function(i) length(which(starting.state==i)))
        pi_0[iter,]<-rdirichlet(1,alpha = 1 + start.points)
        pi<-pi_0[iter,]
      }
    
    
    psi[iter,,]<-psi[iter-1,,]
    
    ##Update Psi|U # Joint Updating
    if(Updating=="Joint"){
      for(s in 1:S){
        last.psi<-psi[iter,,]
        if(iter>200){
        if(iter>300&iter%%100==0){
          constant[s]<- calibrate(rho=mean(ACCEPTANCE[[s]][-(1:(iter-100))]),constant = constant[s])
        }
        Sigma<-((constant[s]^2)/dim(psi)[2])*cov(psi[2:(iter-1),,s])
        }else{
        Sigma<-Sigma.psi[[s]]*0.5}
        proposal<-rmvnorm(1,c(last.psi[1,s],log(last.psi[2:4,s])),Sigma)
        proposal[2:4]<-exp(proposal[2:4])
        if(proposal[3]<proposal[4]){ 
          psi.prop<-last.psi
          psi.prop[,s]<-proposal
          lnum<-sum(sapply(data, function(data) ind.cond.lik(data,st=s,eta[iter-1],gamma[iter-1],psi.prop,log=TRUE)))
          lden<-sum(sapply(data, function(data) ind.cond.lik(data,st=s,eta[iter-1],gamma[iter-1],last.psi,log=TRUE)))
          lprior.num<-dnorm(psi.prop[1,s],mean = 0,sd=1000,log=TRUE)+dgamma(psi.prop[2,s],1,100,log = TRUE) +dunif(psi.prop[3,s],0,psi.prop[4,s],log = TRUE) +dunif(psi.prop[4],0,15,log = TRUE)
          lprior.den<-dnorm(last.psi[1,s],mean = 0,sd=1000,log=TRUE)+dgamma(last.psi[2,s],1,100,log = TRUE) +dunif(last.psi[3,s],0,last.psi[4,s],log = TRUE) + dunif(last.psi[4],0,15,log = TRUE)
          ACCEPT<-exp(sum(lnum)+lprior.num-sum(lden)-lprior.den)
          #### FIXING NUMERICAL PROBLEMS ######
          if(is.na(ACCEPT)|is.nan(ACCEPT)) ACCEPT=0
          if(ACCEPT>runif(1)){
            psi[iter,,]=psi.prop  
            ACCEPTANCE[[s]]<-c(ACCEPTANCE[[s]],1)
          }else{ 
           ACCEPTANCE[[s]]<-c(ACCEPTANCE[[s]],0)
          }
        }else{  ACCEPTANCE[[s]]<-c(ACCEPTANCE[[s]],0) }
      }
    }

    
    if(Updating=="Sequential"){
      for(s in 1:S){
        for (j in 1:2){
          last.psi<-psi[iter,,]
          if(iter>200) Sigma.psi[[s]][j,j]<- sqrt(2.38^2*var(psi[1:(iter-1),j,s]))
          if(j==1)  proposal<-rnorm(1,last.psi[j,s],sd=Sigma.psi[[s]][j,j] )
          if(j==2)  proposal<-exp(rnorm(1,log(last.psi[j,s]),Sigma.psi[[s]][j,j]))
          psi.prop<-psi[iter,,]
          psi.prop[j,s]<-proposal
          lnum<-sum(sapply(data, function(data) ind.cond.lik(data,st=s,eta[iter-1],gamma[iter-1],psi.prop,log=TRUE)))
          lden<-sum(sapply(data, function(data) ind.cond.lik(data,st=s,eta[iter-1],gamma[iter-1],last.psi,log=TRUE)))
          if(j==1){ 
            lprior.num<-dnorm(psi.prop[j,s],mean = 0,sd=10,log = TRUE)
            lprior.den<-dnorm(last.psi[j,s],mean = 0,sd=10,log = TRUE)
          }
          if(j==2){
            lprior.num<-dgamma(psi.prop[j,s],1,1,log = TRUE)
            lprior.den<-dgamma(last.psi[j,s],1,1,log = TRUE)
          }
          ACCEPT<-exp(sum(lnum)+lprior.num-sum(lden)-lprior.den)
          if(ACCEPT>runif(1)) psi[iter,j,s]=proposal 
        }
        #### Update Hawkes part ####
         j=3:4
    #    if(iter>200){Sigma.psi[[s]][j[1],j[1]]<- sqrt(2.38^2*var(psi[1:(iter-1),j[1],s])) ; Sigma.psi[[s]][j[2],j[2]]<- sqrt(2.38^2*var(psi[1:(iter-1),j[2],s])) }
        if(iter>200){Sigma.psi[[s]][j[1],j[1]]<- sqrt(2.38^2*var(psi[1:(iter-1),j[1],s])) ; Sigma.psi[[s]][j[2],j[2]]<- sqrt(2.38^2*var(psi[1:(iter-1),j[2],s])) }
          last.psi<-psi[iter,,]
        proposal<-exp(rnorm(1,log(last.psi[j[1],s]),sd=Sigma.psi[[s]][j,j][1,1]))
        proposal<-c(proposal,rtruncnorm(1,a=proposal,b=Inf,mean = last.psi[j[2],s], sd=Sigma.psi[[s]][j,j][2,2]))
        psi.prop<-psi[iter,,]
        psi.prop[j,s]<-proposal
        #lprior.num<-dunif(psi.prop[j[1],s],0,psi.prop[j[2],s],log = TRUE) + log(dtruncnorm(psi.prop[j[2],s],a=0,b=Inf,mean=2,sd=1))
        #lprior.den<-dunif(last.psi[j[1],s],0,last.psi[j[2],s],log = TRUE) + log(dtruncnorm(last.psi[j[2],s],a=0,b=Inf,mean=2,sd=1))
        lprior.num<-dunif(psi.prop[j[1],s],0,psi.prop[j[2],s],log = TRUE) + dunif(psi.prop[j[2],s],0,100,log = TRUE)
        lprior.den<-dunif(last.psi[j[1],s],0,last.psi[j[2],s],log = TRUE) + dunif(last.psi[j[2],s],0,100,log = TRUE)
        lnum<-sum(sapply(data, function(data) ind.cond.lik(data,st=s,eta[iter-1],gamma[iter-1],psi.prop,log=TRUE) ))
        lden<-sum(sapply(data, function(data) ind.cond.lik(data,st=s,eta[iter-1],gamma[iter-1],last.psi,log=TRUE) ))
        ACCEPT<-exp(sum(lnum)+lprior.num-sum(lden)-lprior.den)
        if(ACCEPT>runif(1)) psi[iter,j,s]=proposal 
      }
    }
    
    #Update gamma
    last.gamma<-gamma[iter-1]
    gamma.prop<-rnorm(1,last.gamma,sigma.gamma)
    lnum<-sapply(1:S, function(s) sum(sapply(data, function(data) ind.cond.lik(data,st=s,eta[iter-1],gamma.prop,psi[iter,,],log=TRUE))))
    lden<-sapply(1:S, function(s) sum(sapply(data, function(data) ind.cond.lik(data,st=s,eta[iter-1],last.gamma,psi[iter,,],log=TRUE))))
    ACCEPT<-exp(sum(lnum)+dnorm(gamma.prop,0,10,log = TRUE)-sum(lden)-dnorm(last.gamma,0,10,log=TRUE)); 
    #### FIXING NUMERICAL PROBLEMS ######
    if(is.na(ACCEPT)|is.nan(ACCEPT)) ACCEPT=0
    if(ACCEPT>runif(1)) gamma[iter]=gamma.prop else gamma[iter]=last.gamma
    
    # Update the baseline 
    last.eta<-eta[iter-1]
    eta.prop<-exp(rnorm(1,log(last.eta),sigma.eta))
    lnum<-sapply(1:S, function(s) sum(sapply(data, function(data) ind.cond.lik(data,st=s,eta.prop,gamma[iter],psi[iter,,],log=TRUE))))
    lden<-sapply(1:S, function(s) sum(sapply(data, function(data) ind.cond.lik(data,st=s,last.eta,gamma[iter],psi[iter,,],log=TRUE))))
    ACCEPT<-exp(sum(lnum)+dnorm(log(eta.prop),0,10,log = TRUE)-sum(lden)-dnorm(log(last.eta),0,10,log=TRUE)); 
    #### FIXING NUMERICAL PROBLEMS ######
    if(is.na(ACCEPT)|is.nan(ACCEPT)) ACCEPT=0
    if(ACCEPT>runif(1)) eta[iter]=eta.prop else eta[iter]=last.eta
  }
  
  burnin<-round(MCMCiter*burn.in)
  # Posterior summary statistics #
  out<-list() 
  out$psi<-psi
  out$gamma<-gamma
  out$eta<-eta
  liks<-liks[,-(1:burnin)]
  sample.size<-dim(liks)[2]
  out$WAIC <- -2*(sum(apply(liks,1,sumlog)-log(1/sample.size))-sum(apply(liks,1,var)))
  #####
  quant <- function(x) c(quantile(x,0.025),quantile(x,0.975))
  psi.posterior<-list()
  psi.post.mean<-t(apply(psi[burnin:(iter-1),,],2:3,mean))
  colnames(psi.post.mean)<-c("mu","theta","alpha","beta")
  psi.post.sd<-t(apply(psi[burnin:(iter-1),,],2:3,sd))
  colnames(psi.post.sd)<-c("mu","theta","alpha","beta")
  c.i.psi<-apply(psi[-(1:burnin),,], 2:3, quant)
  colnames(c.i.psi)<-c("mu","theta","alpha","beta")
  psi.posterior$post.mean<-psi.post.mean
  psi.posterior$post.sd<-psi.post.sd
  psi.posterior$cred.int<-c.i.psi
  gamma.posterior<-list()
  gamma.posterior$post.mean<-mean(gamma[-(1:burnin)])
  gamma.posterior$post.sd<-sd(gamma[-(1:burnin)])
  gamma.posterior$cred.int<-quant(gamma[-(1:burnin)])
  eta.posterior<-list()
  eta.posterior$post.mean<-mean(eta[-(1:burnin)])
  eta.posterior$post.sd<-sd(eta[-(1:burnin)])
  eta.posterior$cred.int<-quant(eta[-(1:burnin)])
  pi_0.posterior<-list()
  pi_0.posterior$post.mean<-apply(pi_0,2,mean)
  pi_0.posterior$post.sd<-apply(pi_0,2,sd)
  #Transition Probability Posterior summary
  if(S>2){
    P.posterior<-list()
    P.posterior$post.mean<-apply(Prs, 2:3, mean)
    P.posterior$post.sd<-apply(Prs, 2:3, sd)
    P.posterior$cred.int<-apply(Prs, 2:3, quant)
    out$P.posterior<-P.posterior
  }
  #Transition rate Posterior summary
  Q.posterior<-list()
  Q.posterior$post.mean<-apply(Q, 2:3, mean)
  Q.posterior$post.sd<-apply(Q, 2:3, sd)
  Q.posterior$cred.int<-apply(Q, 2:3, quant)
  
  out$psi.posterior<-psi.posterior
  out$gamma.posterior<-gamma.posterior
  out$Q.posterior<-Q.posterior
  out$pi_0<-pi_0.posterior
  return(out)
}

# HMM_parametric_MCMCsamplingTIPs=cmpfun(HMM_parametric_MCMCsamplingTIPs)
# simParHM<-HMM_parametric_MCMCsamplingTIPs(data=data,MCMC.strategy="Joint",MCMCiter=5000,S=3)

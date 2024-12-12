# FUNZIONI

#1) lambda(t)
#2) integrale di lambda in [0,T]
#3) integrale di lambda con integrate (per check)
#4) integrale di lambda per individui/tempi con zero catture
#5) likelihood per singola unita'
#6) opposto della likelihood per tutte le unita'
#7) opposto della likelihood condizionata per tutte le unità

#0) funzione integranda per quadratura Chebyshev
library(statmod)
dd=gauss.quad(39,"chebyshev1")
fint=function(x,eta,theta,lower,upper) {
  cv=(upper-lower)*x/2+(lower+upper)/2
  eta*cv^(eta-1)*exp(theta*cv^(eta))*sqrt(1-x^2)*(upper-lower)/2}

####################################################################
#1) lambda(t)

lambda=function(alpha,beta,eta,theta,gamm,mu,XX,individual.hist,tt){ 
  lam=alpha*sum(exp(-beta*(tt-individual.hist[individual.hist<tt])))+
    eta*(tt^(eta-1))*exp((XX%*%gamm)+mu)*exp(-theta*(sum(individual.hist<tt)-(tt^(eta))))
  return(lam)
}

####################################################################
#2) integrale di lambda in [0,T]

lambda.integral=function(alpha,beta,eta,theta,gamm,mu,XX,individual.hist,TT){
  #SS=length(individual.hist)
  #parte self-exciting
  somma.sep=-alpha/beta*sum(exp(-beta*(TT-individual.hist))-1)
  #parte self-correcting
  fuori=exp((XX%*%gamm)+mu)#*eta  ##as.numeric?
  tt2=c(0,individual.hist,TT)
  vett.int=(sapply(1:(length(tt2)-1), function(j) sum(fint(dd$nodes,eta,theta,tt2[j],tt2[j+1])*dd$weights))) ##as.numeric?
  somma.scp=fuori*sum(exp(-theta*0:length(individual.hist))*vett.int)
  return(as.numeric(somma.sep+somma.scp))
}

####################################################################
#4) integrale di lambda per individui/tempi con zero catture

lambda.integral.0capture=function(eta,theta.v,gamm,mu.v,XX,TT,WW){
  XX=as.matrix(XX)
  fuori=sapply(1:length(mu.v), function(j) exp((XX%*%gamm)+mu.v[j]))#*eta ##as.numeric?
  integr=sapply(1:length(theta.v), 
                function(j) sum(dd$weights*fint(dd$nodes,eta,theta.v[j],0,TT)))
  if(!all(dim(fuori)==dim(WW))) warning("Problem with function lambda integral")
  out=as.numeric((fuori*WW)%*%integr) #as.numeric?
  return(out)
}

####################################################################
#4b) integrale di lambda per individui/tempi con zero catture come matrice

lambda.integral.0capture.mat=function(eta,theta.v,gamm,mu.v,XX,TT,CC){
  XX=as.matrix(XX)
  fuori=sapply(1:CC, function(j) exp((XX%*%gamm)+mu.v[j]))#*eta ##as.numeric?
  integr=sapply(1:CC, function(j) sum(fint(dd$nodes,eta,theta.v[j],0,TT)*dd$weights))
  t(integr*t(matrix(fuori,ncol=CC)))
}


####################################################################
#5a) likelihood condizionata a classe c per singola unita'

lik.unit.condC=function(alpha,beta,eta,theta,gamm,mu,XX,individual.hist,TT){
  lik=0
  #if(length(c(individual.hist,TT))>1){ #se c'è almeno una cattura LO ASSUMO!!
  for(i in 1:length(individual.hist))
    lik=lik+log(lambda(alpha,beta,eta,theta,gamm,mu,XX,individual.hist[1:i],individual.hist[i]))
  
  lik=lik-lambda.integral(alpha,beta,eta,theta,gamm,mu,XX,individual.hist,TT)
  # } else { #se non ci sono catture
  #lik=-lambda.integral.0capture(eta,theta,gamm,mu,XX,TT,WW) #NB se usi questa riga devi aggiungere WW tra gli argomenti
  #}
  return((lik)) #as.numeric?
}

####################################################################
#5b) likelihood completa per singola unita'

lik.unit.compl=function(alpha.v,beta.v,eta,theta.v,gamm,mu.v,XX,individual.hist,TT,WW,pi.C){
  lik.v=(sapply(1:length(alpha.v), function(j) lik.unit.condC(alpha.v[j],beta.v[j],eta,theta.v[j],
                                                              gamm,mu.v[j],XX,individual.hist,TT)))
  lik.v*(WW)+(WW)*log(pi.C) #as.numeric(WW)?
}

####################################################################
#6) opposto della likelihood per tutte le unita'

#obs.data è una lista, ogni suo elemento è un individuo con tempi di cattura 
#restituisce meno la verosimiglianza per tutto il dataset
lik.all=function(alpha.v,beta.v,eta,theta.v,gamm,mu.v,XX,obs.data,TT,WW,pi.C){
  sum(unlist(lapply(1:length(obs.data), function(j) -lik.unit.compl(alpha.v,beta.v,eta,
                                                                    theta.v,gamm,mu.v,XX[j,],unlist(obs.data[[j]]),TT,WW[j,],pi.C))))}

####################################################################
#7) opposto della likelihood condizionata per tutte le unità

#obs.data è una lista, ogni suo elemento è un individuo con tempi di cattura 
cond.lik.all=function(alpha.v,beta.v,eta,theta.v,gamm,mu.v,XX,obs.data,TT,WW,pi.C) {
  lik.all(alpha.v,beta.v,eta,theta.v,gamm,mu.v,XX,obs.data,TT,WW,pi.C)+
    log(1-exp(-lambda.integral.0capture(eta,theta.v,gamm,mu.v,XX,TT,WW)))
}

#qui c'è qualcosa da sistemare, dobbiamo assicurarci che a seconda di o, te, b funzioni tutto
#se modello completo, par.vec lungo CCx4 + pp + 1 è il vettore che ha in ordine c(alpha.v, beta.v, eta, theta.v, gamm, mu.v)
#b=F no behavioural -> no alpha, no beta, no theta (as many as CC), par.vec lungo CC + pp + 1 cioè c(eta, gamm, mu.v)
#o=F no covariates -> no gamma (as many as ncol(XX)), par.vec c(alpha.v, beta.v, eta, theta.v, mu.v)
#te=F no time effects -> no eta (1 parameter), par.vec c(alpha.v, beta.v, theta.v, gamm, mu.v)
#
wrap.lik=function(par.vec, XX, pp, obs.data, TT, CC, WW, pi.C, o, te, b) {
  if(b) alpha.v=exp(par.vec[1:CC]) else alpha.v=rep(0,CC)
  if(b) beta.v=exp(par.vec[CC+1:CC])+alpha.v else beta.v=rep(1,CC)
  if(b & te) eta=exp(par.vec[2*CC+1]) else if(!te) eta=1 else if(!b & te) eta=exp(par.vec[1])
  if(b & te) theta.v=exp(par.vec[(2*CC+1)+1:CC]) else if(b & !te) theta.v=exp(par.vec[(2*CC)+1:CC]) else if(!b) theta.v=rep(0,CC)
  if(!o) gamm=rep(0, pp) else gamm=par.vec[(length(par.vec)-(CC-1))-pp:1]
  mu.v=par.vec[(length(par.vec)-(CC-1)):length(par.vec)]
  
  
  sum(unlist(lapply(1:length(obs.data), function(j) 
    -lik.unit.compl(alpha.v,beta.v,eta, theta.v,gamm,mu.v,matrix(XX[j,], nrow = 1),obs.data[[j]],
                    TT,matrix(WW[j,], nrow = 1),pi.C))))+
    sum((log(1-exp(-lambda.integral.0capture(eta,theta.v,gamm,mu.v,XX,TT,
                    WW)[!is.na(lambda.integral.0capture(eta,theta.v,gamm,mu.v,XX,TT,WW))])+1e-30)))
}



#####################################################################
#ALTRE FUNZIONI PER TIP

safe_sumlog = function(x) {

s <- (try(snipEM::sumlog(x), silent=T))
if(inherits(s,"try-error")) {
s <- max(x)}

return(s)}

unconstrain = function(alpha.v, beta.v, eta, theta.v, gamm, mu.v, o, te, b) {
  par.vec = c()

  if (b) par.vec = c(par.vec, log(alpha.v), log(beta.v - alpha.v))
  if (te) par.vec = c(par.vec, log(eta))
  if (b) par.vec = c(par.vec, log(theta.v))
  if (o) par.vec = c(par.vec, gamm)

  par.vec = c(par.vec, mu.v)

  return(par.vec)
}

unconstrainPiC = function(alpha.v, beta.v, eta, theta.v, gamm, mu.v, piC, o, te, b) {
  par.vec = c()

  if (b) par.vec = c(par.vec, log(alpha.v), log(beta.v - alpha.v))
  if (te) par.vec = c(par.vec, log(eta))
  if (b) par.vec = c(par.vec, log(theta.v))
  if (o) par.vec = c(par.vec, gamm)

  par.vec = c(par.vec, mu.v, log(piC[-1] / piC[1]))

  return(par.vec)
}

stDirectVec=function(par.vec,CC,XX,TT,obsData,b,te,o,inits) {

   pp = ncol(XX)
   n  = nrow(XX)
    i = 0

    if(!b) {alpha=inits$alpha
            beta=inits$beta
	    theta=inits$theta}

    if(!te) {eta = inits$eta}

    if(!o) {gamm = inits$gamm}

    if (b) {
      alpha = exp(par.vec[i + 1:CC])
      i = i + CC
      beta = exp(par.vec[i + 1:CC]) + alpha
      i = i + CC
    }

    if (te) {
      eta = exp(par.vec[i+1])
      i = i + 1
    }

    if(b) {
      theta = exp(par.vec[i + 1:CC])
      i = i + CC
    }

    if (o) {
      gamm = par.vec[i + 1:pp]
      i = i + pp
    }

    mu = par.vec[i + 1:CC]
    i = i + CC

    piC = c(1, exp(par.vec[i + 1:(CC - 1)]))
    piC = piC / sum(piC)

    or = order(mu)
    piC = matrix(piC[or], n, CC, byrow = TRUE)
    mu = mu[or]
    alpha = alpha[or]
    beta = beta[or]
    theta = theta[or]

obj=list()
obj$alpha=alpha[or]
obj$beta=beta[or]
obj$eta=eta
obj$theta=theta[or]
obj$gamm=gamm
obj$mu=mu[or]
obj$piC=piC

stDirect(obj, XX, TT, CC, obsData)[1]}

stDirect = function(object, XX, TT, CC, obsData) {
  n = nrow(XX)
  XX = as.matrix(XX)
  pp = ncol(XX)
  ex = expand.grid(1:CC, 1:n)
  WW = matrix(1 / CC, n, CC)
  eta = object$eta
  piC = object$piC
  theta.v = object$theta
  gamm = object$gamm
  mu.v = object$mu
  alpha.v = object$alpha
  beta.v = object$beta
  ma = 1 - exp(-lambda.integral.0capture.mat(eta,theta.v,gamm,mu.v,XX,TT,CC))

  ap = matrix(apply(ex, 1, function(jj) lik.unit.condC(alpha.v[jj[1]], beta.v[jj[1]], eta, theta.v[jj[1]], gamm, mu.v[jj[1]], XX[jj[2], ],
                                                       obsData[[jj[2]]], TT)), nrow = n, byrow = TRUE)

  sa = -log(1 - exp(t(t(-lambda.integral.0capture.mat(eta, theta.v, gamm, mu.v, XX, TT, CC)))))

  if (CC > 1)
    sa = sa + matrix(log(piC), n, CC, byrow = TRUE)

  if (CC == 1)
    sa = sa + log(piC)

  condlik.mat = ap + sa

  # pesi delle classi per ogni individuo
  obj = apply(condlik.mat, 1, safe_sumlog)
  log.WW = condlik.mat - obj
  WW = exp(log.WW)
  wm = apply(WW, 1, which.max)
  Nhat = sum(1 / diag(ma[, wm]))
  vett.condiz = sapply(1:nrow(XX), function(j) 1 - exp(-lambda.integral.0capture(eta, theta.v, gamm, mu.v, matrix(XX[j, ], nrow = 1), TT, matrix(WW[j, ], nrow = 1))))

  Nhat2 = sum(1 / vett.condiz)

  return(c(Nhat, Nhat2))
}


objf=function(par.vec, CC, o, b, te, xxobs, n, TT, obs_data, ex) {
    pp = ncol(xxobs) 
                i = 0
        if (b) {
            alpha.v = exp(par.vec[i + 1:CC])
            i = i + CC
            beta.v = exp(par.vec[i + 1:CC]) + alpha.v
            i = i + CC
        }
        if (te) {
            eta = exp(par.vec[i + 1])
            i = i + 1
        }
        if (b) {
            theta.v = exp(par.vec[i + 1:CC])
            i = i + CC
        }
        if (o) {
            gamm = par.vec[i + 1:pp]
            i = i + pp
        }
        mu.v = par.vec[i + 1:CC]
        i = i + CC
        piC = c(1, exp(par.vec[i + 1:(CC - 1)]))
        piC = piC/sum(piC)
        or = order(mu.v)
        piC = matrix(piC[or], n, CC, byrow = TRUE)
        mu.v = mu.v[or]
        alpha.v = alpha.v[or]
        beta.v = beta.v[or]
        theta.v = theta.v[or]

    ap = matrix(apply(ex, 1, 
                      function(jj) lik.unit.condC(alpha.v[jj[1]], beta.v[jj[1]], eta, theta.v[jj[1]], 
                                                  gamm, mu.v[jj[1]], xxobs[jj[2], ], obs_data[[jj[2]]], TT)), 
                nrow = n, byrow = TRUE)

    sa = -log((1 - exp(t(t(-lambda.integral.0capture.mat(eta, theta.v, gamm, mu.v, xxobs, TT, CC)))))+1e-30)

    if (CC > 1) 
      sa = sa + matrix(log(piC), n, CC, byrow = TRUE)
    condlik.mat = ap + sa
    obj = apply(condlik.mat, 1, safe_sumlog) 
     sum(obj)}

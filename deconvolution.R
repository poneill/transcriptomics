brown <- function(start,time,step){
  "brownian process with exponential killing time"
  if(time <= 0){
    start
  }
  else if (runif(1) < .5){
    brown(start - step, time - 1, step)
  }
  else{
    brown(start + step, time - 1, step)
  }

}
display <- function(x){print(x)
                     x}

rnormlap <- function(mu,sigma.squared,alpha,beta){
  "sample NL distribution"
  mu + sqrt(sigma.squared)*rnorm(1) + rexp(1)/alpha - rexp(1)/beta
}

rnormlap.prime <- function(mu,sigma.squared,alpha,beta){
  mu + sigma.squared*rnorm(1) + rexp(1)/alpha - rexp(1)/beta
}

u <- function(n,mu,sigma.squared,alpha,beta){
  replicate(n,rnormlap(mu,sigma.squared,alpha,beta))
}

prob.vector <- function(m){
  "return a convex combination of m components"
  xs <- runif(m)
  xs/sum(xs)
}

cumulant <- function(a,b,r){
  "compute kr, r > 2, for NL distribution"
  factorial(r-1)*(1/(a^r) + (-1)^r * 1/(b^r))
}

R <- function(z)(1-pnorm(z))/dnorm(z)

dnormlap <- function(y,mu,sigma.squared,alpha,beta){
  "density function for normal-laplace"
  sigma <- sqrt(sigma.squared)
  ans <- alpha*beta/(alpha+beta)*dnorm((y-mu)/sigma)*(R(alpha*sigma-(y-mu)/sigma)+
                                               R(beta*sigma+(y-mu)/sigma))
  ifelse(is.finite(ans),ans,0)
}

log.l <- function(xs,mu,sigma.squared,alpha,beta){
  sum(log(dnormlap(xs,mu,sigma.squared,alpha,beta)))
}

log.l.mixture <- function(xs,mus,sigmas.squared,alphas,betas,pis){
  print(paste("mus:",mus))
  print(paste("sigmas.squared:",sigmas.squared))
  print(paste("alphas:",alphas))
  print(paste("betas:",betas))
  print(paste("pis:",pis))
  f <- function(x) sum(sapply(1:length(pis),function(i) (pis[i] * dnormlap(x,mus[i],
                                                         sigmas.squared[i],
                                                         alphas[i],
                                                         betas[i]))))
  #print(pastef(xs))
  sum(log(f(xs)))
}

foo <- function(y,mu,sigma.squared,alpha,beta){
  sigma <- sqrt(sigma.squared)
  alpha*beta/(alpha+beta)*dnorm((y-mu)/sigma)
}

bar <- function(y,mu,sigma.squared,alpha,beta){
  sigma <- sqrt(sigma.squared)
  R(alpha*sigma-(y-mu)/sigma)
}

baz <- function(y,mu,sigma.squared,alpha,beta){
  sigma <- sqrt(sigma.squared)
  R(beta*sigma+(y-mu)/sigma)
nn}

dnormlap.prime <- function(y,mu,sigma.squared,alpha,beta) foo(y,mu,sigma.squared,alpha,beta)*(bar(y,mu,sigma.squared,alpha,beta) + baz(y,mu,sigma.squared,alpha,beta))

m <- 3
as <- prob.vector(3)
nl1.params <- runif(4)
nl2.params <- runif(4)
nl3.params <- runif(4)

unpack.params <- function(params){
  "accept parameters and return a function for sampling from nl distribution with those parameters"
  mu            <- params[1]
  sigma.squared <- params[2]
  alpha         <- params[3]
  beta          <- params[4]
  function(n){replicate(n,rnormlap(mu,sigma.squared,alpha,beta))}
  }

select.nl <- function(){
  sample(size=1,sapply(list(nl1.params,nl2.params,nl3.params),unpack.params),prob = as)[[1]]
}

moments <- function(xs,ws=replicate(length(xs),1/length(xs))){
  weight <- sum(ws)
  if(is.na(weight) || weight == 0)
    ws=replicate(length(xs),1/length(xs))
  "Compute raw moments of sample xs"
#  print(paste("sum weights: ",sum(ws)))
  ms <- sapply(seq(1,5),function(i)(sum(ws*xs^i)/sum(ws)))
#  print(paste("moments: ",ms))
  ms
}

cumulants <- function(ms){
  "Compute first 5 cumulants, given moments"
  m1 <- ms[1]
  m2 <- ms[2]
  m3 <- ms[3]
  m4 <- ms[4]
  m5 <- ms[5]
  k1 <- m1
  k2 <- m2 - m1^2
  k3 <- 2*m1^3 -3*m1*m2 + m3
  k4 <- -6*m1^4 + 12 * m1^2*m2 - 3 * m2^2 - 4*m1*m3 + m4
  k5 <- 24*m1^5 - 60*m1^3*m2 + 20*m1^2*m3 -10*m2*m3 +5*m1 * (6*m2^2 - m4) + m5
  c(k1,k2,k3,k4,k5)
}

cumulants.from.sample <- function(xs) cumulants(moments(xs))

eq1 <- function(a,b,k3){
  "See eq. 3.12 of Wu 2005"
  k3 - 2 * (a^-3 - b^-3)
}

eq2 <- function(a,b,k4){
  "See eq. 3.12 of Wu 2005"
  k4 - 6 * (a^-4 + b^-4)
}


grad.descent <- function(a,b,k1,k2,k3,k4,epsilon=.0001,delta=.0001,limit=Inf){
  "Perform a primitive psuedo-gradient descent optimization to solve the system eq1, eq2 for a and b, given k3,k4"
  print(c("entering grad descent w/ ",a,b,k1,k2,k3,k4))
  err <- function(a,b){eq1(a,b,k3)^2 + eq2(a,b,k4)^2}
  gen <- 0
  while(err(a,b) > epsilon && gen < limit) {
    deltas <- lapply(list(c(1,1),c(1,-1),c(-1,1),c(-1,-1)),
                     function(x) x * delta + c(a,b))
    winner <- deltas[which.min(lapply(deltas,
                                      function(pair)err(pair[1],pair[2])))][[1]]
    a <- winner[1]
    b <- winner[2]
    if(gen %% 10000 == 0){
#        mu <- k1 - 1/a + 1/b
#        sigma.squared <- k2 - (1/a^2 + 1/b^2)
#      print(paste(mu,sigma.squared,a,b,err(a,b)))}
        print(paste(a,b,err(a,b)))}
    gen <- gen + 1
  }
  c(a,b)
}

mme.estimate.nl <- function(xs,epsilon=.00001,delta=.001){
  "Given a sample xs from an NL distribution, estimate its parameters via method of moments, described in section 3.1 of Wu 2005"
  ms <- moments(xs)
  ks <- cumulants(ms)
  k1 <- ks[1]
  k2 <- ks[2]
  k3 <- ks[3]
  k4 <- ks[4]
  ab <- grad.descent(1,1,k1,k2,k3,k4,epsilon,delta)
  a <- ab[1]
  b <- ab[2]  
  mu <- k1 - 1/a + 1/b
  sigma.squared <- k2 - (1/a^2 + 1/b^2)
  c(mu,sigma.squared,a,b)
}

likelihood <- function(xs,mu,sigma.squared,a,b){
  "Return the sum of the log densities of xs according to NL(mu,sigma.squared,a,b).
The higher the returned value, the greater the likelihood of xs~NL(...)"
  sum(log((function(x)(dnormlap(x,mu,sigma.squared,a,b)))(xs)))
}

square <- function(x)x^2

finitify <- function(xs){
  mu <- mean(xs)
  sapply(xs,function(x) ifelse(is.finite(x),x,mu))
}

evaluate.grouping <- function(xs,params.list,js){
  sum(square(sapply(seq(1:length(params.list)),
                    function(i)likelihood(xs[i==js],
                                          params.list[[i]][1],
                                          params.list[[i]][2],
                                          params.list[[i]][3],
                                          params.list[[i]][4]))))
}

em <- function(ys,num.comps=2){
  "Do E-M for parameter estimation of ys~NL(nu,tau,alpha,beta). Cf. Reed 2004, section 4.3"
  n <- length(ys)
  ws <- rexp(n)
  zs <- rnorm(n)
                                          alphas <- rexp(num.comps)
#  alphas <- c(3,3)
                                          betas <- rexp(num.comps)
 # betas <- c(4,4)
                                         mus <- replicate(num.comps,mean(ys))
  print(mus)
#  mus <- c(0,5)
                                         sigmas <- replicate(num.comps,sqrt(var(ys)))
#  sigmas <- c(1,1)
  logl <- logl.old <- -Inf
  nc <- ifelse(num.comps==1,2,num.comps)
  taus <- (replicate(n,prob.vector(nc))) #randomize intial responsibilities
  pis <- prob.vector(num.comps)
  print(paste("pis:", pis))
  qs <- ps <- denoms <- ws <- wvs <- w.2s <- zs <- zs.2 <- mat.or.vec(nc,n)
  do <- TRUE
  gen <- 0
  while (do || (logl > logl.old) | gen < 5) {
    print(paste("em generation ",gen))
    gen <- gen + 1
    do <- FALSE
                                        #E step    
    for(i in 1:num.comps){
      ys.adjusted <- ys #ys*n*taus[i,]/sum(taus[i,])
      qs[i,] <- betas[i]  * sigmas[i] + (ys.adjusted - mus[i])/sigmas[i]
      ps[i,] <- alphas[i] * sigmas[i] - (ys.adjusted - mus[i])/sigmas[i]
      denoms[i,] <- (R(qs[i,]) + R(ps[i,]))
      ws[i,] <- sigmas[i]*(qs[i,]*R(qs[i,]) - ps[i,]*R(ps[i,]))/denoms[i,]
      wvs[i,] <- sigmas[i]*(1-ps[i,]*R(ps[i,]))/denoms[i,]
      w.2s[i,] <- sigmas[i]^2*((1+qs[i,]^2)*R(qs[i,]) +
                               (1+ps[i,]^2)*R(ps[i,]) -
                               sigmas[i]*(alphas[i]+betas[i]))/denoms[i,]
      zs[i,] <- ys.adjusted - ws[i,]
      zs.2[i,] <- ys.adjusted^2 - 2 * ys.adjusted * ws[i,] + w.2s[i,]
  }
                                        #M step
    fs <- sapply(1:num.comps,function(i){
      function(x)pis[i] * dnormlap(x,mus[i],sigmas[i]^2,alphas[i],betas[i])
    })
    f <- function(x) sum(sapply(1:num.comps,function(i) fs[[i]](x)))
    for(i in 1:num.comps){
      print(mus)
      print(paste("about to assign taus with parameters: ",mus[i],
            sigmas[i],alphas[i],betas[i]))
      taus[i,] <- sapply(ys,function(y)fs[[i]](y)/f(y))
      pis[i] <- sum(taus[i,])/n
      weighted.ks <- cumulants(moments(ys,taus[i,]))
      weighted.mean <- weighted.ks[1]
      mus[i] <- weighted.mean - 1/alphas[i] + 1/betas[i]#solving for (17)#mean(zs[i,])
      sigmas[i] <- sqrt(abs((sum(taus[i,]*(ys-weighted.mean)^2)/sum(taus[i,])
                    - 1/alphas[i]^2
                    - 1/betas[i]^2)))
                                        #sqrt(mean(zs.2[i,]) - mean(zs[i,])^2)
#      sigmas[i] <- sqrt(mean((zs[i,] - mean(zs[i,]))^2))
      ab <- grad.descent(alphas[i],betas[i],weighted.ks[1],weighted.ks[2],
                         weighted.ks[3],weighted.ks[4])
      alphas[i] <- ab[1]
      betas[i] <- ab[2]      
      ## A <- mean(wvs[i,])
      ## B <- A - mean(ws[i,])
      ## alphas[i] <- 1/(A + sqrt(A*B))
      ## betas[i] <- 1/(B + sqrt(A*B))
    }    
    logl.old <- logl
    logl <- log.l.mixture(ys,mus,sigmas^2,alphas,betas,pis)
    print(paste("pis:",pis))
    print(logl)
  }
  c(mus,sigmas^2,alphas,betas,pis)
}

em2 <- function(ys,num.comps=2){
  "Do E-M for parameter estimation of ys~NL(nu,tau,alpha,beta). Cf. Reed 2004, section 4.3"
  n <- length(ys)
  ws <- rexp(n)
  zs <- rnorm(n)
  alphas <- rexp(num.comps)
#  alphas <- c(3,3)
  betas <- rexp(num.comps)
 # betas <- c(4,4)
  mus <- replicate(num.comps,mean(ys)) + rnorm(num.comps)
#  mus <- c(0,5)
  sigmas <- replicate(num.comps,sqrt(var(ys))) + rexp(num.comps)
#  sigmas <- c(1,1)
  logl <- logl.old <- -Inf
  nc <- ifelse(num.comps==1,2,num.comps)
  taus <- (replicate(n,prob.vector(nc))) #randomize intial responsibilities
  pis <- prob.vector(num.comps)
  print(paste("pis:", pis))
  do <- TRUE
  gen <- 0
  plot(density(ys))
  while (do || (logl > logl.old) | gen < 5) {
    print(paste("em generation ",gen))
    do <- FALSE
                                        #E step
    fs <- sapply(1:num.comps,function(i){
      function(x)pis[i] * dnormlap(x,mus[i],sigmas[i]^2,alphas[i],betas[i])
    })
    f <- function(x) sum(sapply(1:num.comps,function(i) fs[[i]](x)))
    for(i in 1:num.comps){
      curve(fs[[i]](x),add=TRUE,type='line',col=ifelse(i%%2,'red','blue'))
      print(paste("about to assign taus for component", i, "with parameters: ",mus[i],
            sigmas[i],alphas[i],betas[i]))
      taus[i,] <- sapply(ys,function(y)fs[[i]](y)/f(y))
      pis[i] <- sum(taus[i,])/n
    }
   
                                        #M step
      for(i in 1:num.comps){
        print("eming")
        ys.adjusted <- ys*n*taus[i,]/sum(taus[i,])
        ms <- moments(ys,taus[i,])
        expectation <- ms[1]
        variance <- ms[2] - ms[1] ^2
        if(gen %% 2 == 0){ #then use alpha, beta to estimate mu, sigma
        print(paste("updating mu,sigma for component ",i,"because gen: ",gen))
          mus[i] <- expectation - 1/alphas[1] + 1/betas[1]
          sigmas[i] <- sqrt(abs(variance - (1/alphas[1]^2 + 1/betas[1]^2)))
       }
        else{#use mu, sigma to estimate alpha, beta
          print(paste("updating alpha,beta for component",i,"because gen: ",gen))
          print(paste("using values mu: ",mus[i],"sigma",sigmas[i],"expectation",expectation,"variance",variance))
          ab <- recover.ab(mus[i],sigmas[i],expectation,variance)
          alphas[i] <- ab[1]
          betas[i] <- ab[2]
          print(paste("updated alpha, beta to :",alphas[i],betas[i]))
        }
    }
    
    logl.old <- logl
    logl <- log.l.mixture(ys,mus,sigmas^2,alphas,betas,pis)
    print(paste("pis:",pis))
    print(paste("logl: ", logl))
    gen <- gen + 1
  }
  c(mus,sigmas^2,alphas,betas,pis)
}

recover.ab <- function(mu,sigma,expectation,variance,old.alpha=1,old.beta=1){
  print(paste("entering recover.ab: ",mu,sigma,expectation,variance))
  err1 <- function(mu,expectation,a,b){
    mu - (expectation - 1/a + 1/b)
  }
  err2 <- function(sigma,variance,a,b){
    print(paste("sigma: ",sigma))
   sigma^2  - (variance - (1/a^2 + 1/b^2))
  }
  err.total <- function(ab){
    a <- ab[1]
    b <- ab[2]
    err1(mu,expectation,a,b)^2 + err2(sigma,variance,a,b)^2
  }
  nlm(err.total,c(old.alpha,old.beta))$estimate
}

recover.ab.prime <- function(k3,k4,old.alpha=1,old.beta=1){
  err1 <- function(a,b){
    k3 - 2*(1/a^3 - 1/b^3)
  }
  err2 <- function(a,b){
    k4 - 6*(1/a^4 + 1/b^4)
  }

  err.total <- function(ab){
    a <- ab[1]
    b <- ab[2]
    err1(a,b)^2 + err2(a,b)^2
  }
  nlm(err.total,c(old.alpha,old.beta))$estimate
}

#xs <- replicate(1000000,rnormlap(0,2,3,4))
#ks <- cumulants(moments(xs))
#k3 <- ks[3]; k4 <- ks[4]

f1 <- function(a,b)   k3 - 2 * (a^-3 - b^-3)
f2 <- function(a,b)     k4 - 6 * (a^-4 + b^-4)
f <- function(p){
  a <- p[1]
  b <- p[2]
  f1(a,b)^2 + f2(a,b)^2
}

dwl <- function(x,mu,a,b){
  ifelse(x >= mu,(a*b/(a + b))*exp(a*(mu-x)),(a*b/(a + b))*exp(b*(x-mu)))
  
}

em3 <- function(ys,num.comps=2){
  "Do E-M for parameter estimation of ys~NL(nu,tau,alpha,beta). Cf. Reed 2004, section 4.3"
  n <- length(ys)
  ws <- rexp(n)
  zs <- rnorm(n)
  alphas <- rexp(num.comps)
#  alphas <- c(3,3)
  betas <- rexp(num.comps)
 # betas <- c(4,4)
  mus <- replicate(num.comps,mean(ys)) + rnorm(num.comps)
#  mus <- c(0,5)
  sigmas <- replicate(num.comps,sqrt(var(ys))) + rexp(num.comps)
#  sigmas <- c(1,1)
  logl <- logl.old <- -Inf
  nc <- ifelse(num.comps==1,2,num.comps)
  taus <- (replicate(n,prob.vector(nc))) #randomize intial responsibilities
  pis <- prob.vector(num.comps)
  print(paste("pis:", pis))
  do <- TRUE
  gen <- 0
  plot(density(ys))
  while (do || (logl > logl.old) | gen < 10) {
    print(paste("em generation ",gen))
    do <- FALSE
                                        #E step
    fs <- sapply(1:num.comps,function(i){
      function(x)pis[i] * dnormlap(x,mus[i],sigmas[i]^2,alphas[i],betas[i])
    })
    f <- function(x) sum(sapply(1:num.comps,function(i) fs[[i]](x)))
    for(i in 1:num.comps){
      print(paste("plotting component: ",i))
      curve(fs[[i]](x),add=TRUE,type='line',col=ifelse(i%%2,'red','blue'))
      print(paste("about to assign taus for component", i, "with parameters: ",mus[i],
            sigmas[i],alphas[i],betas[i]))
      taus[i,] <- sapply(ys,function(y)fs[[i]](y)/f(y))
      pis[i] <- sum(taus[i,])/n
    }
   
                                        #M step
      for(i in 1:num.comps){
        print("eming")
        ys.adjusted <- ys*n*taus[i,]/sum(taus[i,])
        ms <- moments(ys,taus[i,])
        expectation <- ms[1]
        variance <- ms[2] - ms[1] ^2

        if(gen %% 3 == 0){ 
        print(paste("updating mu,sigma for component ",i,"because gen: ",gen))
        mus[i] <- expectation - 1/alphas[1] + 1/betas[1]
        sigmas[i] <- sqrt(abs(variance - (1/alphas[1]^2 + 1/betas[1]^2)))
        }

        if(gen %% 3 > 1){ 
        print(paste("updating alpha,beta for component ",i,"because gen: ",gen))
        ab <- recover.ab(mus[i],sigmas[i],expectation,variance)
        if (gen %% 3 == 2)
          alphas[i] <- ab[1]
        else
          betas[i] <- ab[2]
        }

        
        ## if(gen %% 4 == 2){ 
        ##   print(paste("updating alpha for component ",i,"because gen: ",gen))

        ## }

        ##         if(gen %% 4 == 3){ 
        ## print(paste("updating beta for component ",i,"because gen: ",gen))

        ## }

    }
    
    logl.old <- logl
    logl <- log.l.mixture(ys,mus,sigmas^2,alphas,betas,pis)
    print(paste("pis:",pis))
    print(paste("logl: ", logl))
    gen <- gen + 1
  }
  c(mus,sigmas^2,alphas,betas,pis)
}

recover.params <- function(xs,params=c(0,1,1,1)){
  ll <- function(params){
  mu <- params[1]
  sigma.2 <- params[2]
  alpha <- params[3]
  beta <- params[4]
  -log.l(xs,mu,sigma.squared,alpha,beta) #to be minimized
}
  nlm(ll,params)
}

concentrated.likelihood <- function(abt,ys){
  "Eq 31. Reed 2004"
  a <- abt[1]
  b <- abt[2]
  t <- abt[3]
  y.bar <- mean(ys)
  n <- length(ys)
  (n * (log(a) + log(b) - log(a + b)) +
   sum(sapply(1:n,function(i)dnorm((ys[i]-y.bar + 1/a - 1/b)/t))) +
   sum(sapply(1:n,function(i)log(R(a*t - (ys[i]-y.bar + 1/a - 1/b)/t) +
                                  R(b*t - (ys[i]-y.bar + 1/a - 1/b)/t)))))
}

mle.params <- function(ys,mme.guess=c(0,1,1,1)){
  ks <- cumulants(moments(ys))
  print(ks)
  ab.mme <- recover.ab.prime(ks[3],ks[4])
  print(ab.mme)
  #approximate starting values for numeric MLE with MME
  mme.guess <- c(ab.mme,1)
  l.hat <- function(abt){
    -concentrated.likelihood(abt,ys)
  }
  # do MLE
  abt <- nlm(l.hat,mme.guess)
  print(abt)
  a <- abt[1]
  b <- abt[2]
  t <- abt[3]
  n <- mean(ys) - 1/a + 1/b
  c(n,t,a,b)
}
  
wilcox.solver <- function(ys,params=c(0,1,1,1)){
  n <- length(ys)
  ns <- rnorm(n)
  as <- rexp(n)
  bs <- rexp(n)
  f.wilcox <- function(params){
   mu <- params[1]
   sigma.2 <- sqrt(params[2])
   alpha <- params[3]
   beta <- params[4]
   1 - wilcox.test(xs,mu+sqrt(sigma.2)*ns+as/alpha-bs/beta)$p.value
 }
  nlm(f.wilcox,params)
}

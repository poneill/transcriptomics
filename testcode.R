library(lattice)
mle <- function(xs,epsilon=.00001,delta=.00001,m1=1,m2=1,verbose=FALSE){
  c1 <- mean(log(1+exp(xs)))
  c2 <- mean(log(1+exp(-xs)))
  lhs1 <- function(x,y)-(digamma(y) - digamma(x + y))
  lhs2 <- function(x,y)-(digamma(x) - digamma(x + y))
  err1 <- function(x,y) (lhs1(x,y) - c1)^2
  err2 <- function(x,y) (lhs2(x,y) - c2)^2
  err <- function(x,y) err1(x,y) + err2(x,y)
  gen <- 0
  old.err <- 0
  new.err <- err(m1,m2)
  while(err(m1,m2) > epsilon && (old.err != new.err)){
    deltas <- lapply(list(c(1,1),c(1,-1),c(-1,1),c(-1,-1)),
                     function(x) x * delta + c(m1,m2))
    winner <- deltas[which.min(lapply(deltas,
                                      function(pair)err(pair[1],pair[2])))][[1]]
    m1 <- winner[1]
    m2 <- winner[2]
    old.err <- new.err
    new.err <- err(m1,m2)
    if(verbose && gen %% 1000 == 0){
      print(paste(m1,m2,old.err,new.err,gen))
    }
    gen <- gen + 1
  }
  c(m1,m2)
}

grad.descent <- function(xs,epsilon=.001,delta=.00001,m1=1,m2=1){
  c1 <- mean(log(1+exp(xs)))
  c2 <- mean(log(1+exp(-xs)))
  lhs1 <- function(x,y)-(digamma(y) - digamma(x + y))
  lhs2 <- function(x,y)-(digamma(x) - digamma(x + y))
  x <- c(m1,m2)
  G <- function(x) matrix(c(lhs1(x[1],x[2]) - c1,lhs2(x[1],x[2]) - c2))
  F <- function(x) t(G(x))%*%G(x)
  Jg <- function(x) t(matrix(c(trigamma(x[1] + x[2]),
                                -trigamma(x[2]) + trigamma(x[1] + x[2]),
                                -trigamma(x[1]) + trigamma(x[1] + x[2]),
                                trigamma(x[1] + x[2])),ncol=2))
  grad.F <- function(x) t(Jg(x))%*%G(x)
  gen <- 0
  while(F(x) > epsilon){
    grad <- grad.F(x)
    hi <- delta
    lo <- 0
    while(hi!=lo){
#      print(c(lo,hi))
      d <- (hi+lo)/2
      x.prime <- x - d * grad
      if(F(x.prime) > F(x)){
        hi <- d
      }
      else{
        lo <- d
      }
    }
    x <- x.prime
    if(gen %% 1000 == 0)
      print(paste(gen,x))
  }
  x
    
}
meta.mle <- function(xs,epsilon=1,delta=1,m1=1,m2=2,factor=2){
  while(epsilon>0){
    ms <- mle(xs,epsilon,delta,m1,m2)
    print(ms)
    m1 <- ms[1]
    m2 <- ms[2]
    epsilon <- epsilon/factor
    delta <- delta/factor
  }
}
make.chain <- function(n,x,dtarget,rproposal,dproposal){
  chain <- array(x)
  while(length(chain) < n){
    y <- rproposal(x)
    a <- runif(1)
    p <- dtarget(x)
    p.prime <- dtarget(y)
    q <- dproposal(x,y)
    q.prime <- dproposal(y,x)
    if(a < (p.prime/p)*(q/q.prime)){
      x <- y
    }
    chain <- append(chain,x)
  }
  chain
}

f <- function(w,m1,m2) exp(w*m1)*(1+exp(w))^(-(m1+m2))/beta(m1,m2)
dprentice <- function(x,m1,m2){1/beta(m1,m2) * 1/((1+exp(-x))^m1*(1+exp(x))^m2)}
rprentice <- function(n,m1,m2){#simulate n draws from Prentice(m1,m2)
                               #via accept-reject sampling
  

}

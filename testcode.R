library(lattice)
mle <- function(c1,c2,epsilon=.00001,delta=.00001,m1=1,m2=1){
  lhs1 <- function(x,y){digamma(x) - digamma(x + y)}
  lhs2 <- function(x,y){digamma(y) - digamma(x + y)}
  err1 <- function(x,y) (lhs1(x,y) - c1)^2
  err2 <- function(x,y) (lhs2(x,y) - c2)^2
  err <- function(x,y) err1(x,y) + err2(x,y)
  gen <- 0
  while(err(m1,m2) > epsilon){
    deltas <- lapply(list(c(1,1),c(1,-1),c(-1,1),c(-1,-1)),
                     function(x) x * delta + c(m1,m2))
    winner <- deltas[which.min(lapply(deltas,
                                      function(pair)err(pair[1],pair[2])))][[1]]
    m1 <- winner[1]
    m2 <- winner[2]
    if(gen %% 1000 == 0) print(paste(m1,m2,err(m1,m2),gen))
    gen <- gen + 1
  }
  c(m1,m2)
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

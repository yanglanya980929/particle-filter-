library(ggplot2)
library(coda)
library(MASS)

# sample M-1 samples from 1:M
systematic_sampling <- function(weights, size) {
  # arguments: weights must be normalised weights
  # size: the number of samples
  M <- length(weights)
  cumulative_sum <- cumsum(weights)
  step <- 1 / size
  start <- runif(1, min = 0, max = step)
  positions <- seq(start, 1, by = step)
  new_samples <- integer(size)
  i <- 1
  for (j in 1:(size)){
    while (positions[j] > cumulative_sum[i]){
      i <- i + 1
    }
    new_samples[j] <- i
  }
  return(new_samples)
}

PF.sys <- function(simX0, simXt, loglike, simX0_params, simXt_params, loglike_params, ys, M, N, d){
  # create the objects
  X <- matrix(nrow=M, ncol=d, data=NA)
  X.filtering <- array(NA, dim=c(M, d, N)) # N slices
  ancestors <- matrix(NA, nrow = M, ncol = N) 
  ess <- numeric(N)
  
  ancestors[, 1] <- 1:M # initialising
  
  # sampling forwards
  for(i in 1:N){
    if (i == 1){
      X[,] <- simX0(simX0_params, M, d)
    }
    else{
      X[,] <- simXt(X.tilda, simXt_params, d)
    }
    # compute the normalised weights
    logW <- loglike(X, loglike_params, ys[i]) 
    logW.max <- max(logW)
    wstars <- exp(logW - logW.max)
    
    y.loglikely <- log(mean(wstars)) + logW.max # maximum marginal loglikelhood
    w.nor <- wstars/sum(wstars)
    ess[i] <- 1/sum(w.nor^2) # ESS
    
    index <- systematic_sampling(w.nor) # resampling
    if (i < N){
      ancestors[, i+1] <- index # store the ancestor of the samples at the next time point
    }
    X.tilda <- X[index,]
    X.filtering[, ,i] <- X.tilda # store the result sampled from the filtering distribution
  }
  return(list(X.filtering = X.filtering, X=X, weights=w.nor, ancestors=ancestors, y.loglikely=y.loglikely, ess=ess))
}

conPF.sys <- function(simX0, simXt, loglike, simX0_params, simXt_params, loglike_params, ys, M, N, d, path0){
  # create the objects
  X <- matrix(nrow=M-1, ncol=d, data=NA) # length M-1
  X.tilda <- matrix(nrow=M-1, ncol=d, data=NA) # length M-1
  X.predict <- array(NA, dim=c(M, d, N)) # N slices
  ancestors <- matrix(NA, nrow = M, ncol = N) # M*N
  
  X.predict[1, , ] <- path0 # assign path0 to x^1_{1:T}
  ancestors[, 1] <- 1:M # initialising
  ancestors[1, ] <- 1 # ancestors of the existing line 
  
  # sampling forwards
  for(i in 1:N){
    if (i == 1){
      X[,] <- simX0(simX0_params, M-1, d) # length M-1
    }
    else{
      X[,] <- simXt(X.tilda, simXt_params, M-1, d) # length M-1
    }
    X.predict[2:M, , i] <- X # length M-1
    X.whole <- matrix(X.predict[, , i], nrow = M, ncol = d) # make sure the slice is a matrix
    
    # compute the normalised weights
    logW <- loglike(X.whole, loglike_params, ys[i]) 
    logW.max <- max(logW)
    wstars <- exp(logW - logW.max) # length M
    
    y.loglikely <- log(mean(wstars)) + logW.max # maximum marginal loglikelhood
    w.nor <- wstars/sum(wstars)
    
    indexes <- systematic_sampling(w.nor, M-1) # resampling and get (M-1) indexes
    if (i < N){
      ancestors[2:M, i+1] <- indexes # store the ancestor of the samples at the next time point
    }
    X.tilda[,] <- X.whole[indexes,] # length M-1
  }
  return(list(X.predict=X.predict, weights=w.nor, ancestors=ancestors, y.loglikely=y.loglikely))
}

# conPF.test <- conPF.sys(simX0_2, simXt_2, loglike_2, simX0_params_2, simXt_params_2, loglike_params_2_low, ys.2.low, 10, 20, 1, xs.2)

########### only 'a' is unknown ############# 
ParticleGibbs <- function(simPara, simX0, simXt, loglike, simX0_params, simXt_params, loglike_params,ys,path0, M,num){
  d <- nrow(path0)
  N <- ncol(path0)
  Theta <- numeric(num+1) 
  a <- simXt_params[1] # initial value of 'a'
  b <- simXt_params[2]
  vx <- simXt_params[3]
  Theta[1] <- a
  Paths <- array(NA, dim=c(d, N, num+1))
  Paths[, ,1] <- path0
  path <- path0
  for (i in 1:num){
    a <- simPara(b, vx, path) # simulate 'a'
    Theta[i+1] <- a
    simXt_params[1] <- a
    # run a conditional SMC
    conPF <- conPF.sys(simX0, simXt, loglike, simX0_params, simXt_params, loglike_params, ys, M, N, d, path)
    weights <- conPF$weights
    k <- systematic_sampling(weights, 1)
    # get the new path back from k
    path <- matrix(NA, nrow = d, ncol = N)
    for(j in 1:N){
      x_k <- conPF$X.predict[k, ,N+1-j]
      path[ ,N+1-j] <- matrix(x_k, nrow = d, ncol = 1) # X.predict <- array(NA, dim=c(M, d, N))
      k <- conPF$ancestors[k,N+1-j]
    }
    Paths[, ,i+1] <- path
  }
  return(list(Theta=Theta, Paths=Paths))
}

############# Dynamic linear Model ########################
simX0_2 <- function(params, M, d){
  v0 <- params[1]
  x <- rnorm(M, mean = 0, sd = sqrt(v0))
  return(x)
}

simXt_2 <- function(x, params, M, d){
  a <- params[1]
  b <- params[2]
  vx <- params[3]
  xt <- rnorm(M, mean = a*x+b , sd = sqrt(vx))
  xt
}

loglike_2 <- function(xt, params, y){
  vy <- params[1]
  dnorm(y, mean = xt, sd = sqrt(vy), log=TRUE)
}

# simulate "a" exactly from the true posterior distribution
simPara.a <- function(b, vx, path){
  N <- length(path)
  m <- (sum(c((path - b)[-1],0)*path))/(sum(path**2)-path[N]**2)
  sigma2 <- vx/(sum(path**2)-path[N]**2)
  return(rnorm(1, mean = m, sd = sqrt(sigma2)))
}

d <- 1
M <- 100
N <- 100

# parameters
v0 <- 1 # starting from a point far away from zero
vx <- 1
vy.low <- 0.25
vy.high <- 4
d <- 1
a <- 0.8
b <- 0

# arguments 
simX0_params_2 <- c(v0)
simXt_params_2 <- c(a, b, vx)
loglike_params_2_low <- c(vy.low)
loglike_params_2_high <- c(vy.high)

# simulation 2
ts <- seq(from = 1, to = N, by = 1)
xs.2 <- rep(NA, N)
set.seed(123456)
xs.2[1] <- simX0_2(simX0_params_2, M, d)
for(i in 1:(N-1)){
  xs.2[i+1] <- simXt_2(xs.2[i], simXt_params_2, d)
}

# simulate two observed processes with low and high variance
set.seed(123456)
ys.2.low <- rnorm(N, mean = xs.2, sd = sqrt(vy.low))
set.seed(123456)
ys.2.high <- rnorm(N, mean = xs.2, sd = sqrt(vy.high))

plot(ts, ys.2.high, type="b",pch=1,lwd=2,xlab="t",ylab="x",col="orange") # high
lines(ts, xs.2, type="b",pch=4,lwd=2,xlab="t",ylab="x",col="blue") # hidden process
lines(ts, ys.2.low, type="b",pch=1,lwd=2,xlab="t",ylab="x",col="red") # low
abline(h=0,col="black",lty=2)
legend("bottomleft",c("Yt.low","Yt.high","Xt"),col=c("red","orange","blue"), lty = c(1, 1), lwd = c(2, 2))

########################### test of PG ########################
# path0 <- matrix(rep(1, N), nrow = 1, ncol = N)
path0 <- matrix(xs.2, nrow = 1, ncol = N)
num <- 2000
a0 <- 0.2
simXt_params_2 <- c(a0, b, vx) # suppose a is unknown and we start from a random value a0

PG.1 <- ParticleGibbs(simPara.a, simX0_2, simXt_2, loglike_2, simX0_params_2, simXt_params_2, loglike_params_2_low, ys.2.low, path0, M,num)
Theta.1 <- PG.1$Theta
Paths.1 <- PG.1$Paths

plot(Theta.1, type = "l", xlab = "Iteration", ylab = "Theta", main = "Trace Plot of Theta") # theta
plot(Paths.1[,N,], type = "l", xlab = "Iteration", ylab = "X_T", main = "Trace Plot of X_T") # X_T
plot(Paths.1[,1,], type = "l", xlab = "Iteration", ylab = "X_1", main = "Trace Plot of X_1") # X_1
a.hat <- mean(Theta.1[100:length(Theta.1)])

# 95% CI
se <-  sd(Theta.1[100:length(Theta.1)])/sqrt(length(Theta.1[100:length(Theta.1)]))
alpha <- 0.05
z_value <- qnorm(1 - alpha/2)  # Critical value for 95% CI from the normal distribution
CI.1 <- a.hat + c(-1,1)*se*z_value
CI.1

############### vx is unknown ############


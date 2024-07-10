library(ggplot2)
library(coda)

# We use PF is to get the approximation of two things: p(x_t|y_{1:t}, theta) and p(y_{1:T}|theta)
# so the outputs should be: X_T, w_T and \hat{p}(y_{1:T}|theta)
# In order to prepare for the Particle Gibbs, we also need the ancestor matrix

# Particle Filter using multinomial resampling
PF.multi <- function(simX0, simXt, loglike, simX0_params, simxt_params, loglike_params, ys, M, N, d){
  # create the objects
  X <- array(NA, dim = c(d, M, N))
  ancestors <- array(NA, dim =c(M, N))
  w.nor <- rep(1/M, M)
  
  #initialise the objects                   
  X[, , 1] <- simX0(simX0_params, M, d)
  ancestors[, 1] <- 1:M
  y.loglikely <- 0
  
  # update
  for(i in 1:(N-1)){
    index <- sample(1:M, size = M, replace = TRUE, prob = w.nor) # resample the index 1:M based on weights[, i]
    ancestors[, i+1] <- index # the ancestor of the (i+1)th generation of the particles
    x.tilda <- X[, , i][index] # get the resampled particles
    X[, , i+1] <- simXt(x.tilda, simxt_params, d) # sample from the transition kernel based on x.tilda
    
    logW <- loglike(X[, , i+1], loglike_params, ys[i]) # compute the log weights for X[, , i+1]
    logW.max <- max(logW)
    wstars <- exp(logW - logW.max)
    
    y.loglikely <- y.loglikely + log(mean(wstars)) + logW.max # maximum marginal loglikelhood
    w.nor <- wstars/sum(wstars) # normalised weight
  }
  return(list(X=X, weights=w.nor, ancestors=ancestors, y.loglikely=y.loglikely))
}

# Particle Filter using systematic resampling
systematic_sampling <- function(samples, weights) {
  # arguments: weights must be normalised weights
  M <- length(samples)
  cumulative_sum <- cumsum(weights)
  step <- 1 / M
  start <- runif(1, min = 0, max = step)
  positions <- seq(start, 1, by = step)
  new_samples <- numeric(M)
  i <- 1
  for (j in 1:M){
    while (positions[j] > cumulative_sum[i]){
      i <- i + 1
    }
    new_samples[j] <- samples[i]
  }
  return(new_samples)
}

PF.sys <- function(simX0, simXt, loglike, simX0_params, simxt_params, loglike_params, ys, M, N, d){
  # create the objects
  X <- array(NA, dim = c(d, M, N))
  ancestors <- array(NA, dim =c(M, N))
  w.nor <- rep(1/M, M)
  ess <- numeric(N)
  
  #initialise the objects                   
  X[, , 1] <- simX0(simX0_params, M, d)
  ancestors[, 1] <- 1:M
  y.loglikely <- 0
  ess[1] <- M
  
  # update
  for(i in 1:(N-1)){
    index <- systematic_sampling(1:M, w.nor) # resample the index 1:M based on w.nor
    ancestors[, i+1] <- index # the ancestor of the (i+1)th generation of the particles
    x.tilda <- X[, , i][index] # get the resampled particles
    X[, , i+1] <- simXt(x.tilda, simxt_params, d) # sample from the transition kernel based on x.tilda
    
    logW <- loglike(X[, , i+1], loglike_params, ys[i]) # compute the log weights for X[, , i+1]
    logW.max <- max(logW)
    wstars <- exp(logW - logW.max)
    
    y.loglikely <- y.loglikely + log(mean(wstars)) + logW.max # maximum marginal loglikelhood
    w.nor <- wstars/sum(wstars) # normalised weight
    
    ess[i+1] <- 1/sum(w.nor^2)
  }
  return(list(X=X, weights=w.nor, ancestors=ancestors, y.loglikely=y.loglikely, ess=ess))
}

# visualisation
quantiles <- function(x) {
  quantile(x, probs = c(0.05, 0.5, 0.95))
}

############### Example 1: stochastic-volatility model #########################
simX0_1 <- function(params, M, d){
  sigma <- params[1]
  phi <- params[2]
  x <- rnorm(M, mean = 0, sd = sqrt(sigma^2/(1-phi^2)))
  return(x)
}

loglike_1 <- function(xt, params, y){
  K <- params[1]
  dnorm(y, mean = 0, sd = K*sqrt(exp(xt)), log=TRUE)
}

simXt_1 <- function(x, params, d){
  phi <- params[1]
  sigma <- params[2]
  epsilon <- rnorm(M, mean = 0, sd = 1)
  xt <- phi*x + sigma*epsilon
  xt
}

# parameters
d <- 1
M <- 10
N <- 15
sigma <- 0.4
phi <- 0.95
K <- 0.6

### arguments ###################
# sigma, phi, K

simX0_params_1 <- c(sigma, phi)
simxt_params_1 <- c(phi, sigma)
loglike_params_1 <- c(K)

# simulation 1
ts <- seq(from = 1, to = N, by = 1)
xs <- rep(NA, N)
set.seed(1234567)
xs[1] <- simX0_1(simX0_params_1, 1)
for(i in 1:(N-1)){
  xs[i+1] <- simXt_1(xs[i], simxt_params_1)
}
# print(xs)
ys_1 <- rnorm(N, mean = 0, sd = K*sqrt(exp(xs)))
plot(ts, xs,type="b",pch=1,lwd=2,xlab="t",ylab="x",col="blue") # hidden process
plot(ts, ys_1,xlab="t",ylab="y",type="b",lwd=2,pch=4)
abline(h=0,col="black",lty=2)

######################## PF test1 ######################################
set.seed(1234567)
pf.1 <- PF.multi(simX0_1, simXt_1, loglike_1, simX0_params_1, simxt_params_1, loglike_params_1, ys_1, M, N, d)
X.1 <- pf.1$X
weights.1 <- pf.1$weights
ancestors.1 <- pf.1$ancestors
y.loglikely.1 <- pf.1$y.loglikely # -23.00247

# visualisation
qs <- matrix(nrow=3, ncol=N, data = NA)
for(i in 1:N){
  qs[,i] <- quantiles(X.1[1, ,i])
}

Ts=(1:N)
qlo=min(qs); qhi=max(qs)
plot(Ts,qs[2,],type="l",lwd=2,col="red",xlab="time",ylab="position",
     ylim=c(qlo,qhi))
lines(Ts,qs[1,],lwd=2,lty=2,col="red")
lines(Ts,qs[3,],lwd=2,lty=2,col="red")
points(Ts,ys_1,pch=4,lwd=2)
lines(Ts,xs,col="blue",lty=1,lwd=3)
title(main = "stochastic-volatility model")
legend("bottomright",c("q50","q5","q95","xt"),col=c("red","red","red","blue"),
       lty=c(1,2,2,1),lwd=c(2,2,2,3))

############################ Example 2 Dynamic linear model ###################################
simX0_2 <- function(params, M, d){
  v0 <- params[1]
  x <- rnorm(M, mean = 0, sd = sqrt(v0))
  return(x)
}

simXt_2 <- function(x, params, d){
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

d <- 1
M <- 10
N <- 20

# parameters
v0 <- 50 # starting from a point far away from zero
vx <- 0.2
vy.low <- 0.05
vy.high <- 5
d <- 1
a <- 0.3
b <- 0.1

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
print(xs.2[1])
for(i in 1:(N-1)){
  xs.2[i+1] <- simXt_2(xs.2[i], simXt_params_2, d)
}

# simulate two observed processes with low and high variance
set.seed(123456)
ys.2.low <- rnorm(N, mean = xs.2, sd = sqrt(vy.low))
ys.2.high <- rnorm(N, mean = xs.2, sd = sqrt(vy.high))

plot(ts, ys.2.high, type="b",pch=1,lwd=2,xlab="t",ylab="x",col="orange") # high
lines(ts, xs.2, type="b",pch=4,lwd=2,xlab="t",ylab="x",col="blue") # hidden process
lines(ts, ys.2.low, type="b",pch=1,lwd=2,xlab="t",ylab="x",col="red") # low
abline(h=0,col="black",lty=2)
legend("topright",c("Yt.low","Yt.high","Xt"),col=c("red","orange","blue"), lty = c(1, 1), lwd = c(2, 2))

######################## PF test2 ######################################
# # PF.multi
# set.seed(123456)
# pf.2 <- PF.multi(simX0_2, simXt_2, loglike_2, simX0_params_2, simXt_params_2, loglike_params_2, ys.2, M, N, d)
# X.2 <- pf.2$X
# weights.2 <- pf.2$weights
# ancestors.2 <- pf.2$ancestors
# y.loglikely.2 <- pf.2$y.loglikely 
# 
# qs.2 <- matrix(nrow=3, ncol=N, data = NA)
# for(i in 1:N){
#   qs.2[,i] <- quantiles(X.2[1, ,i])
# }
# 
# Ts=(1:N)
# qlo=min(qs.2); qhi=max(qs.2)
# plot(Ts,qs.2[2,],type="l",lwd=2,col="red",xlab="time",ylab="position",
#      ylim=c(qlo,qhi))
# lines(Ts,qs.2[1,],lwd=2,lty=2,col="red")
# lines(Ts,qs.2[3,],lwd=2,lty=2,col="red")
# points(Ts,ys.2,pch=4,lwd=2)
# lines(Ts,xs.2,col="blue",lty=1,lwd=3)
# title(main = "dynamic linear model(multi)")
# legend("bottomright",c("q50","q5","q95","xt"),col=c("red","red","red","blue"),
#        lty=c(1,2,2,1),lwd=c(2,2,2,3))

########### PF.sys ###############
set.seed(1234567)
pf.2.low <- PF.sys(simX0_2, simXt_2, loglike_2, simX0_params_2, simXt_params_2, loglike_params_2_low, ys.2.low, M, N, d)
pf.2.high <- PF.sys(simX0_2, simXt_2, loglike_2, simX0_params_2, simXt_params_2, loglike_params_2_high, ys.2.high, M, N, d)

# low vy
X.2.low <- pf.2.low$X
weights.2.low <- pf.2.low$weights
ancestors.2.low <- pf.2.low$ancestors
y.loglikely.2.low <- pf.2.low$y.loglikely 
ess.2.low <- pf.2.low$ess

qs.2.low <- matrix(nrow=3, ncol=N, data = NA)
for(i in 1:N){
  qs.2.low[,i] <- quantiles(X.2.low[1, ,i])
}

# high vy
X.2.high <- pf.2.high$X
weights.2.high <- pf.2.high$weights
ancestors.2.high <- pf.2.high$ancestors
y.loglikely.2.high <- pf.2.high$y.loglikely 
ess.2.high <- pf.2.high$ess

qs.2.high <- matrix(nrow=3, ncol=N, data = NA)
for(i in 1:N){
  qs.2.high[,i] <- quantiles(X.2.high[1, ,i])
}

######### visualisation #########
# ESS
Ts=(1:N)
plot(Ts, ess.2.low, type="b",lwd=2,pch=1,col="red",xlab="time",ylab="ess")
lines(Ts, ess.2.high, type="b",lwd=2,pch=4,lty=2,col="blue")
title(main = "ESS for dynamic linear model")
legend("topright",c("low","high"),col=c("red","blue"),
       lty=c(1,2,2,1),lwd=c(2,2,2,3))

#low
Ts=(1:N)
qlo=min(qs.2.low); qhi=max(qs.2.low)
plot(Ts,qs.2.low[2,],type="l",lwd=2,col="red",xlab="time",ylab="position",
     ylim=c(qlo,qhi))
lines(Ts,qs.2.low[1,],lwd=2,lty=2,col="red")
lines(Ts,qs.2.low[3,],lwd=2,lty=2,col="red")
points(Ts,ys.2.low,pch=4,lwd=2)
lines(Ts,xs.2,col="blue",lty=1,lwd=3)
title(main = "dynamic linear model(low)")
legend("bottomright",c("q50","q5","q95","xt"),col=c("red","red","red","blue"),
       lty=c(1,2,2,1),lwd=c(2,2,2,3))

# high
qlo=min(qs.2.high); qhi=max(qs.2.high)
plot(Ts,qs.2.high[2,],type="l",lwd=2,col="red",xlab="time",ylab="position",
     ylim=c(qlo,qhi))
lines(Ts,qs.2.high[1,],lwd=2,lty=2,col="red")
lines(Ts,qs.2.high[3,],lwd=2,lty=2,col="red")
points(Ts,ys.2.high,pch=4,lwd=2)
lines(Ts,xs.2,col="blue",lty=1,lwd=3)
title(main = "dynamic linear model(high)")
legend("bottomright",c("q50","q5","q95","xt"),col=c("red","red","red","blue"),
       lty=c(1,2,2,1),lwd=c(2,2,2,3))

# the function of obtaining the num of ancestors 
aces_num <- function(aces_matrix){
  ancestors_num <- rep(NA, N)
  ancestors_num[N] <- M
  aces <- aces_matrix[,N]
  for(i in 1:(N-1)){
    index <- unique(aces) # make the index unique
    ancestors_num[N-i] <- length(index) # the number of ancestors
    aces <- aces_matrix[index,N-i] # get the index of the previous ancestors
  }
  return(ancestors_num)
}

Ts=(1:N)
aces.low <- aces_num(ancestors.2.low)
aces.high <- aces_num(ancestors.2.high)
plot(Ts, aces.low, type="b",lwd=2,pch=1,col="red",xlab="time",ylab="num of ancestors")
lines(Ts, aces.high, type="b",lwd=2,pch=4,lty=2,col="blue")
title(main = "num of ancestors for dynamic linear model")
legend("topright",c("low","high"),col=c("red","blue"),
       lty=c(1,2,2,1),lwd=c(2,2,2,3))

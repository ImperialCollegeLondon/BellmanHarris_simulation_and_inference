set.seed(123)
theta1 = 3
theta2 = 1
rbase = function(n) rgamma(n,theta1,theta2)
dbase = function(n) dgamma(n,theta1,theta2)
pbase = function(n) pgamma(n,theta1,theta2)

Rt=function(x){ 
  R0 + sin(0.25*x)
}
dx=1
max_time <- 50
R0 <- 1.3
xx=seq(1,max_time,dx)
nn=length(xx)

g = dbase(xx)
g=g/sum(g)

A = matrix(0,nrow=nn,ncol=nn)
C = matrix(NA,nrow=nn,ncol=nn)
G = matrix(NA,nrow=nn,ncol=nn)
for(i in 1:nn){
  A[,i] = Rt(pmax(0,xx-xx[i]))
  C[i,] = g
  G[,i] = (1-cumsum(g))[i]
  
}
A[upper.tri(A)] <- 0


A2 = matrix(0,nrow=nn,ncol=nn)
R=Rt(xx)
for(i in 1:nn){
  for(j in 1:i)
    A2[,i] = R[i-j+1]
}


B1 = matrix(0,nrow=nn,ncol=nn)
B1[,1]=G[,1]
B1[,2]=G[,2]+C[,(2-1):1]*A[,1:(2-1)]*B1[,1:(2-1)]

#all further steps
library(tictoc)
tic()
pb <- txtProgressBar(min = 0, max = nn, style = 3)
for(i in 3:nn){
  convolution1=rowSums(C[,(i-1):1]*A[,1:(i-1)]*B1[,1:(i-1)],na.rm=T)

  B1[,i] = G[,i] + convolution1  
  setTxtProgressBar(pb, i)
}
toc()
plot(xx,diag(B1),type='l')


x=xx
n=length(x)
response1=diag(B1) + rnorm(n,0,0.4*diag(B1))


library(splines)
b=10
Basis =bs(x,df=b,degree=3)
library(cmdstanr)
library(posterior)
library(bayesplot)
file <- file.path('branching_process.stan')

mod <- cmdstan_model(file,cpp_options=list("STAN_NO_RANGE_CHECKS"=TRUE))

Truncation=20
stan_data = list(N=n,prev=ceiling(response1),g=g,G=G,C=C,x=x,Basis=Basis,B=b,T=Truncation)

fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  refresh = 50,
  iter_sampling =200,
  iter_warmup=200
)
stanfit <- rstan::read_stan_csv(fit$output_files())
library(rstan)
par(mfrow=c(1,2))
e=extract(stanfit,'output1')
plot(x,colMeans(e$output1),type='l',lwd=2)
lines(x,stan_data$prev,col='red',type='p',pch=16)
lines(x,diag(B1),col='blue',type='l',pch=16)

e=extract(stanfit,'output2')
plot(x,colMeans(e$output2),type='l')
lines(x,Rt(x),col='red')

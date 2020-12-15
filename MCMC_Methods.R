## The code in this module shows some of the mainstream MCMC methods used in research.
###SIMULATION OF RANDOM VARIABLES:LAB 2
n=1000
x=runif(n,0,1) ##simulate 100 realizations from U(0,1)
x
hist(x)
###implement inverse method to simulate exponential 
#1. simulate from uniform
#2. x=-(1/lambda)*log(U)
lam=2
n=1000; U=runif(n)
E=-(1/lam)*log(U)
hist(E)
##simulate Bernoulli(p) from a uniform dist.
p=0.5
n=1000; U=runif(n)
hist(U)
X=ifelse(U<(1-p),0,1)
X
##simulate a gamma distribution from a uniform distribution
#x=c(1,2,3); sum(x)
alpha=2;beta=1;n=1000
U=matrix(runif(n*alpha),alpha,n)
E=-log(U) ##exp(1)
G=beta*colSums(E)
hist(G)
###chi-sq and beta dist simulation for HW
###make a custom function to simulate gamma's 
gamma_sim=function(nsim=1000,alpha=3,beta=2){
  U=matrix(runif(nsim*alpha),alpha,nsim)
  E=-log(U) ##exp(1)
  G=beta*colSums(E) 
  return(G)}
x=gamma_sim(10000,4,0.5)
hist(x)
x=gamma_sim(alpha=4,beta=0.5)
x=gamma_sim()
hist(x)
#############################################################################################
## SIMULATION OF A NORMAL DISTRIBUTION USING A BOX-MUELLER METHOD, AND THE ACCEPT REJECT METHOD
##Normal distribution using Box-Mueller method
##u1, u2~U(0,1); x1=sqrt(-2log u1)cos (2piu2); x2=..
nsim=5000; ns=ceiling(nsim/2)
u1=runif(ns,0,1);u2=runif(ns,0,1)
x1=(sqrt(-2*log(u1)))*cos(2*pi*u2)
x2=(sqrt(-2*log(u1)))*sin(2*pi*u2)
x=c(x1,x2) ##nsim realizations of a N(0,1) distribution
length(x)
hist(x)
normal_sim=function(nsim){
  ns=ceiling(nsim/2);  u1=runif(ns,0,1);u2=runif(ns,0,1)
  x1=(sqrt(-2*log(u1)))*cos(2*pi*u2)
  x2=(sqrt(-2*log(u1)))*sin(2*pi*u2)
  x=c(x1,x2)[1:nsim]; return(x)}
normal_sim(15)
###once we have a N(0,1) we can simulate any univariate 
#or multivariate normal distribution
##to simulate N(mu,si^2)
#x~N(0,1)---> z=mu+si*x ~ N(mu,si^2)
norm_sim=function(nsim,mu=0,var=1){ x=normal_sim(nsim); 
z=mu+(sqrt(var))*x;  return(z)}
z=norm_sim(1000,mu=3,var=2)
hist(z)
###multivariate normal distribution
##suppose we want to simulate N_p(mu,SI)
##to generate z=(z_1,...,z_p)~N(mu,Si)
#x=(x_1,...,x_p)~N(0,1); z=mu+Ga*x
nsim=100; mu=c(0,0); sig=matrix(c(1,0.5,0.5,1),2,2)
sig
multi_normal=function(nsim,mu,sig){
  p=length(mu)
  x=matrix(normal_sim(nsim*p),nsim,p)
  Gamma=chol(sig)
  z=mu+t(Gamma)%*%t(x)
  return(t(z))}
x=multi_normal(1000,mu,sig)
cor(x)
###Accept reject algorithm 
##N(0,1)--> target distribution
##DE(al) --> instrumental distribution
#1. x~g; u~u(0,1)
#2. y=x, u <= f(x)/Mg(x)
#3. return to 1 otherwise
#function to simulate from the instrument
DE_sim=function(nsim,al=1){
  u1=runif(nsim,0,1);u2=runif(nsim,0,1)
  e1=-(1/al)*log(u1); e2=-(1/al)*log(u2)
  x=e1-e2;return(x)}
DE_sim(10,2)
bound=function(M,x,al){
  num=(1/sqrt(2*pi))*exp(-(x^2)/2)
  den=M*(al/2)*exp(-al*abs(x))
  ratio=num/den
  return(ratio)}
al=1
M=(sqrt(2/pi))*(1/al)*exp(al^2/2)
M
bound(M,2,al)
bound(M,-1,al)
nsim=5000
ysim=c()
while(length(ysim)<nsim){
  x=DE_sim(1,al); u=runif(1)
  y=ifelse(u<=bound(M,x,al),x,NA)
  if(is.na(y)==F){ysim=c(ysim,y)}}
length(ysim)
hist(ysim)
#############################################################################################
##MARSAGLIA'S POLAR METHOD, Envelope Accept Reject(EAR), MC approximations.
##normal simulation using marsaglia's polar method
nsim=5000;x1=c();x2=c();its=0
for(j in 1:10^8){
  u1=runif(1,-1,1); u2=runif(1,-1,1)
  s=u1^2+u2^2
  if(s<1){x1[j]=sqrt(-2*log(s))*(u1/sqrt(s))
  x2[j]=sqrt(-2*log(s))*(u2/sqrt(s))
  its=its+1}
  if(its>=ceiling(nsim/2)){break}}
x1=x1[is.na(x1)==F]; x2=x2[is.na(x2)==F]
x=c(x1,x2)
length(x); hist(x)
##make a custom function of the above code
normal_sim=function(nsim){
  x1=c();x2=c();its=0
  for(j in 1:10^8){
    u1=runif(1,-1,1); u2=runif(1,-1,1)
    s=u1^2+u2^2
    if(s<1){x1[j]=sqrt(-2*log(s))*(u1/sqrt(s))
    x2[j]=sqrt(-2*log(s))*(u2/sqrt(s))
    its=its+1}
    if(its>=ceiling(nsim/2)){break}}
  x1=x1[is.na(x1)==F]; x2=x2[is.na(x2)==F]
  x=c(x1,x2)[1:nsim]
  return(x)}
normal_sim(50)
###Envelope accept reject###
##1/2pi(1-1/x^2)<=fx <= M(DE)
###Monte carlo approximations
#integrate x from 0 to 2 (true value : 2)
nsim=5000; mc=c(); v=c()
upper=c();lower=c()
u=runif(nsim,0,2)
for(j in 1:nsim){
  mc[j]=2*mean(u[1:j])
  v[j]=4*(1/j)*var(u[1:j])
  upper[j]=mc[j]+1.96*sqrt(v[j])
  lower[j]=mc[j]-1.96*sqrt(v[j])}
plot(mc,type="l")
abline(h=2)
#install.packages("ggplot2")
library(ggplot2)
values=c(mc,upper,lower); type=c(rep("mc",nsim),rep("upper",nsim),                                 rep("lower",nsim))
itr=seq(1,nsim)
data=data.frame(values=values,type=type,itr=itr)
ggplot(data,aes(itr,values))+geom_line(aes(color=type))
##approximate P(X<-1); X~N(0,1)
nsim=5000
x=normal_sim(nsim)
hist(x)
y=ifelse((x<(-1)),1,0)
###<- as an =
mc=c(); v=c()
upper=c();lower=c()
for(j in 1:nsim){
  mc[j]=mean(y[1:j])
  v[j]=(1/j)*var(y[1:j])
  upper[j]=mc[j]+1.96*sqrt(v[j])
  lower[j]=mc[j]-1.96*sqrt(v[j])}
plot(mc,type="l")
abline(h=0.16)
mean(x<0)
###############################################################################################
##MC APPROXIMATIONS, IMPORTANCE SAMPLING: LAB 5
###MC approximation
##Q. Approximate P(0<X<3), X~cauchy(0,1) using MC approx
##X1~N(0,1); X2~N(0,1); X1/X2~Cauchy(0,1)
## sum_{j=1}^{nsim}1[0<X_j<3]/nsim
nsim=1e4
x1=normal_sim(nsim);x2=normal_sim(nsim)
x=x1/x2;hist(x)
indicator=function(x){
  y=ifelse((x>0 & x<3),1,0)
  return(y)}
#indicator(4);indicator(2)
mc=c(); v=c();upper=c();lower=c()
for(j in 1:nsim){
  mc[j]=mean(indicator(x[1:j]))  
  v[j]=(j^{-1})*var(indicator(x[1:j]))
  upper[j]=mc[j]+1.96*sqrt(v[j])
  lower[j]=mc[j]-1.96*sqrt(v[j])}
plot(mc,type="l")
#plot(mc)
##Q2. Do the same question using the importance
#sampling estimate using U(0,3) as an instrument.
## h*f/g
cdensity=function(x){y=1/(pi*(1+(x^2)));return(y)}
cdensity(1)
nsim=1e4
x=runif(nsim,0,3)
mc=mean(3*cdensity(x))
mc=c(); v=c();upper=c(); lower=c()
for(j in 1:nsim){mc[j]=mean(3*cdensity(x[1:j]))
v[j]=(j^{-1})*var(3*cdensity(x[1:j]))
upper[j]=mc[j]+1.96*sqrt(v[j])
lower[j]=mc[j]-1.96*sqrt(v[j])}
plot(mc,type="l")
#################################################################################################
### NEWTON RHAPSON UPDATE, EM ALGORITHM
##Newton Rhapson update
##approximate the value of sqrt(2)
##sqrt(2) is the solution of x^2-2=0
##maximization of x^3-2x
g=function(x){x^2-2}
gprime=function(x){2*x}
#options(digits=10)
tol=1e-10;xold=1;err=1
maxits=3000;its=1
while(err>tol & its<maxits){
  xnew=xold-g(xold)/gprime(xold)
  err=abs(xnew-xold);xold=xnew;its=its+1}
##using the for loop instead of while loop
tol=1e-10;xold=1;err=1
maxits=3000;its=1
for(i in 1:maxits){
  xnew=xold-g(xold)/gprime(xold)
  err=abs(xnew-xold);xold=xnew;its=its+1
  if(err<tol){break}}
###NR update for GLM: HW 
###EM algorithm for mixture of normals
##simulate data
n=500; pr=0.3; mu1=-1;mu2=3;si1=0.5;si2=1
z1=rbinom(n,1,pr);z2=1-z1
x1=rnorm(n,mu1,si1);x2=rnorm(n,mu2,si2)
x=(z1*x1)+(z2*x2)
hist(x,breaks=25)
x
###from here we assume that we only observe x
##implement the em algorithm to estimate
##pr, mu1,mu2,si1,si2
##implementing the update from notes

f=function(x,mu,si){
  fx=(1/(si*sqrt(2*pi)))*exp(-((x-mu)^2)/((2*(si^2))))
  return(fx)}
f(1,0,1)
f(0,0,1)
f(-1,0,1)

oldpr=0.5;oldmu1=-1;oldmu2=1;oldsi1sq=1;oldsi2sq=1
tol=1e-5;err=1;maxits=3000;its=1
n=length(x)
while(err>tol & its<maxits){
  ##update for probability
  num=oldpr*f(x,oldmu1,sqrt(oldsi1sq))
  den=oldpr*f(x,oldmu1,sqrt(oldsi1sq))+
    (1-oldpr)*f(x,oldmu2,sqrt(oldsi2sq))  
  oldz1=num/den
  oldz2=1-oldz1
  newpr=sum(oldz1)/n 
  ##update for means
  newmu1=sum(x*oldz1)/sum(oldz1)
  newmu2=sum(x*oldz2)/sum(oldz2)
  ##update for variances
  newsi1sq=sum(oldz1*((x-newmu1)^2))/sum(oldz1)
  newsi2sq=sum(oldz2*((x-newmu2)^2))/sum(oldz2)
  err=max(abs(c(newpr-oldpr,newmu1-oldmu1,
                newmu2-oldmu2,newsi1sq-oldsi1sq,
                newsi2sq-oldsi2sq)))
  its=its+1
  oldpr=newpr; oldmu1=newmu1;
  oldmu2=newmu2;oldsi1sq=newsi1sq;
  oldsi2sq=newsi2sq
  
}

newpr
newmu1
newmu2
newsi1sq
newsi2sq
################################################################################################
## SIMULATED ANNEALING: LAB
###Example from class: cauchy MLE
x=c(-4.8,-2.8,-1.35,-0.02,0.70,0.98,2.92,5.50)
##assuming these are realizations from a C(alpha,0.1)
###MLE of alpha using SA
cauchy.l=function(alpha){
  y=-sum(log(0.01+(x-alpha)^2))
  return(y)}
cauchy.l(1)
al=seq(-10,10,length.out=500)
y=sapply(as.list(al),cauchy.l)
plot(al,y,type="l")
#implement SA to find MLE of alpha
Tt=function(t){y=1/((1+t)^(1/5));return(y)}
Tt=function(t){y=1/(log(1+t))^{1/2};return(y)}
Tt(1)
Tt(2)
Tt(3)
N=25000
alpha.old=-5;r=0.5;updates=alpha.old
for(i in 1:N){
  alpha.can=runif(1,alpha.old-r,alpha.old+r)
  f.can=cauchy.l(alpha.can); f.old=cauchy.l(alpha.old)
  if(f.can>f.old){alpha.new=alpha.can}else{
    rho=exp((f.can-f.old)/Tt(i)); b=rbinom(1,1,rho)
    if(b==1){alpha.new=alpha.can}
    if(b==0){alpha.new=alpha.old}}
  updates=c(updates,alpha.new)
  alpha.old=alpha.new}
plot(updates,type="l")
alpha.new
f.can
f.old
alpha.can
##assuming these are realizations from a C(alpha,beta)
###MLE of alpha using SA: beta>0
cauchy.l=function(alpha,beta=0.5){
  y=-sum(log(beta^2+(x-alpha)^2))
  return(y)}
cauchy.l(1,1)
al=seq(-10,10,length.out=500)
y=sapply(as.list(al),cauchy.l)
plot(al,y,type="l")
#implement SA to find MLE of alpha
Tt=function(t){y=1/((1+t)^(1/5));return(y)}
Tt=function(t){y=1/(log(1+t))^{1/2};return(y)}
Tt(1)
Tt(2)
Tt(3)
N=25000
alpha.old=-5;r=0.5;updates=alpha.old
for(i in 1:N){
  alpha.can=runif(1,alpha.old-r,alpha.old+r)
  f.can=cauchy.l(alpha.can); f.old=cauchy.l(alpha.old)
  if(f.can>f.old){alpha.new=alpha.can}else{
    rho=exp((f.can-f.old)/Tt(i)); b=rbinom(1,1,rho)
    if(b==1){alpha.new=alpha.can}
    if(b==0){alpha.new=alpha.old}}
  updates=c(updates,alpha.new)
  alpha.old=alpha.new}
plot(updates,type="l")
alpha.new

f.can
f.old
alpha.can

######## METROPOLIS HASTINGS ALGORITHM
##Suppose we want to simulate a markov chain with a stationary dist
##f(x) propto x^3-x, 1<x<5.
f=function(x){
  if(1<x & x<5){y=x^3-x}else{y=0}
  return(y)}
f(1);f(2);f(10)
N=10000
##implement : independent MH:
#intialize: 
x=c();x[1]=1
##simulate a candidate
##proposal density: U(0,5)
for(j in 1:N){
  y=runif(1,0,5)
  ##decide which value to use as the next step in the chain
  ##rho=min{(f(y)/f(x))*(q(x)/q(y))}
  rho=min((f(y)/f(x[j])),1)
  if(rho==1){x[j+1]=y}
  if(rho!=1){coin=rbinom(1,1,rho)
  if(coin==1){x[j+1]=y}
  if(coin==0){x[j+1]=x[j]}}}
x[1:10];hist(x,breaks=25)
###NOTE: to get rid of the initializers effect, 
#it is common to have a burn in period. 
#(start the chain after a certain number of realizations have passed.)
x.burn=x[5001:N]
hist(x.burn,breaks=25)
###using random walk MH
##target: Normal distribution (0,1)
f=function(x){y=exp(-x^2/2);return(y)}
f(1);f(2);f(10)
N=10000
##implement : random walk MH
#intialize: 
x=c();x[1]=1
##simulate a candidate
##proposal density: U(0,5)
for(j in 1:N){
  y=x[j]+runif(1,-1,1)
  ##decide which value to use as the next step in the chain
  ##rho=min{(f(y)/f(x))*(q(x)/q(y))}
  rho=min((f(y)/f(x[j])),1)
  if(rho==1){x[j+1]=y}
  if(rho!=1){coin=rbinom(1,1,rho)
  if(coin==1){x[j+1]=y}
  if(coin==0){x[j+1]=x[j]}}}
x[1:10];hist(x,breaks=25)
###NOTE: to get rid of the initializers effect, 
#it is common to have a burn in period. 
#(start the chain after a certain number of realizations have passed.)
x.burn=x[5001:N]
hist(x.burn,breaks=25)
####We can also simulate mutidimensional distributions
#e.g simulate a bivariate normal ((0,0), Sigma)
mu=c(0,0); sigma=matrix(c(1,0.5,0.5,1),2,2)
sigma

f=function(x){z=exp(-t(x)%*%solve(sigma)%*%x/2);return(z)}
f(c(1,2))
##implement : random walk MH
#intialize: 
N=10000
x=matrix(0,2,N); 
x[,1]=1
##simulate a candidate
##proposal density: U(0,5)
for(j in 1:(N-1)){
  y=x[,j]+c(runif(1,-1,1),runif(1,-1,1))
  ##decide which value to use as the next step in the chain
  ##rho=min{(f(y)/f(x))*(q(x)/q(y))}
  rho=min((f(y)/f(x[,j])),1)
  if(rho==1){x[,j+1]=y}
  if(rho!=1){coin=rbinom(1,1,rho)
  if(coin==1){x[,j+1]=y}
  if(coin==0){x[,j+1]=x[,j]}}}

hist(x[1,])
hist(x[2,])

x[1:10];hist(x,breaks=25)
###NOTE: to get rid of the initializers effect, 
#it is common to have a burn in period. 
#(start the chain after a certain number of realizations have passed.)
x.burn=x[5001:N]
hist(x.burn,breaks=25)
##### SINGLE COMPONENT MH
###single component MH:

##Simulate from a Beta-binomial distribution using SCMH
alpha=2;beta=1;n=5
##proposals used in lecture notes: discrete unif (0,...,n) --> 1st comp, unif--> 2nd comp.

##f(x|y)--> Binomial(n,y)
##f(y|x)--> Beta(x+alpha,n-x+beta)

N=20000 ##length of the chain
##initialize
x=c();y=c();x[1]=1; y[1]=0.5
for(j in 2:N){z1=sample((0:n),1);z2=runif(1)
##calculate acceptance probabilities rho1 and rho 2
num=(y[j-1]^z1)*(1-y[j-1])^(n-z1)
den=(y[j-1]^x[j-1])*(1-y[j-1])^(n-x[j-1])
rho1=num/den
num2=(z2^(x[j-1]+alpha-1))*((1-z2)^(n-x[j-1]+beta-1))
den2=(y[j-1]^(x[j-1]+alpha-1))*((1-y[j-1])^(n-x[j-1]+beta-1))
rho2=num2/den2
if(rho1>=1){x[j]=z1}else{coin1=rbinom(1,1,rho1)
if(coin1==1){x[j]=z1}else{x[j]=x[j-1]}}
if(rho2>=1){y[j]=z2}else{coin2=rbinom(1,1,rho2)
if(coin2==1){y[j]=z2}else{y[j]=y[j-1]}}}
x[1:100];y[1:100]
###burn in period
pd=5000
z=cbind(x,y)
dim(z)
z.burn=z[((pd+1):N),]
hist(z.burn[,2])
nrow(z);ncol(z);length(pd+1:n)

### Homeworks:
##1
## Q.4 Generate a uniform, convert to a bernoulli then to a binomial
p=0.25;n=1000
U=runif(3,0,1)
Bernoulli=(U<p)
Bin=sum(Bernoulli) #this is one binomial random variable
#we use a for loop to generate 1000 such variables
Binomial=c() #open column vector
for (i in 1:1000){
  U=runif(3,0,1)
  Bernoulli=(U<p)
  Binomial[i]<-sum(Bernoulli)
}
hist(Binomial)

##Q5. Y=-beta*Sum from i to alpha(log Uj), then Y~Ga(alpha, beta). For question 5, we set 
#a=4,and Beta=1
alpha=4;beta=1;n=1000
U=matrix(runif(n*alpha),alpha,n)
E=-log(U) ##exp(1)
G=beta*colSums(E)
hist(G)

## Q7For this number, we proceed as follows: We generate U~U(0,1) then x=1/lambda*log(1-U) 
# if we set lambda=1/2, then x~EXP(1/2)~Ga(1,1/2). WE also know that k*Ga(a,b)~Ga(ka,b)
#so set k=1/2 and you a Ga(1/2,1/2) Variable. Which is what we need to generate B(1/2,1/2).
lambda=1/2;k=1/2;n=1000
U=runif(n)
E=-(1/lambda)*log(U) #1000 realizations of Exponential(lambda=1/2)
hist(E)
Ga=k*E #This is Gamma(1/2,1/2)
Gb=k*E #This is also Gamma(1/2,1/2)
o=Ga/(Ga+Gb)
hist(o)

##Q8, Box Mueller method.Generate X~N(0,1), use mu+x*sigma to get N(mu, sigma)
mu=3
sig_squared=4
ptm=proc.time()
sims=2000; ns=ceiling(sims/2)
u1=runif(ns,0,1);u2=runif(ns,0,1)
x1=(sqrt(-2*log(u1)))*cos(2*pi*u2)
x2=(sqrt(-2*log(u1)))*sin(2*pi*u2)
x=c(x1,x2)
t=x*sqrt(sig_squared)+mu
hist(t)
time=proc.time()-ptm
time# elaspsed time for this is 0.22 seconds
#Masarglia's polar method. Generate u1,u2~U(-1,1), s=u1^2+u2^2<1
#x=sqrt(-2*log(s))*us=u1/sqrt(s)~N(0,1)
ptm_=proc.time()
x1=list(); x2=list()
Z=0
for (i in 1:1000){
  u1=runif(1,-1,1); u2=runif(1,-1,1);
  s=u1^2+u2^2
  if(s<1){
    x1[[i]]=(sqrt(-2*log(s)))*(u1/sqrt(s))
    x2[[i]]=(sqrt(-2*log(s)))*(u1/sqrt(s))
    Z=Z+1
  }
}
normalx1=unlist(x1); normalx2=unlist(x2)
Y=c(normalx1,normalx2)
Y_=Y*sqrt(sig_squared)+mu
time_=proc.time()-ptm_
time_  #elapsed time for this is 0.33 seconds

#maybe my implementation of the Masarglia's polar method is not very efficient because it takes
#longer,



## 2
## Question 1c.This is the inverse method,generate a Uniform~0,1, 
n=100
inv_vector=vector(length=n)
j=1
while (j<n+1){
  u=runif(1,0,1)
  if (u<1/2){y=sqrt(2*runif(1,0,1))-1}  # as long as U<1/2, then y=the first condition,if not, set
  #y=2nd Case.
  else{y=1-sqrt(2*runif(1,0,1))}
  inv_vector[j]=y
  j=j+1
  
}
inv_vector #print the 100 values

## Question 1d.## The AR method.
AR_vector=vector(length=n)
i=1
while (i<n+1){
  repeat{
    y=runif(1,-1,1)
    G=1+y
    G2=1-y
    u=runif(1,0,1)
    if (u<G & u<G2)
      break
  }
  AR_vector[i]=y
  i=i+1
}
AR_vector

##Q2a. This code is adopted from Professor Abhishek's lab session, MC simulations with n=10^4
sims=10000
x=rnorm(sims,0,1)
mc=c();v=c(); upper=c();lower=c()
indicator=function(x){
  y=ifelse((x>3),1,0) #1 if x>3, 0 otherwise.
  return(y)
}
for (j in 1:sims){
  mc[j]=mean(indicator(x[1:j]))
  v[j]=(j^(-1))*var(indicator(x[1:j]))
  upper[j]=mc[j]+1.96*sqrt(v[j])
  lower[j]=mc[j]-1.96*sqrt(v[j])
}
library(ggplot2)
variables=c(mc,upper,lower)
plot=c(rep("approximating P(Z>3)",sims),rep("upper_bound",sims),rep("lower_bound",sims))
iterations=rep(seq(1:sims),3)
data=data.frame(x=variables,gd=plot,sample=iterations)
graph=ggplot(data,aes(sample,x))+ geom_point(na.rm=TRUE,
                                             aes(colour=factor(gd)))
print(graph + ggtitle("MC approx.of P(Z>3) from a N(0,1) "))

#Q2.B #Shifted exponential and the importance sampler estimate, estimate with n=10^4
n=10000
mc=c();v=c(); upper=c();lower=c()
u=runif(n,0,1) # generate a uniform, then from the uniform generate an exponential
exponential=-(log(u)) #exponential 
x=exponential+3 #shifted exponential,by 3
for(i in 1:n){
  mc[i]=mean((2*pi)^(-1/2)*(exp (((-x[1:i]^2)/2) + x[1:i] - 3))) #h(x)/g(x), normal/exp
  v[i]=(i^{-1})*var((mc[1:i]))
  upper[i]=mc[i]+1.96*sqrt(v[i])
  lower[i]=mc[i]-1.96*sqrt(v[i])}

library(ggplot2)
variables=c(mc,upper,lower)
type=c(rep("P(Z>3) approximation",n),rep("upper_bound",n),rep("lower_bound",n))
iterations=rep(seq(1:n),3)
data=data.frame(p=variables,type=type,sample_size=iterations)
graph2=ggplot(data,aes(sample_size,p))+geom_point(na.rm=TRUE,
                                                  aes(colour=factor(type)))
print(graph2+ggtitle("Importance Sampling Estimate of P(z>3)via ~N(0,1), using a shifted exp"))
##############################################################################################
#Q3a.X~Cauchy(0,1): What we want p=P(0<Z<2), we can use the normal_sim function that we built
#in earlier labs to simulate two normals, then:X~cauchy(0,1) using MC approx
##X1~N(0,1); X2~N(0,1); X1/X2~Cauchy(0,1)
sims=1000; ns=ceiling(sims/2)
normal_sim=function(sims){
  ns=ceiling(sims/2);  u1=runif(ns,0,1);u2=runif(ns,0,1)
  x1=(sqrt(-2*log(u1)))*cos(2*pi*u2)
  x2=(sqrt(-2*log(u1)))*sin(2*pi*u2)
  z=c(x1,x2)[1:sims]; return(z)}
x1=normal_sim(sims);x2=normal_sim(sims)
x=x1/x2 #this is cauchy

mc=c();v=c();upper=c();lower=c()
indicator=function(x){
  y=ifelse((x>0 & x<2),1,0)
  return(y)
}
for (j in 1:sims){
  mc[j]=mean(indicator(x[1:j]))
  v[j]=(1/j)*var(indicator(x[1:j]))
  upper[j]=mc[j]+1.96*sqrt(v[j])
  lower[j]=mc[j]-1.96*sqrt(v[j])
}
library(ggplot2)
values=c(mc,upper,lower)
type=c(rep("estimate of p=P(0<Z<2)",sims),rep("upper
bound",sims),rep("lower bound",sims))
iter=rep(seq(1:sims),3)
data=data.frame(p=values,type=type,sample_size=iter)
graph_3=ggplot(data,aes(sample_size,p))+geom_point(na.rm=TRUE,
                                                   aes(colour=factor(type)))
print(graph_3+ggtitle("MC approximation of p=P(0<Z<2) via Cauchy(0,1)"))

##Q3b use U~(0,2) as an instrument, to compute the importance Sampler, n=10^3

x=runif(sims,0,2)
cauchy_pdf=function(x){
  y=1/(pi*(1+x^2))
  return (y)
}
#looks a bit similar to lab 5, change the range of the z tail
indicator=function(x){
  z=ifelse((x>0 & x<2),1,0)
  return(z)
}
mc=c();v=c();upper=c();lower=c()
for (j in 1:sims){
  mc[j]=mean(indicator(x[1:j])*2*cauchy_pdf(x[1:j]))
  v[j]=(j^{-1})*var(indicator(x[1:j])*2*cauchy_pdf(x[1:j]))
  upper[j]=mc[j]+1.96*sqrt(v[j])
  lower[j]=mc[j]-1.96*sqrt(v[j])
}

library(ggplot2)
variables=c(mc,upper,lower)
type=c(rep("estimate of p=P(0<Z<2)",sims),rep("upper_bound",sims),rep("lower_bound",sims))
iterations=rep(seq(1:sims),3)
data=data.frame(p=variables,type=type,sample_size=iterations)
graph_4=ggplot(data,aes(sample_size,p))+geom_point(na.rm=TRUE,
                                                   aes(colour=factor(type)))
print(graph_4 + ggtitle("Importance Sampling estimate of p=P(0<Z<2) from cauchy(0,1)
using U(0,2)"))
###############################################################################
##Q4b generate from X in part a above
sims=1000
normal_=rnorm(sims,2,9)
x=(3*normal_-1)
mc=c(); v=c(); upper=c(); lower=c()
for (j in 1:sims){
  mc[j]=mean(x[1:j])
  v[j]=(sims^(-1))*var(x[1:j])
  upper[j]=mc[j]+1.96*sqrt(v[j])
  lower[j]=mc[j]-1.96*sqrt(v[j])
}
library(ggplot2)
variables=c(mc, upper, lower)
type=c(rep("E[3X-1]", sims),rep("upper_bound",sims),rep("lower_bound", sims))
iterations=rep(seq(1:sims),3)
data=data.frame(EV=variables, tp=type, sample_size=iterations)
plot5=ggplot(data,aes(sample_size,EV))+
  geom_point(na.rm=TRUE,aes(colour=factor(tp)))
print(plot5 + ggtitle("MC estimate of E[3X-1] for"))

###############################################################################################
#p1=57/76,(1-p1)=19/76,p2=43/76, (1-p2)=33/76, n=76
#Q5b.Using a multinomial to examine association between the variables.
LR=c()
sims=1000
n=76
for (i in 1:sims){
  p1=57/n
  #(1-p1)=19/76 # this assignment is not valid
  p2=43/n
  #(1-p2)=33/76 #this assignment is not valid
  probabilities=c(p1*p2, p1*(1-p2),p2*(1-p1),(1-p1)*(1-p2))
  y=rmultinom(1,n,probabilities) #1 realization of size and probability=probabilities
  p1_hat=(y[1]+y[2])/n;p2_hat=(y[1]+y[3])/n # under the null model
  p_hat=c(p1_hat*p2_hat,p1_hat*(1-p2_hat), (1-p1_hat)*p2_hat,(1-p1_hat)*(1-p2_hat))
  num=-2*sum(y*log(p_hat))#likelihood ratio statistic numerator
  p_estimated=y/n
  den=-2*sum(y*log(p_estimated)) #denominator
  LR[i]=num-den
}
cutoff=quantile(LR,0.95)
cutoff
#########################################################################################
#Q6.Evaluate P(X>5) a truncated normal, using a shifted exponential.The question does not give 
#the number of simulations that we should use, I'll try different values.
sims=1000
u=runif(sims,0,1)
exponential=-(log(u)) #this is exponential
x=exponential+1
normal_=function(x){
  z=exp(-(x^2)/2)
  return (z)
}
indicator=function(x){
  t=ifelse((x>1.5),1,0)
  return (t)
}
m=c(); v=c(); upper=c(); lower=c()
for(i in 1:sims){
  m[i]=sum(indicator(x[1:i])*(normal_(x[1:i])/exp(1-x[1:i])))/sum(normal_(x[1:i])/exp(1-x[1:i]))
  v[i]=(i^(-1))*var((m[1:i]))
  upper[i]=m[i]+1.96*sqrt(v[i])
  lower[i]=m[i]-1.96*sqrt(v[i])                                                          
}
library(ggplot2)
values=c(m,upper,lower)
type=c(rep("P(X>1.5)",sims),rep("upper",sims),rep("lower",sims))
iter=rep(seq(1:sims),3)
data=data.frame(p=values,tp=type,sample_size=iter)
plot6=ggplot(data,aes(sample_size,p))+geom_point(na.rm=TRUE,aes(colour=
                                                                  factor(tp)))
print(plot6 + ggtitle("Evaluation of P(X>1.5) for X ~truncated normal
using shifed exponential"))

##############################################################################
#Question 7a Cholesky decomposition to generate 1000 trivariate,mu=(-2,4,3)
sims=1000
V_mat=matrix(nrow=sims,ncol=3)
for (i in 1:sims){
  cov_mat=matrix(c(2,-1,0.5,-1,4,1,0.5,1,5),3,3)
  mu_mat=matrix(c(-2,4,3))
  z=matrix(c(rnorm(3,0,1)))
  Dec=chol(cov_mat)
  V_=mu_mat+Dec%*%z
  V_mat[i,]<-t(V_)
}
mu_mat=colMeans(V_mat)
variance=matrix(c(var(V_mat[,1]),var(V_mat[,2]),var(V_mat[,3])))
mu_mat
variance
##  Q7b. obtain an MC approximation of the variance-covariance matrix
sims1=1000
V=matrix(nrow=sims1,ncol=3)
for (i in 1:sims1){
  cov_mat1=matrix(c(2,-1,0.5,-1,4,1,0.5,1,5),3,3)
  mu=matrix(c(-2,4,3))
  z=matrix(c(rnorm(3,0,1)))
  dec=chol(cov_mat1)
  Vs=mu+dec%*%z
  V[i,]<-t(Vs) 
}
a=(V[1,]-mu)%*%t(V[1,]-mu)
index=vector("list",length=sims1)
for (i in 1:sims1){
  index[[i]]=(1/(sims1-1))*(V[i,]-mu)%*%t(V[i,]-mu)
}
func_=function(x) {Reduce("+", x)}
func_(index)

##4
#Metropolis Hastings Algorithm

MH=function(n=1000,sigma=0.025,x0=-1){ #Defaulted with the first values of n,sigma and x0
  x=vector(length=n)
  x[1]=x0
  for (i in 1:(n-1)){
    y=rnorm(1,x[i],sigma)
    u=runif(1,0,1)
    rho=min(1,(exp(-y^2))*(2+sin(5*y)+sin(2*y))/(exp(-x[i]^2)*(2+sin(5*x[i])+sin(2*x[i]))))
    x[i+1][u<=rho]=y
    x[i+1][u>rho]=x[i]}
  x
}
y1=MH()
y2=MH(1000,1,-1)
y3=MH(1000,50,-1)

den=function(x){exp(-x^2)*(2+sin(5*x)+sin(2*x))}
int=integrate(den,lower=-Inf,upper=Inf)$value
s=seq(-2,2,.01)
fs=exp(-s^2)*(2+sin(5*s)+sin(2*s))/int

num=c(1:1000)
hist(y1,freq=FALSE,ylim=c(0,1), xlim =c(-2,2),main="Hist and density of f(x),sigma=0.025 & x0=-1")
lines(s,fs)

num=c(1:1000)
hist(y2,freq=FALSE,ylim=c(0,1), xlim =c(-2,2),main="Hist and density of f(x),sigma=1 & x0=-1")
lines(s,fs)

num=c(1:1000)
hist(y3,freq=FALSE,ylim=c(0,1), xlim =c(-2,2), main="Hist and density of f(x),sigma=50 & x0=-1")
lines(s,fs)

## MH1
MH1=function(n=15000,sigma=0.2,x0=-1){ #defaulted to 15000,0.2 and -1
  x=vector(length=n)
  x[1]=x0
  for (i in 1:(n-1)){
    y=rnorm(1,x[i],sigma)
    u=runif(1,0,1)
    rho=min(1,(exp(-y^2))*(2+sin(5*y)+sin(2*y))/(exp(-x[i]^2)*(2+sin(5*x[i])+sin(2*x[i]))))
    x[i+1][u<=rho]=y
    x[i+1][u>rho]=x[i]}
  x
}
y1=MH1()
y2=MH1(x0=-3)
y3=MH1(x0=5)

den=function(x){exp(-x^2)*(2+sin(5*x)+sin(2*x))} # denominator
int=integrate(den,lower=-Inf,upper=Inf)$value
s=seq(-2,2,.01)
fs=exp(-s^2)*(2+sin(5*s)+sin(2*s))/int


num=c(1:15000)
hist(y1[2000:15000],freq=FALSE,ylim=c(0,1), xlim =c(-2,2),main="Hist and density of f(x),sigma=0.2 & x0=-1")
lines(s,fs)


num=c(1:15000)
hist(y2[2000:15000],freq=FALSE,ylim=c(0,1), xlim =c(-2,2),main="Hist and density of f(x),sigma=0.2 & x0=-3")
lines(s,fs)

num=c(1:15000)
hist(y3[2000:15000],freq=FALSE,ylim=c(0,1), xlim =c(-2,2), main="Hist and density of f(x),sigma=0.2 & x0=5")
lines(s,fs)


MH_lag=function(n=1000,sigma=50,x0=-1,lag=1){ ## lag defaulted to 1
  x=vector(length=n)
  l=c()
  x[1]=x0
  p=seq(1,n,lag+1)
  for (i in 1:(n-1)){
    y=rnorm(1,x[i],sigma)
    u=runif(1,0,1)
    rho=min(1,(exp(-y^2))*(2+sin(5*y)+sin(2*y))/(exp(-x[i]^2)*(2+sin(5*x[i])+sin(2*x[i]))))
    x[i+1][u<=rho]=y
    x[i+1][u>rho]=x[i]}
  l=x[p]
}
y1=MH_lag(n=2000)
y2=MH_lag(n=11000,lag=10)
y1=MH_lag(n=101000,lag=100)

den=function(x){exp(-x^2)*(2+sin(5*x)+sin(2*x))} # denominator
int=integrate(den,lower=-Inf,upper=Inf)$value
s=seq(-2,2,0.01)
fs=exp(-s^2)*(2+sin(5*s)+sin(2*s))/int

num=c(1:1000)  ## Lag=1.
hist(y1[200:1000],freq=FALSE,ylim=c(0,1), xlim =c(-2,2),main="Hist and density of f(x),lag=1")
lines(s,fs)

num=c(1:1000) #Lag=10
hist(y2[200:1000],freq=FALSE,ylim=c(0,1), xlim =c(-2,2),main="Hist and density of f(x),lag=10")
lines(s,fs)

num=c(1:1000) # Lag=100
hist(y3[200:1000],freq=FALSE,ylim=c(0,1), xlim =c(-2,2), main="Hist and density of f(x),lag=100")
lines(s,fs)

## Truncated Normal
trn=function(n,d){
  x=vector(length=n)
  z=vector(length=n)
  x[1]=4
  z[1]=0.03
  for (i in 2:n){
    z[i]=runif(1,0,exp(-((x[i-1]^2)/2)))
    x[i]=runif(1,d,sqrt(-2*log(z[i])))
  }
  x    
}
x1=trn(5000,2) #n=5000, d=2
x1_bar=mean(x1[2001:5000]) #burnout=2000 
var_x1=var(x1[2001:5000])
x1_bar;var_x1

x2=trn(10000,3)# n=10000, d=3
x2_bar=mean(x2[2001:5000])
var_x2=var(x2[2001:5000])

x2_bar;var_x2
## Gibbs Sampler
gibbs=function (n){
  x=matrix(ncol=2,nrow = n)
  rho=0.6
  x[1,1]=1
  x[1,2]=1
  for (i in 2:n) {
    x[i,2]=rnorm(1,rho*x[i-1,1],sqrt(1-rho^2))
    x[i,1]=rnorm(1,rho*x[i,2],sqrt(1-rho^2))
  }
  x
}
y=gibbs(15000);mean(y);var(y) # the mean is a bout 0.00xxxxx, which is not as close to zero as I would like.

#The variance covariance matrix however estimates the variance the original covariance matrix better.


## Inverse Gamma.
IG=function(N,a,b,mu,tau_sq,x){
  n=length(x)
  xbar=mean(x)
  theta=matrix(ncol=2,nrow = N)
  theta[1,1]=0
  theta[1,2]=1000
  
  for(i in 2:N){
    theta[i,1]=rnorm(1,((n*tau_sq*xbar)+mu)/((n*tau_sq)+1),sqrt((theta[i-1,2]*tau_sq)/((n*tau_sq)+1)))
    theta[i,2]=1/rgamma(1,((n+1)/2)+a,1/(((sum((x-theta[i,1])^2))/2)+(((theta[i,1]-mu)^2)/(2*tau_sq))+(1/b)))
  }
  theta
  return(mu)
}
N=5000;a=2;b=1;mu=0;tau_sq=10000
x=c(910,1089,857,1713,504,970,557,1043,609,929,693,1384,727,764,1195,803)
theta=IG(N,a,b,mu,tau_sq,x);mean(theta) 















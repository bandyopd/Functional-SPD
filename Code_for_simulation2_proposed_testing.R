library(fda.usc)
library(expm)
library(tmvtnorm)
library(matrixcalc)

source('C:/Users/USER/Downloads/Functions_for_simulation2.R')

f1=function(x)
{
  A=matrix(c(1.5,exp(-sin(pi*x)),exp(-2*sin(pi*x)),exp(-sin(pi*x)),1.5,exp(-sin(pi*x)),exp(-2*sin(pi*x)),exp(-sin(pi*x)),1.5),3,3)
  return(A%*%t(A))
}

f2=function(x)
{
  A=matrix(c(log(2+x),0.1,0.2,0.1,log(2+2*x^2),0.1,0.2,0.1,log(2+3*x^3)),3,3)
  return(A%*%t(A))
}

beta1=function(s)
{
  A=matrix(c(0,sin(pi*s),sin(2*pi*s),sin(pi*s),0,sin(pi*s),sin(2*pi*s),sin(pi*s),0),3,3)
  return(A)
}

n=50 # change this for other n
p=1
d=2
T=2
sigma=matrix(c(1,0.5,0.5,1),d,d)
Y_domain <- seq(0, 1/2, length = T)
M=100
statistics=array(,dim=c(M,p,T,9))
test_result=array(,dim=c(M,p,T))

for(m in 1:M)
{
  print(m)
  
  set.seed(m)
  
  # data generation
  
  e=array(,dim=c(n,T,3,3))
  for(i in 1:n)
  {
    gp=rproc2fdata(6,t=Y_domain,sigma="vexponential",par.list=list("scale"=0.01))$data
    e[i,,1,1]=gp[1,]
    e[i,,2,1]=gp[2,]
    e[i,,3,1]=gp[3,]
    e[i,,2,2]=gp[4,]
    e[i,,3,2]=gp[5,]
    e[i,,3,3]=gp[6,]
    e[i,,1,2]=e[i,,2,1]
    e[i,,1,3]=e[i,,3,1]
    e[i,,2,3]=e[i,,3,2]
  }
  
  Z=rtmvnorm(n = n, mean = rep(0.5, d), sigma = sigma, lower=rep(0, length = d), upper=rep(1, length = d))
  
  X=matrix(as.double(rbinom(n,1,prob=0.25)),n,1)
  
  Y=array(,dim=c(n,T,9))
  for(i in 1:n)
  {
    for(t in 1:T)
    {
      Y[i,t,]=c(X[i,1]*beta1(Y_domain[t])+logm((2+sin(2*pi*Y_domain[t]))*f1(Z[i,1]))+logm((2+cos(2*pi*Y_domain[t]))*f2(Z[i,2]))+e[i,t,,])
    }
  }
  
  # finding optimal bandwidths
  CBS_result=optimal_h_SBF_PL(Z,Z,X,Y,cv=5)
  
  # finding beta_hat, X_tilde, Y_tilde
  SBF_result=SBF_PL_testing(Z,X,Z,X,Y,h_vector=CBS_result$optimal_smoothing)
  beta_hat=SBF_result$betahat
  XY_tilde=SBF_result$ZY_tilde
  
  # for estimating E(epsilon^2)
  Y_hat=array(,dim=c(n,T,9))
  for(i in 1:n)
  {
    Y_hat[i,,]=SBF_PL(matrix(Z[i,],ncol=d),matrix(X[i,],ncol=p),Z[-i,],matrix(X[-i,],ncol=p),Y[-i,,],h_vector=CBS_result$optimal_smoothing)$yhat
  }
  
  # mean of X_tilde product t(X_tilde)
  tilde_prod=matrix(0,p,p)
  for(i in 1:n)
  {
    tilde_prod=tilde_prod+XY_tilde[i,1:p,1,1]%*%t(XY_tilde[i,1:p,1,1])
  }
  
  # computing test statistics
  for(ell in 1:p)
  {
    for(t in 1:T)
    {
      for(v in c(1,4,5,7,8,9))
      {
        statistics[m,ell,t,v]=sqrt(n)*beta_hat[ell,t,v]/sqrt(solve(tilde_prod)[ell,ell]*sum((Y[,t,v]-Y_hat[,t,v])^2))
      }
    }
  }
  
  # multiple testing
  for(ell in 1:p)
  {
    max_value=c()
    for(t in 1:T)
    {
      max_value[t]=max(abs(na.omit(statistics[m,ell,t,])))
      if(max_value[t]>qnorm(1-0.025/6)) test_result[m,ell,t]=1
      else test_result[m,ell,t]=0
    }
  }
}

# number of rejects
sum(test_result[,1,1])
sum(test_result[,1,2])
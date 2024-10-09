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
n0=100
d=2
T=2
sigma=matrix(c(1,0.5,0.5,1),d,d)
Y_domain <- seq(0, 1/2, length = T) # S={0,1/2}
M=100
mean_error=c()

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
  
  e_test=array(,dim=c(n0,T,3,3))
  for(i in 1:n0)
  {
    gp=rproc2fdata(6,t=Y_domain,sigma="vexponential",par.list=list("scale"=0.01))$data
    e_test[i,,1,1]=gp[1,]
    e_test[i,,2,1]=gp[2,]
    e_test[i,,3,1]=gp[3,]
    e_test[i,,2,2]=gp[4,]
    e_test[i,,3,2]=gp[5,]
    e_test[i,,3,3]=gp[6,]
    e_test[i,,1,2]=e_test[i,,2,1]
    e_test[i,,1,3]=e_test[i,,3,1]
    e_test[i,,2,3]=e_test[i,,3,2]
  }
  
  Z_test=rtmvnorm(n = n0, mean = rep(0.5, d), sigma = sigma, lower=rep(0, length = d), upper=rep(1, length = d))
  
  X_test=matrix(as.double(rbinom(n0,1,prob=0.25)),n0,1)
  
  Y_test=array(,dim=c(n0,T,9))
  for(i in 1:n0)
  {
    for(t in 1:T)
    {
      Y_test[i,t,]=c(X_test[i,1]*beta1(Y_domain[t])+logm((2+sin(2*pi*Y_domain[t]))*f1(Z_test[i,1]))+logm((2+cos(2*pi*Y_domain[t]))*f2(Z_test[i,2]))+e_test[i,t,,])
    }
  }
  
  # prediction
  
  CBS_result=optimal_h_SBF_PL(Z_test,Z,X,Y,cv=5)
  
  SBF_result=SBF_PL(Z_test,X_test,Z,X,Y,h_vector=CBS_result$optimal_smoothing)
  
  Y_hat=SBF_result$yhat
  
  Y_test_matrix=array(,dim=c(n0,T,3,3))
  Y_hat_matrix=array(,dim=c(n0,T,3,3))
  for(i in 1:n0)
  {
    for(t in 1:T)
    {
      Y_test_matrix[i,t,,]=matrix(Y_test[i,t,],3,3)
      Y_hat_matrix[i,t,,]=matrix(Y_hat[i,t,],3,3)
    }
  }
  
  # evaluation
  
  error=c()
  for(i in 1:n0)
  {
    distance=c()
    for(t in 1:T)
    {
      distance[t]=frobenius.norm(Y_hat_matrix[i,t,,]-Y_test_matrix[i,t,,])
    }
    error[i]=sum(distance^2)
  }
  mean_error[m]=mean(error)
}

mean(mean_error)
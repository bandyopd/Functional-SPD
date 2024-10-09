library(pracma)
library(fda.usc)
library(expm)
library(tmvtnorm)
library(matrixcalc)

source('C:/Users/USER/Downloads/Functions_for_simulation1.R')

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

n=50 # change this for other n
n0=100 # size of test set
d=2
T=101 # number of grid for S=[0,1]
sigma=matrix(c(1,0.5,0.5,1),d,d) # covariance matrix for generating Z
Y_domain <- seq(0, 1, length = T) # grid for S=[0,1]
M=100 # number of Monte-Carlo samples
mean_error=c()

for(m in 1:M)
{
  print(m)
  
  set.seed(m)
  
  # generate epsilon
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
  
  # generate Z
  Z=rtmvnorm(n = n, mean = rep(0.5, d), sigma = sigma, lower=rep(0, length = d), upper=rep(1, length = d))
  
  # generate Y
  Y=array(,dim=c(n,T,9))
  for(i in 1:n)
  {
    for(t in 1:T)
    {
      Y[i,t,]=c(logm((2+sin(2*pi*seq(0,1,length=T)[t]))*f1(Z[i,1]))+logm((2+cos(2*pi*seq(0,1,length=T)[t]))*f2(Z[i,2]))+e[i,t,,])
    }
  }
  
  # generate epsilon_test
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
  
  # generate Z_test
  Z_test=rtmvnorm(n = n0, mean = rep(0.5, d), sigma = sigma, lower=rep(0, length = d), upper=rep(1, length = d))
  
  # generate Y_test
  Y_test=array(,dim=c(n0,T,9))
  for(i in 1:n0)
  {
    for(t in 1:T)
    {
      Y_test[i,t,]=c(logm((2+sin(2*pi*seq(0,1,length=T)[t]))*f1(Z_test[i,1]))+logm((2+cos(2*pi*seq(0,1,length=T)[t]))*f2(Z_test[i,2]))+e_test[i,t,,])
    }
  }
  
  # finding optimal bandwidth via 5-fold cross-validation
  CBS_result=optimal_h_SBF(Z_test,Z,Y,Y_domain,cv=5) 

  # fitting with optimal bandwidth (CBS_result$optimal_smoothing)
  SBF_result=SBF(Z_test,Z,Y,Y_domain,h_vector=CBS_result$optimal_smoothing)

  # get predicted Y and convert its form
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
  
  # computing APE_m
  error=c()
  for(i in 1:n0)
  {
    distance=c()
    for(t in 1:T)
    {
      distance[t]=frobenius.norm(Y_hat_matrix[i,t,,]-Y_test_matrix[i,t,,])
    }
    error[i]=trapz(seq(0,1,length=T),distance^2)
  }
  mean_error[m]=mean(error)
}

mean(mean_error) # computing APE

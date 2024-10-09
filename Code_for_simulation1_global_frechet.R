library(pracma)
library(fda.usc)
library(expm)
library(tmvtnorm)
library(matrixcalc)

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
n0=100
d=2
T=101
sigma=matrix(c(1,0.5,0.5,1),d,d)
Y_domain <- seq(0, 1, length = T)
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
  
  Y_sym=array(,dim=c(n,T,3,3))
  for(i in 1:n)
  {
    for(t in 1:T)
    {
      Y_sym[i,t,,]=logm((2+sin(2*pi*seq(0,1,length=T)[t]))*f1(Z[i,1]))+logm((2+cos(2*pi*seq(0,1,length=T)[t]))*f2(Z[i,2]))+e[i,t,,]
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
  
  Y_sym_test=array(,dim=c(n0,T,3,3))
  for(i in 1:n0)
  {
    for(t in 1:T)
    {
      Y_sym_test[i,t,,]=logm((2+sin(2*pi*seq(0,1,length=T)[t]))*f1(Z_test[i,1]))+logm((2+cos(2*pi*seq(0,1,length=T)[t]))*f2(Z_test[i,2]))+e_test[i,t,,]
    }
  }
  
  # prediction
  
  Z_center=Z
  for(j in 1:d)
  {
    Z_center[,j]=Z[,j]-mean(Z[,j])
  }
  cov_hat=matrix(0,d,d)
  for(i in 1:n)
  {
    cov_hat=cov_hat+matrix(Z_center[i,],ncol=1)%*%matrix(Z_center[i,],nrow=1)/n
  }
  cov_hat_inv=solve(cov_hat)
  Y_sym_hat=array(,dim=c(n0,T,3,3))
  for(k in 1:n0)
  {
    weight=c()
    for(i in 1:n)
    {
      weight[i]=1+(matrix(Z_center[i,],nrow=1)%*%cov_hat_inv)%*%matrix(Z_test[k,]-colMeans(Z),ncol=1)
    }
    lower=sum(weight)
    for(t in 1:T)
    {
      upper=matrix(0,3,3)
      for(i in 1:n)
      {
        upper=upper+weight[i]*Y_sym[i,t,,]
      }
      Y_sym_hat[k,t,,]=upper/lower
    }
  }
  
  # evaluation
  
  error=c()
  for(k in 1:n0)
  {
    distance=c()
    for(t in 1:T)
    {
      distance[t]=frobenius.norm(Y_sym_hat[k,t,,]-Y_sym_test[k,t,,])
    }
    error[k]=trapz(seq(0,1,length=T),distance^2)
  }
  mean_error[m]=mean(error)
}

mean(mean_error)
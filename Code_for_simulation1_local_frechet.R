library(pracma)
library(fda.usc)
library(expm)
library(tmvtnorm)
library(matrixcalc)
library(pdist)

# function for computing weighted mean for log_Euclidean transformed SPD matrices
weighted_mean_LogE_logm_fast <- function(logS,w_vec)
{
  w_Log_S=array(,dim=c(3,3,dim(logS)[3]))
  for(i in 1:3)
  {
    for(j in 1:3)
    {
      w_Log_S[i,j,]=w_vec*logS[i,j,]
    }
  }
  rowSums(w_Log_S,dims=2)
}

initial_h_LogE_logm=function(x.out,X,nfolds)
{
  set.seed(1)
  n=nrow(X)
  s=sample(n)
  X=X[s,]
  folds=cut(1:n,breaks=nfolds,labels=FALSE)
  distance1=c()
  for(k in 1:nfolds)
  {
    X.training=X[-which(folds==k),]
    X.test=X[which(folds==k),]
    new.distance=c()
    for(p in 1:nrow(X.test))
    {
      new.distance[p]=min(as.matrix(pdist(X=X.test[p,],Y=X.training)))
    }
    distance1[k]=max(new.distance)
  }
  distance2=c()
  for(k in 1:nrow(x.out))
  {
    distance2[k]=min(as.matrix(pdist(X=x.out[k,],Y=X)))
  }
  h_initial=max(distance1,distance2)+1
  return(h_initial)
}

# function for computing optimal bandwidth
optimal_h_LogE_logm_fast=function(x.out,X,logY,nfolds,h_add,h_length,h_initial)
{
  set.seed(1)
  n=nrow(X)
  s=sample(n)
  X=X[s,]
  logY=logY[,,s]
  folds=cut(1:n,breaks=nfolds,labels=FALSE)
  h_vector=seq(h_initial,h_initial+h_add,length=h_length)
  error.h=c()
  for(j in 1:h_length)
  {
    error.fold=c()
    for(k in 1:nfolds)
    {
      X.training=X[-which(folds==k),]
      X.test=X[which(folds==k),]
      logY.training=logY[,,-which(folds==k)]
      logY.test=logY[,,which(folds==k)]
      logY.test.hat=array(,dim=c(3,3,nrow(X.test)))
      error.fold[k]=0
      for(p in 1:nrow(X.test))
      {
        logY.test.hat[,,p]=weighted_mean_LogE_logm_fast(logS=logY.training,w_vec=K_weight(X.test[p,],X.training,h_vector[j]))
        error.fold[k]=error.fold[k]+(frobenius.norm(logY.test.hat[,,p]-logY.test[,,p]))^2
      }
    }
    error.h[j]=sum(error.fold)
  }
  return(min(h_vector[which.min(error.h)]))
}

K=function(x)
{
  3/4*(1-x^2)*dunif(x,-1,1)*2
}

K_weight=function(x,X,h)
{
  K_values=c(K(as.matrix(pdist(X=x,Y=X))/h))
  K_values/sum(K_values)
}

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
  
  Y_sym=array(,dim=c(T,3,3,n))
  for(i in 1:n)
  {
    for(t in 1:T)
    {
      Y_sym[t,,,i]=logm((2+sin(2*pi*seq(0,1,length=T)[t]))*f1(Z[i,1]))+logm((2+cos(2*pi*seq(0,1,length=T)[t]))*f2(Z[i,2]))+e[i,t,,]
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
  
  Y_sym_test=array(,dim=c(T,3,3,n0))
  for(i in 1:n0)
  {
    for(t in 1:T)
    {
      Y_sym_test[t,,,i]=logm((2+sin(2*pi*seq(0,1,length=T)[t]))*f1(Z_test[i,1]))+logm((2+cos(2*pi*seq(0,1,length=T)[t]))*f2(Z_test[i,2]))+e_test[i,t,,]
    }
  }
  
  # prediction
  
  Y_sym_hat=array(,dim=c(T,3,3,n0))
  h_initial=initial_h_LogE_logm(Z_test,Z,nfolds=5)
  for(t in 1:T)
  {
    opt_h=optimal_h_LogE_logm_fast(Z_test,Z,Y_sym[t,,,],nfolds=5,h_add=2,h_length=21,h_initial)
    for(k in 1:n0)
    {
      Y_sym_hat[t,,,k]=weighted_mean_LogE_logm_fast(logS=Y_sym[t,,,],w_vec=K_weight(Z_test[k,],Z,opt_h))
    }
  }
  
  # evaluation
  
  error=c()
  for(k in 1:n0)
  {
    distance=c()
    for(t in 1:T)
    {
      distance[t]=frobenius.norm(Y_sym_hat[t,,,k]-Y_sym_test[t,,,k])
    }
    error[k]=trapz(seq(0,1,length=T),distance^2)
  }
  mean_error[m]=mean(error)
}

mean(mean_error)
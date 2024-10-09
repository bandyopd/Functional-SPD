dyn.load('C:/Users/USER/Downloads/CBS_vector_valued_L2_new_grid.dll')
dyn.load('C:/Users/USER/Downloads/SBF_vector_valued_L2_new_grid.dll')

library(fields)
library(matrixStats)

# function for dividing data for cross-validation
new_cut=function(n,cv)
{
  folds=c()
  size=floor(n/cv)
  for(k in 1:cv)
  {
    if(k<cv) folds[((k-1)*size+1):(k*size)]=k
    if(k==cv) folds[((cv-1)*size+1):n]=k
  }
  return(folds)
}

# function for finding the smallest bandwidths that are valid for SBF
smoothing_options=function(W_training,epsilon,g_length,cv,smoothing_length,smoothing_add)
{
  n=nrow(W_training)
  d=ncol(W_training)
  folds=new_cut(n,cv)
  distance=matrix(,nrow=cv,ncol=d)
  for(k in 1:cv)
  {
    W_training_training=matrix(W_training[-which(folds==k),],ncol=d)
    W_training_test=matrix(W_training[which(folds==k),],ncol=d)
    for(j in 1:d)
    {
      grid_add=seq(0,1,length=g_length)
      total_grid=c(W_training_test[,j],grid_add)
      distance[k,j]=max(rowMins(rdist(total_grid,W_training_training[,j])))
    }
  }
  smoothing_matrix=matrix(,nrow=max(smoothing_length),ncol=d)
  for(j in 1:d)
  {
    smoothing_matrix[1:smoothing_length[j],j]=seq(max(distance[,j])+epsilon,max(distance[,j])+epsilon+smoothing_add[j],length=smoothing_length[j])
  }
  return(smoothing_matrix)
}

# function for finding the smallest bandwidths that are valid for SBF
smoothing_options_for_test=function(W_test,W_training,epsilon,g_length)
{
  d=ncol(W_training)
  distance=c()
  for(j in 1:d)
  {
    grid_add=seq(0,1,length=g_length)
    total_grid=c(W_test[,j],grid_add)
    distance[j]=max(rowMins(rdist(total_grid,W_training[,j])))
  }
  return(distance+epsilon)
}

# function for finding optimal bandwidths in additive regression (for the case where S is an interval)
# W_test: n0 by d numeric matrix with values in [0,1] (corresponding to Z_test)
# W_training: n by d numeric matrix with values in [0,1] (corresponding to Z)
# d should be >1
# Y_training: n by T by vec_length numeric array (vec_length is nrow * ncol of SPD matrix)
# cv: number of folds of cross-validation
# output includes optimal bandwidths
optimal_h_SBF=function(W_test,W_training,Y_training,Y_domain,cv)
{
  d_x_1=ncol(W_training)
  
  # some tuning options
  max_sbf_iteration=as.integer(50) # maximum number of SBF iteration
  max_cbs_iteration=as.integer(10) # maximum number of CBS iteration
  epsilon=10^-4 # stopping tolerance of SBF iteration
  smoothing_length=rep(as.integer(21),d_x_1)
  smoothing_add=rep(0.5,d_x_1)
  # seq(min=a_j,max=a_j+smoothing_add[j],length=smoothing_length[j]) will be the jth bandwidth grid, where a_j is defined in the supplement
  g_length=as.integer(101) # number of grid for numerical integration over [0,1] via trapezoidal rule
  
  # do not change below (some objects are for more complex models)
  n0=nrow(W_test)
  n=nrow(W_training)
  
  shuffle=sample(1:n)
  Y_training=Y_training[shuffle,,]
  W_training=W_training[shuffle,]
  
  cv=as.integer(cv) 
  d_x_2=as.integer(0)
  d_x_3=as.integer(0)
  W_training_2=matrix(0,nrow=n,ncol=2*d_x_2)
  W_training_3=matrix(0,nrow=n,ncol=3*d_x_3)
  max_g_length=as.integer(g_length)
  g_length_2=as.integer(0)
  g_length_3=as.integer(0)
  g_2=matrix(0,g_length_2,2)
  g_3=matrix(0,g_length_3,3)
  gg_2=matrix(0,g_length_2,2*d_x_2)
  
  class=rep(as.integer(1),d_x_1)
  class_2=array(as.integer(0),dim=d_x_2)
  axis_length_2=array(0,dim=2*d_x_2)
  class_3=array(as.integer(0),dim=d_x_3)
  axis_length_3=array(0,dim=3*d_x_3)
  
  d_u=as.integer(0)
  d_v=as.integer(0)
  uv_cardinality=array(as.integer(0),dim=d_u+d_v)
  v_values=matrix(0,d_v,max_g_length)
  
  max_smoothing_length=max(smoothing_length)
  cbs_iteration=as.integer(0)
  T=dim(Y_training)[2]
  vec_length=dim(Y_training)[3]
  Y_domain=as.double(Y_domain)
  
  W_test_2=matrix(0,nrow=n0,ncol=2*d_x_2)
  W_test_3=matrix(0,nrow=n0,ncol=3*d_x_3)
  
  minimum_bandwidth_1_for_cbs=smoothing_options(W_training,epsilon,g_length,cv,smoothing_length,smoothing_add)[1,]
  minimum_bandwidth_1_for_sbf=smoothing_options_for_test(W_test,W_training,epsilon,g_length)
  optimal_smoothing=pmax(minimum_bandwidth_1_for_cbs,minimum_bandwidth_1_for_sbf)
  
  smoothing_matrix=matrix(0,nrow=max(smoothing_length),ncol=d_x_1+d_x_2+d_x_3)
  for(j in 1:(d_x_1+d_x_2+d_x_3))
  {
    smoothing_matrix[1:smoothing_length[j],j]=seq(optimal_smoothing[j],optimal_smoothing[j]+smoothing_add[j],length=smoothing_length[j])
  }
  
  cv_test_size=as.integer(floor(n/cv))
  sbf_iteration=array(as.integer(0),dim=c(max_cbs_iteration,d_x_1+d_x_2+d_x_3,max_smoothing_length,cv))
  phat1=array(1,dim=c(d_x_1+d_x_2+d_x_3,max_g_length+n,max_cbs_iteration,d_x_1+d_x_2+d_x_3,max_smoothing_length,cv))
  
  # apply fortran dll file for fitting
  cbs_result=.Fortran('CBS_vector_valued_L2_new_grid',n=n,d_x_1=d_x_1,d_x_2=d_x_2,d_x_3=d_x_3,d_u=d_u,d_v=d_v,T=T,m_pt=vec_length,W_training=W_training,W_training_2=W_training_2,
                      W_training_3=W_training_3,Y_training=Y_training,cv=cv,cv_test_size=cv_test_size,max_g_length=max_g_length,g_length=g_length,
                      g_length_2=g_length_2,g_length_3=g_length_3,g_2=g_2,g_3=g_3,gg_2=gg_2,class=class,class_2=class_2,class_3=class_3,
                      axis_length_2=axis_length_2,axis_length_3=axis_length_3,uv_cardinality=uv_cardinality,v_values=v_values,
                      Y_domain=Y_domain,max_smoothing_length=max_smoothing_length,smoothing_length=smoothing_length,
                      smoothing_matrix=smoothing_matrix,epsilon=epsilon,max_cbs_iteration=max_cbs_iteration,
                      max_sbf_iteration=max_sbf_iteration,cbs_iteration=cbs_iteration,sbf_iteration=sbf_iteration,
                      phat1=phat1,optimal_smoothing=optimal_smoothing)
  
  return(cbs_result)
}

# function for fitting additive models (for the case where S is an interval)
# h_vector: d-vector consists of bandwidths
# output includes predicted Y
SBF=function(W_test,W_training,Y_training,Y_domain,h_vector)
{
  # some tuning options
  max_sbf_iteration=as.integer(50) # maximum number of SBF iteration
  epsilon=10^-4 # stopping tolerance of SBF iteration
  g_length=as.integer(101) # should be the number used in optimal_h_SBF
  
  # do not change below (some objects are for more complex models)
  d_x_1=ncol(W_training)
  n0=nrow(W_test)
  n=nrow(W_training)
  
  d_x_2=as.integer(0)
  d_x_3=as.integer(0)
  W_training_2=matrix(0,nrow=n,ncol=2*d_x_2)
  W_training_3=matrix(0,nrow=n,ncol=3*d_x_3)
  max_g_length=as.integer(g_length)
  g_length_2=as.integer(0)
  g_length_3=as.integer(0)
  g_2=matrix(0,g_length_2,2)
  g_3=matrix(0,g_length_3,3)
  gg_2=matrix(0,g_length_2,2*d_x_2)
  
  class=rep(as.integer(1),d_x_1)
  class_2=array(as.integer(0),dim=d_x_2)
  axis_length_2=array(0,dim=2*d_x_2)
  class_3=array(as.integer(0),dim=d_x_3)
  axis_length_3=array(0,dim=3*d_x_3)
  
  d_u=as.integer(0)
  d_v=as.integer(0)
  uv_cardinality=array(as.integer(0),dim=d_u+d_v)
  v_values=matrix(0,d_v,max_g_length)
  
  T=dim(Y_training)[2]
  vec_length=dim(Y_training)[3]
  Y_domain=as.double(Y_domain)
  
  W_test_2=matrix(0,nrow=n0,ncol=2*d_x_2)
  W_test_3=matrix(0,nrow=n0,ncol=3*d_x_3)
  
  sbf_iteration=as.integer(0)
  mhat=array(0,dim=c(n0,d_x_1+d_x_2+d_x_3,T,vec_length))
  yhat=array(0,dim=c(n0,T,vec_length))
  smoothing_vector=as.vector(h_vector)
  
  sbf_result=.Fortran('SBF_vector_valued_L2_new_grid',n0=n0,n=n,d_x_1=d_x_1,d_x_2=d_x_2,d_x_3=d_x_3,d_u=d_u,d_v=d_v,T=T,m_pt=vec_length,
                      W_test=W_test,W_test_2=W_test_2,W_test_3=W_test_3,W_training=W_training,W_training_2=W_training_2,
                      W_training_3=W_training_3,Y_training=Y_training,max_g_length=max_g_length,g_length=g_length,g_length_2=g_length_2,
                      g_length_3=g_length_3,g_2=g_2,g_3=g_3,gg_2=gg_2,class=class,class_2=class_2,class_3=class_3,axis_length_2=axis_length_2,
                      axis_length_3=axis_length_3,uv_cardinality=uv_cardinality,v_values=v_values,Y_domain=Y_domain,
                      smoothing_vector=smoothing_vector,epsilon=epsilon,max_sbf_iteration=max_sbf_iteration,
                      sbf_iteration=sbf_iteration,mhat=mhat,yhat=yhat)
  
  return(sbf_result)
}
dyn.load('C:/Users/USER/Downloads/CBS_vector_valued_Euclidean_new_grid_PL.dll')
dyn.load('C:/Users/USER/Downloads/SBF_vector_valued_Euclidean_new_grid_PL.dll')
dyn.load('C:/Users/USER/Downloads/SBF_vector_valued_Euclidean_new_grid_PL_tilde_out.dll')

source('C:/Users/USER/Downloads/Functions_for_simulation1.R')

# function for finding optimal bandwidths in partially linear additive regression (for the case where S is finite)
# V_training: n by p numeric matrix (corresponding to X_training)
# output includes optimal bandwidths
optimal_h_SBF_PL=function(W_test,W_training,V_training,Y_training,cv)
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
  d_z=ncol(V_training)
  
  shuffle=sample(1:n)
  Y_training=Y_training[shuffle,,]
  V_training=matrix(V_training[shuffle,],n,d_z)
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
  
  cbs_result=.Fortran('CBS_vector_valued_Euclidean_new_grid_PL',n=n,d_x_1=d_x_1,d_x_2=d_x_2,d_x_3=d_x_3,d_u=d_u,d_v=d_v,d_z=d_z,T=T,m_pt=vec_length,W_training=W_training,W_training_2=W_training_2,
                      W_training_3=W_training_3,Z_training=V_training,Y_training=Y_training,cv=cv,cv_test_size=cv_test_size,max_g_length=max_g_length,g_length=g_length,
                      g_length_2=g_length_2,g_length_3=g_length_3,g_2=g_2,g_3=g_3,gg_2=gg_2,class=class,class_2=class_2,class_3=class_3,
                      axis_length_2=axis_length_2,axis_length_3=axis_length_3,uv_cardinality=uv_cardinality,v_values=v_values,
                      max_smoothing_length=max_smoothing_length,smoothing_length=smoothing_length,
                      smoothing_matrix=smoothing_matrix,epsilon=epsilon,max_cbs_iteration=max_cbs_iteration,
                      max_sbf_iteration=max_sbf_iteration,cbs_iteration=cbs_iteration,sbf_iteration=sbf_iteration,
                      phat1=phat1,optimal_smoothing=optimal_smoothing)
  
  return(cbs_result)
}

# function for fitting partially linear additive models (for the case where S is finite)
# V_test: n0 by p numeric matrix (corresponding to X_test)
# h_vector: d-vector consists of bandwidths
# output includes predicted Y
SBF_PL=function(W_test,V_test,W_training,V_training,Y_training,h_vector)
{
  # some tuning options
  max_sbf_iteration=as.integer(50) # maximum number of SBF iteration
  epsilon=10^-4 # stopping tolerance of SBF iteration
  g_length=as.integer(101) # should be the number used in optimal_h_SBF
  
  # do not change below (some objects are for more complex models)
  d_x_1=ncol(W_training)
  n0=nrow(W_test)
  n=nrow(W_training)
  d_z=ncol(V_training)
  
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
  
  W_test_2=matrix(0,nrow=n0,ncol=2*d_x_2)
  W_test_3=matrix(0,nrow=n0,ncol=3*d_x_3)
  
  sbf_iteration=as.integer(0)
  mhat=array(0,dim=c(n0,d_x_1+d_x_2+d_x_3,T,vec_length))
  yhat=array(0,dim=c(n0,T,vec_length))
  smoothing_vector=as.vector(h_vector)
  betahat=array(0,dim=c(d_z,T,vec_length))
  betahat_0=array(0,dim=c(T,vec_length))
  
  sbf_result=.Fortran('SBF_vector_valued_Euclidean_new_grid_PL',n0=n0,n=n,d_x_1=d_x_1,d_x_2=d_x_2,d_x_3=d_x_3,d_u=d_u,d_v=d_v,d_z=d_z,T=T,m_pt=vec_length,
                      W_test=W_test,W_test_2=W_test_2,W_test_3=W_test_3,Z_test=V_test,W_training=W_training,W_training_2=W_training_2,
                      W_training_3=W_training_3,Z_training=V_training,Y_training=Y_training,max_g_length=max_g_length,g_length=g_length,g_length_2=g_length_2,
                      g_length_3=g_length_3,g_2=g_2,g_3=g_3,gg_2=gg_2,class=class,class_2=class_2,class_3=class_3,axis_length_2=axis_length_2,
                      axis_length_3=axis_length_3,uv_cardinality=uv_cardinality,v_values=v_values,
                      smoothing_vector=smoothing_vector,epsilon=epsilon,max_sbf_iteration=max_sbf_iteration,
                      sbf_iteration=sbf_iteration,mhat=mhat,yhat=yhat,betahat=betahat,betahat_0=betahat_0)
  
  return(sbf_result)
}

# function for testing in partially linear additive regression (for the case where S is finite)
# output includes beta_hat, X_tilde, Y_tilde
SBF_PL_testing=function(W_test,V_test,W_training,V_training,Y_training,h_vector)
{
  # some tuning options
  max_sbf_iteration=as.integer(50) # maximum number of SBF iteration
  epsilon=10^-4 # stopping tolerance of SBF iteration
  g_length=as.integer(101) # should be the number used in optimal_h_SBF
  
  # do not change below (some objects are for more complex models)
  d_x_1=ncol(W_training)
  n0=nrow(W_test)
  n=nrow(W_training)
  d_z=ncol(V_training)
  
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
  
  W_test_2=matrix(0,nrow=n0,ncol=2*d_x_2)
  W_test_3=matrix(0,nrow=n0,ncol=3*d_x_3)
  
  sbf_iteration=as.integer(0)
  mhat=array(0,dim=c(n0,d_x_1+d_x_2+d_x_3,T,vec_length))
  yhat=array(0,dim=c(n0,T,vec_length))
  smoothing_vector=as.vector(h_vector)
  betahat=array(0,dim=c(d_z,T,vec_length))
  betahat_0=array(0,dim=c(T,vec_length))
  ZY_tilde=array(0,dim=c(n,d_z+1,T,vec_length))
  
  sbf_result=.Fortran('SBF_vector_valued_Euclidean_new_grid_PL_tilde_out',n0=n0,n=n,d_x_1=d_x_1,d_x_2=d_x_2,d_x_3=d_x_3,d_u=d_u,d_v=d_v,d_z=d_z,T=T,m_pt=vec_length,
                      W_test=W_test,W_test_2=W_test_2,W_test_3=W_test_3,Z_test=V_test,W_training=W_training,W_training_2=W_training_2,
                      W_training_3=W_training_3,Z_training=V_training,Y_training=Y_training,max_g_length=max_g_length,g_length=g_length,g_length_2=g_length_2,
                      g_length_3=g_length_3,g_2=g_2,g_3=g_3,gg_2=gg_2,class=class,class_2=class_2,class_3=class_3,axis_length_2=axis_length_2,
                      axis_length_3=axis_length_3,uv_cardinality=uv_cardinality,v_values=v_values,
                      smoothing_vector=smoothing_vector,epsilon=epsilon,max_sbf_iteration=max_sbf_iteration,
                      sbf_iteration=sbf_iteration,mhat=mhat,yhat=yhat,betahat=betahat,betahat_0=betahat_0,ZY_tilde=ZY_tilde)
  
  return(sbf_result)
}
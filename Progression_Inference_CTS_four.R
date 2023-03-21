Progression_Inference_CTS_four=function(DATA,frac,grade,os){
  #-----Progression_Inference-------
  ngene=dim(DATA)[2];
  npatient=dim(DATA)[1];
  grade=grade;
  
  #-----Cauculate Gaussian similarity function-----
  cosSim <- function(x,y){ dot_product <- sum(x * y) #计算两个向量的点积 
  norm_x <- sqrt(sum(x^2)) #计算向量x的范数 
  norm_y <- sqrt(sum(y^2)) #计算向量y的范数 
  cos_sim <- dot_product / (norm_x * norm_y) #计算余弦相似度 
  return (cos_sim) }
  
  
  omega=matrix(rep(0,npatient*npatient),ncol=npatient,nrow=npatient);
  gamma=omega;S=omega;cosin=omega;temp=omega;Ks=omega;
  for (i in 1:npatient)
    for(j in 1:npatient){
      omega[i,j]=(1+abs(grade[i]-grade[j]));
    }
  
  k=10; #number of neighbors considered in knn
  library(RANN)
  
  for (x in 1:dim(DATA)[3])
  {
    data=t(DATA[,,x]);
    knnres=nn2(t(data),t(data),k+1); 
    idx_k=knnres$nn.idx[,k+1];
    dists=knnres$nn.dists;
    sigma=dists[,k+1];
    sigma2=sigma^2;
    for (i in 1:npatient)
      for(j in 1:npatient){
        gamma[i,j]=(sigma2[i]+sigma2[j]);
      }
    for (i in 1:npatient)
      for(j in 1:npatient){
        temp[i,j]=(sum((data[,i]-data[,j])^2))/gamma[i,j]; 
        #gaussian similarity function
        
      }
    S=S+temp 
    
  }
  
  for (i in 1:npatient)
    for(j in 1:npatient){
      cosin[i,j]=cosSim(frac[i,],frac[j,])
      S[i,j]=exp(-(omega[i,j])*S[i,j]/cosin[i,j])
      
    }
  
  #-----Cauculate transition matrix-----
  D=rep(0,npatient);
  H=S;E=D;P=S;
  for(i in 1:npatient)
    D[i]=sum(S[i,]);
  for(i in 1:npatient)
    for(j in 1:npatient)
      H[i,j]=S[i,j]/D[i]/D[j];
  for(i in 1:npatient)
  {H[i,i]=0;E[i]=sum(H[i,]);}
  for(i in 1:npatient)
    for(j in 1:npatient)
      P[i,j]=E[i]^(-0.5)*H[i,j]*E[j]^(-0.5); # transition matrix
  
  phi0=E/sqrt(sum(E^2)); # largest eigenvector of P
  phi0=as.matrix(phi0);
  library(MASS)
  Q=ginv((diag(npatient)-P+phi0%*%t(phi0)))-diag(npatient); # accumulated transition matrix
  
  #---Select root patient
  grade=as.numeric(grade);
  Ind_max=subset(1:npatient,grade==max(grade));
  b=length(Ind_max);
  #b=1;
  rn=ceiling(runif(1)*b);
  x_ref=Ind_max[rn];
  
  depth_to_root=function(Q,s){
    nQ=nrow(Q);
    dpt=rep(0,nQ);
    for(i in 1:nQ){
      dpt[i]=sqrt(sum((Q[i,]-Q[s,])^2));
    }
    return(dpt);
  }
  
  TPD_to_xref=depth_to_root(Q,x_ref);
  Ind_min=subset(1:npatient,grade==min(grade));
  Ind_sort=order(TPD_to_xref,decreasing=TRUE);
  
  for(i in 1:length(Ind_sort)){
    if(grade[Ind_sort[i]]==1 & os[Ind_sort[i]]>mean(os)) {
      root=Ind_sort[i];break;
    }
    #select the root with maximum TPD to x_ref in the minimum grade subgroup
  }
  #-----Transfer into a smoothed trajectory----
  #root=58;
  PPD=depth_to_root(Q,root);
  PPD_order=order(PPD);
  
  
  return(list(Accumulated_Transition_Matrix=Q,Temporal_Progression=PPD,Order=PPD_order));
  
}
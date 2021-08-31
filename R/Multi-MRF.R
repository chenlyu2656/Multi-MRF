null.Multi.MRF<-function(Y,distance,X=NULL,correlation=c('exchangeable','distance-based'),out_type=c("normal","beta")) 
{
  ##Preliminary
  Y<-as.matrix(Y);n<-nrow(Y);index.cluster<-unique(distance[,1])
  X1<-cbind(rep(1,n),X);p<-ncol(X1)-1;V<-diag(1,n);
  S.subject<-matrix(0,n,n)
  
  for (j in index.cluster)
  {
    index<-which(distance[,1]==j);ni<-length(index)
    if (correlation == 'exchangeable'){
      S.subject[index,index]<-1
    }
    if (correlation == 'distance-based'){
      dis <- distance[index,2]
      S.dis<-outer(dis,dis,FUN='-');
      td <- 10000 # can be changed depending on total length of the corresponding region for CpG sites; 
      S.dis<-exp(-abs(S.dis/td*3));
      S.subject[index,index]<-S.dis;
    }
    diag(S.subject)<-0
  }
  
  #Iteration to get beta and eta
  step.error<-Inf;beta<-rep(0,p+1);eta<-0;iter<-0 #iteration time
  while (step.error>0.0001)
  {
    step.beta<-Inf;mu<-rep(mean(Y),n);Y.res<-Y-mu;
    while (step.beta>0.00001)
    {
      if (out_type=="normal")
      {delta<-rep(1,n);A<-t(X1)%*%V%*%X1;B<-t(X1)%*%V%*%Y.res}
      if (out_type=="beta")
      {
        delta<-as.vector(mu*(1-mu));
        A<-t(delta^{1/2}*X1)%*%V%*%(delta^{1/2}*X1);
        B<-t(delta^{1/2}*X1)%*%V%*%(delta^{-1/2}*Y.res)
      }
      #calculating residual
      beta.new<-beta+solve(A)%*%B;
      if (out_type=="normal"){mu<-X1%*%beta.new;Y.res<-Y-mu}
      if (out_type=="beta"){
        mu<-exp(X1%*%beta.new)/(1+exp(X1%*%beta.new));
        Y.res<-Y-mu
      }
      #calculate step beta error
      step.beta<-sum((beta.new-beta)^2);beta<-beta.new
      beta.new
      step.beta
    }
    #Generate RF psuedo variables
    x.field<-S.subject%*%Y.res;
    eta.new<-lm(Y.res~0+x.field)$coefficients
    #calculate step error
    step.error<-sum(c(beta.new-beta,eta.new-eta)^2);
    if (iter >= 49) {
      warning(paste("Warning: Model cannot converge for type ", out_type, " and correlation ", correlation, sep = '')); break
    } else {iter<-iter+1};
    beta<-beta.new;eta<-eta.new
    V<-diag(1,n)-eta*S.subject
  }
  return(list(Y.res=Y.res,distance=distance,X1=X1,eta=eta,beta=beta,delta=delta,V=V))
}

Multi.MRF<-function(Z,result.null,similarity='GR', weights=1, resample=5000)
{
  if (resample!=0){
    p.value<-test.Multi.MRF.B(Z,result.null,similarity='GR',weights=weights,resample=resample);
  }else{
    p.value<-test.FGRF(Z,result.null,similarity='GR')  
  }
  p.value
}

test.Multi.MRF<-function(Z,result.null,similarity='GR') #Z: n*t matrix, result.null: result of null.Multi.MRF
{
  # recall data
  Y.res<-result.null$Y.res;
  distance<-result.null$distance
  X1<-result.null$X1;
  eta<-result.null$eta
  delta<-result.null$delta;
  V<-result.null$V;
  # tests
  # generate IBS pseudo variable
  if (similarity=='IBS'){G<-IBS_pseudo(Z)}
  # centered genotype
  if (similarity=='GR'){G<-apply(Z,2,center)}
  # preliminary
  dimY<-nrow(Y.res);
  m<-length(unique(distance[,1])); 
  dimX<-ncol(X1); 
  dimG<-ncol(G) 
  index.cluster<-unique(distance[,1])
  # generalized score test
  dS<-matrix(0,2*dimG+dimX,dimX); 
  S1<-matrix(0,2*dimG,1); 
  S0<-matrix(0,dimX,1); 
  D<-matrix(0,2*dimG+dimX,2*dimG+dimX)
  S.mean<-t(cbind(t(Y.res)%*%V%*%G,t(Y.res)%*%G,t(delta^(-1/2)*Y.res)%*%V%*%(delta^{1/2}*X1)))/m 
  for (j in index.cluster)
  {
    index<-which(distance[,1]==j);
    ni<-length(index)
    Vi<-V[index,index] 
    Gi<-G[index,]; 
    Xi<-X1[index,]; 
    Yi<-as.matrix(Y.res[index,])
    delta_i<-delta[index] 
    if (ni==1) {Gi<-t(Gi);Xi<-t(Xi)}
    G.temp<-cbind(Vi%*%Gi,Gi,(delta_i^{-1/2}*Vi)%*%(delta_i^{1/2}*Xi)) 
    # calculate cluster based statistics
    Si<-t(G.temp)%*%Yi; 
    dSi<--t(G.temp)%*%(delta_i*Xi)
    Di<-(Si-S.mean)%*%t(Si-S.mean)
    # calculate overall statistics
    S1<-S1+Si[1:(2*dimG),]; 
    dS<-dS+dSi; 
    D<-D+Di 
  }
  D<-D/m;S1<-S1/sqrt(m);dS<-dS/sqrt(m)
  # test statistic
  Q<-t(S1[(dimG+1):(2*dimG),])%*%S1[1:dimG,]
  # null distribution
  P<-dS[1:(2*dimG),]%*%solve(dS[(2*dimG+1):(2*dimG+dimX),])
  D11<-D[1:(2*dimG),1:(2*dimG)];D10<-D[1:(2*dimG),(2*dimG+1):(2*dimG+dimX)];D00<-D[(2*dimG+1):(2*dimG+dimX),(2*dimG+1):(2*dimG+dimX)]
  M<-D11-P%*%t(D10)-D10%*%t(P)+P%*%D00%*%t(P)
  M.rotate<-cbind(M[,(dimG+1):(2*dimG)],M[,1:dimG])
  eigen.M<-Re(eigen(M.rotate,symmetric=FALSE,only.value=TRUE)$values)
  # calculate p.value
  p.value<-davies(2*Q, eigen.M, acc=10^(-6))$Qq
  
  return(list(pvalue=p.value))
}



test.Multi.MRF.B<-function(Z,result.null,similarity='GR', weights=1,resample=5000) #Z: n*t matrix, result.null: result of null.Multi.MRF
{
  set.seed(20131018)
  # recall data
  Y.res<-result.null$Y.res;
  distance<-result.null$distance
  X1<-result.null$X1;
  eta<-result.null$eta
  delta<-result.null$delta;
  V<-result.null$V;
  W <- diag(sqrt(weights),ncol(Z),ncol(Z))
  # tests
  # generate IBS pseudo variable
  if (similarity=='IBS'){G<-IBS_pseudo(Z); G <- G%*%W}
  # centered genotype
  if (similarity=='GR'){G<-apply(Z,2,center); G <- G%*%W}
  # preliminary
  dimY<-nrow(Y.res); 
  m<-length(unique(distance[,1])); 
  dimX<-ncol(X1); 
  dimG<-ncol(G) 
  index.cluster<-unique(distance[,1])
  # generalized score test
  dS<-matrix(0,2*dimG+dimX,dimX); 
  S1<-matrix(0,2*dimG,1); 
  Si.A<-matrix(0,2*dimG,m) 
  S0<-matrix(0,dimX,1);
  D<-matrix(0,2*dimG+dimX,2*dimG+dimX) 
  S.mean<-t(cbind(t(Y.res)%*%V%*%G,t(Y.res)%*%G,t(delta^(-1/2)*Y.res)%*%V%*%(delta^{1/2}*X1)))/m
  for (j in 1:length(index.cluster))
  {
    subj <- index.cluster[j]
    index<-which(distance[,1]==subj);
    ni<-length(index)
    Vi<-V[index,index] 
    Gi<-G[index,]; 
    Xi<-X1[index,]; 
    Yi<-as.matrix(Y.res[index,]) 
    delta_i<-delta[index] 
    if (ni==1) {Gi<-t(Gi);Xi<-t(Xi)}
    G.temp<-cbind(Vi%*%Gi,Gi,(delta_i^{-1/2}*Vi)%*%(delta_i^{1/2}*Xi))
    # calculate cluster based statistics
    S.temp <- t(G.temp)%*%Yi/sqrt(m)
    Si.A[,j]<- S.temp[1:(2*dimG),] 
  }
  S1<-as.matrix(apply(Si.A,1,sum)) 
  # test statistic
  Q<-t(S1[(dimG+1):(2*dimG),])%*%S1[1:dimG,]
  # null distribution
  perturb.coef<-matrix(rbinom(m*resample,1,0.5)*2-1,m,resample)
  perturb.S1<-Si.A%*%perturb.coef
  # calculate bootstrap based mean, variance and kurtosis
  perturb.Q<-NULL
  for (k in 1:resample){perturb.Q<-c(perturb.Q,t(perturb.S1[(dimG+1):(2*dimG),k])%*%perturb.S1[1:dimG,k])}
  perturb.mean<-mean(perturb.Q)
  perturb.variance<-var(perturb.Q)
  perturb.kurtosis<-mean((perturb.Q-perturb.mean)^4)/perturb.variance^2-3
  # calculate p.value
  if (perturb.kurtosis>0){df<-12/perturb.kurtosis}else{
    df<-100000
  }
  p.value<-pchisq((Q-perturb.mean)*sqrt(2*df)/sqrt(perturb.variance)+df,df,lower.tail = F)
  return(list(pvalue=p.value))
}

# generate IBS pseudo variable
IBS_pseudo<-function(x)
{
  pseudo1<-(0<=x & x<0.5)*sqrt(2)+(0.5<=x & x<1.5)*sqrt(2)/2
  pseudo2<-(0.5<=x & x<1.5)*sqrt(2)/2+(1.5<=x & x<2.5)*sqrt(2)
  pseudo3<-(0.5<=x & x<1.5)
  pseudo<-cbind(pseudo1,pseudo2,pseudo3)
  return(pseudo)
}

call_Kernel_IBS<-function(Z,n,p)
{
  #Kernel_IBS(double * Z, int * pn, int * pp, double * Kernel)
  K<- matrix(rep(0,n*n),nrow = n, ncol = n)  
  temp<-.C("Kernel_IBS",as.integer(as.vector(t(Z))),as.integer(n), as.integer(p),as.double(as.vector(K)))[[4]]
  matrix(temp,nrow=n)
}

Get.sqrt<-function(A){
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){
    stop("Error to obtain matrix square!")
  }
  a.sqrt <- a.eig$vectors[,ID1] %*% diag(sqrt(a.eig$values[ID1])) %*% t(a.eig$vectors[,ID1])
  return(a.sqrt)
}

#distance similarity
Check.same<-function(x,y)
{return(as.numeric(x==y))}

Check.near<-function(x,y)
{return(as.numeric(abs(x-y)==1))}

#Centering
center<-function(x)
{return(x-mean(x))}
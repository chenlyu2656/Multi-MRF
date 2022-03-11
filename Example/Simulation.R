library(MASS)
# input dataset: GT (genotype data); info (information for genotype)
# input parameter: size (sample size); causal (proportion of causal mQTLs variants); correlation structures (i.e. exchangeable, autoregressive and distance-based); causal structures between mQTLs and methylation traits (i.e. “unique”, “half-shared”, and “all-shared”)

# example
GT <- read.table("Genotypes.txt", header=T,sep='\t',stringsAsFactors=F)
info <- read.table("info.info", header=T,sep='\t',stringsAsFactors=F)
size<-100; causal<-0.1;causal_structure <- "unique"; correlation="exchangeable";
set.seed(20131018);

row_n <- sample(nrow(GT),size = size, replace = F)
GT_n <- GT[row_n,]
MAF_n <- colMeans(GT_n)/2
MAF_n <- MAF_n[MAF_n != 0]
snpname<-rownames(as.data.frame(MAF_n))

GT_n <- GT_n[,snpname]

# distance dataset
# Assuming 10 CpG sites for each ID, each CpG corrresponds to a distance
# To make sure the CpGs are within the same gene region with genotype, distance range from min(info$posi) to max(info$posi)
# each time, random choose 10kb as the gene region
# within 10kb, random selected 10 CpG sites
distance <- as.data.frame(matrix(rep(rownames(GT_n),each=10),byrow=T,ncol=1))
colnames(distance) <- "ID"
distance$distance <- 0
s <- floor(runif(1, min=min(info$posi), max=max((info$posi))))
e <-s+10000
cpg <- floor(runif(10, min=s, max=e))
for (j in 1:size){
  distance[(1+10*(j-1)):(10*j),2] <- cpg
}

# updated GT_n within the region
snpname2 <- info$SNP[info$posi >= s & (info$posi <= e)]
snpname2 <- snpname2[snpname2 %in% snpname]
GT_n <- GT_n[,snpname2]
MAF_n <- colMeans(GT_n)/2

# generate variance-covariance matrix based on correlation structure
if (correlation == "exchangeable"){vcov<-matrix(0.7,10,10); diag(vcov) <-1}

if (correlation == "autoregressive"){
  vcov<-matrix(1,10,10);
  for (k in 1:10){
    for (j in 1:10){
      di <- min(k,j)
      vcov[k,j] <- 0.7^(k-di)*0.7^(j-di)
    }}}

if (correlation == "exponential"){
  vcov<-matrix(1,10,10);
  dis <- distance[1:10,2]
  td <- 10000 # 10kb as region, can use average distance
  for (k in 1:10){
    for (j in 1:10){
      vcov[k,j] <- exp(-abs(dis[k]-dis[j])/td*3)
    }
  }
}

# generate multivariate error based on vcov matrix
error <- mvrnorm(n=size, rep(0,10), vcov)
error <- error/(10*max(error))

EFF<- -1/(MAF_n*(1-MAF_n))*0.01;
EFF1<-EFF;
EFF2<-EFF;

# wt: weights for snps
wt<-dbeta(MAF_n,1,25)
loc<-1:length(MAF_n);

nloc1<-floor(ncol(GT_n)*causal);
nloc2<-floor(nloc1*0.5);

# select causal snps for 10 traits
# varying causal structure
if (causal_structure == "unique"){
  sn1_10<-sample(loc,size=nloc1*10,replace=F);
  sn2_10<-sample(sn1_10,size=nloc2*10,replace=F);
}

if (causal_structure == "all shared"){
  sn1_10<-rep(sample(loc,size=nloc1,replace=F),10);
  sn2_10<-rep(sample(sn1_10,size=nloc2,replace=F),10);
}

if (causal_structure == "half shared"){
  sn1_5<-sample(loc,size=nloc1*6,replace=F);
  shared1<-sample(sn1_5,size=nloc1,replace=F);
  sn1_10<-c(rep(shared1,5),sn1_5[! sn1_5 %in% shared1]);
  sn2_5<-sample(sn1_10,size=nloc2*6,replace=F);
  shared2<-sample(1:length(sn2_5),size=nloc2,replace=F);
  sn2_10<-c(rep(shared2,5),sn2_5[-shared2]);
}

# generate methylation trait
# y0: no causal mQTLs
# y.l0: y0 in logit transform
# y1: in association with 10% causal mQTLs in one direction
# y.l1: y1 in logit transform
# y2: in association with 10% causal mQTLs in bi-direction
# y.l2: y2 in logit transform
y0 <- as.data.frame(matrix(rep(0,each=10*size),byrow=T,ncol=1));colnames(y0) <- "Y";
y1 <- y0;y2 <- y0;y.l0 <- y0; y.l1 <- y0; y.l2 <- y0;

for (j in 1:10){
  sn1<-sn1_10[(1+nloc1*(j-1)):(nloc1*j)];
  sn2<-sn2_10[(1+nloc2*(j-1)):(nloc2*j)];
  
  effect1<-EFF1;
  effect1[-sn1]<-0;
  effect2<-effect1;
  effect2[sn2]<-effect2[sn2]*(-1);
  
  N <- size
  
  eta0 <- log(0.1/(1-0.1))+error[,j];
  min0 <- min(eta0)
  max0 <- max(eta0)
  eta0 <- eta0*(2.19-2)/(max0-min0)+(-2*min0+2.19*max0)/(min0-max0)
  eta0[eta0 > 2.19] <- 2.19
  
  # eta1 for one direction
  # eta1 = beta1*x1+beta2*x2+...+betak*xk, make min = -2.19, max = c, c to control power
  c <- 1.5
  eta1<-rowSums(GT_n*matrix(effect1,nrow=nrow(GT_n),ncol=ncol(GT_n),byrow=T))+error[,j]
  min1 <- min(eta1)
  max1 <- max(eta1)
  eta1 <- eta1*(2.19+c)/(max1-min1)+(c*min1+2.19*max1)/(min1-max1)
  eta1[eta1 > 2.19] <- 2.19
  
  # eta2 for bi-direction
  # eta2 = beta1*x1+beta2*x2+...+betak*xk, make min = -2.19, max = d, d to control power
  d <- 2.5
  eta2<-rowSums(GT_n*matrix(effect2,nrow=nrow(GT_n),ncol=ncol(GT_n),byrow=T))+error[,j]
  min2 <- min(eta2)
  max2 <- max(eta2)
  eta2 <- eta2*(2.19+d)/(max2-min2)+(d*min2+2.19*max2)/(min2-max2)
  eta2[eta2 > 2.19] <- 2.19
  
  mu0<-exp(eta0)/(1+exp(eta0));
  mu1<-exp(eta1)/(1+exp(eta1));
  mu2<-exp(eta2)/(1+exp(eta2));
  
  phi <- 30
  
  a0 <- mu0*phi
  a1 <- mu1*phi
  a2 <- mu2*phi
  
  b0 <- (1-mu0)*phi
  b1 <- (1-mu1)*phi
  b2 <- (1-mu2)*phi
  
  tmp.y0<-rbeta(N,a0,b0);
  tmp.y1<-rbeta(N,a1,b1);
  tmp.y2<-rbeta(N,a2,b2);
  
  tmp.y.l0<-log(tmp.y0/(1-tmp.y0));
  tmp.y.l1<-log(tmp.y1/(1-tmp.y1));
  tmp.y.l2<-log(tmp.y2/(1-tmp.y2));
  
  y0[seq(j, nrow(y0), 10),1] <- tmp.y0;
  y.l0[seq(j, nrow(y0), 10),1] <- tmp.y.l0;
  y1[seq(j, nrow(y0), 10),1] <- tmp.y1;
  y.l1[seq(j, nrow(y0), 10),1] <- tmp.y.l1;
  y2[seq(j, nrow(y0), 10),1] <- tmp.y2;
  y.l2[seq(j, nrow(y0), 10),1] <- tmp.y.l2;

}
    
# covariates
# each individual has 10 repeated measures
x1<-rnorm(nrow(GT_n))
x2<-sample(c(0,1),size=nrow(GT_n),replace = T,prob=c(0.5,0.5));

# output
# methylation trait data
Y <- cbind(distance$ID, y0, y1,y2,y.l0,y.l1,y.l2)
colnames(Y) <- c("ID","y0","y1","y2","y.l0","y.l1","y.l2")
write.table(Y,"Simulated Methylation Traits.txt",row.names=F,col.names=T,quote=F,sep='\t')

# distance data
write.table(distance,"Simulated Distance.txt",row.names=F,col.names=T,quote=F,sep='\t')

# covariates data
X <- as.data.frame(cbind(x1,x2))
rownames(X) <- rownames(GT_n)
write.table(X,"Simulated Covariates.txt",row.names=T,col.names=T,quote=F,sep='\t')

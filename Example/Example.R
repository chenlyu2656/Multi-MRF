library(CompQuadForm)
library(betareg)
library(MASS)
library(glasso)
library(PearsonDS)

source('Multi-MRF.R')

# Trait data (n=100, nCpG=10)
# Long format
# y0, y1 and y2 represent methylation traits (beta values), ranging from [0, 1] with each following a beta distribution
# y0 was simulated with no causal association with genotypes
# y1 was simulated to causally associated with 10% of SNPs in one direction (positive associations)
# y2 was simulated to causally associated with 10% of SNPs in two directions (50% positive associations and 50% negative associations)
trait <-read.table("Traits.txt", header=T)

# Distance data (n=100, nCpG=10)
# Long format
# distance column represents the BP information for each CpG site
distance <- read.table("Distance.txt", header=T)

# Genotype data (n=100, nsnp=56)
# Mixture of common and rare variants
geno <- read.table("Genotypes.txt", row.names = 1, header=T)

# Covariate data
# X1 is a quantitative covaiate, while X2 is a binary covariate
cov <- read.table("Covariates.txt", row.names = 1, header=T)

# Generate weights
# To detect rare variants of relatively large effect, higher weights are given to rare variants
MAF <- colMeans(geno)/2
wt<-dbeta(MAF,1,25)

# Align with the long format of traits data
GT <- apply(geno, 2, rep, each=10)
COV <- apply(cov, 2, rep, each=10)

# Assuming methylation traits follow a beta distribution and the correlation structure is distance-based (exchangeable correlation can also be applied)
obj.dist.b0<-null.Multi.MRF(Y=trait$y0, distance=distance,X=COV,correlation='distance-based', out_type='beta');
obj.dist.b1<-null.Multi.MRF(Y=trait$y1, distance=distance,X=COV,correlation='distance-based', out_type='beta');
obj.dist.b2<-null.Multi.MRF(Y=trait$y2, distance=distance,X=COV,correlation='distance-based', out_type='beta');

p.dgrf.dist.b0 <-Multi.MRF(Z=GT, obj.dist.b0, weights=(wt^2))$pvalue
p.dgrf.dist.b1 <-Multi.MRF(Z=GT, obj.dist.b1, weights=(wt^2))$pvalue
p.dgrf.dist.b2 <-Multi.MRF(Z=GT, obj.dist.b2, weights=(wt^2))$pvalue

# Assuming methylation trait follows a normal distribution after logit transformation and the correlation is distance-based (exchangeable correlation can also be applied)
y.l0<-log(trait$y0/(1-trait$y0))
y.l1<-log(trait$y1/(1-trait$y1))
y.l2<-log(trait$y2/(1-trait$y2))

obj.dist.l0<-null.Multi.MRF(Y=y.l0, distance=distance,X=COV,correlation='distance-based', out_type='normal');
obj.dist.l1<-null.Multi.MRF(Y=y.l1, distance=distance,X=COV,correlation='distance-based', out_type='normal');
obj.dist.l2<-null.Multi.MRF(Y=y.l2, distance=distance,X=COV,correlation='distance-based', out_type='normal');

p.dgrf.dist.l0 <-Multi.MRF(Z=GT, obj.dist.l0, weights=(wt^2))$pvalue
p.dgrf.dist.l1 <-Multi.MRF(Z=GT, obj.dist.l1, weights=(wt^2))$pvalue
p.dgrf.dist.l2 <-Multi.MRF(Z=GT, obj.dist.l2, weights=(wt^2))$pvalue

# Multi-trait Methylation Random Field Method
The multi-trait methylation random field (multi-MRF) method is developed to detect methylation quantitative trait loci (mQTLs) by evaluating the joint association between a set of CpG sites and a set of genetic variants within a genomic region. 

The proposed method has several advantages:
1) It is a multi-trait method that allows flexible correlation structures between neighboring CpG sites (e.g. distance-based correlation). 
2) It is also a multi-locus method that integrates the effect of multiple common and rare genetic variants. 
3) It models the methylation traits with a beta distribution to characterize their bimodal and interval properties. 

# Multi-MRF Function
Available [here](./R/Multi-MRF.R)

# Tutorial
Sample data availability
- [Methylation traits](./Example/Data/Traits.txt)
- [Distance info for each CpG site](./Example/Data/Distance.txt)
- [Genotypes](./Example/Data/Genotypes.txt)
- [Information for genotypes](./Example/Data/info.info)
- [Covariates](./Example/Data/Covariates.txt)

Sample code availability
- [Sample code](./Example/Example.R)

Simulation code availability
- [Simulation code](./Example/Simulation.R)

The detailed user's guide is [available](./Multi-MRF%20User's%20Guide.pdf)

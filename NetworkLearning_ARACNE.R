# sample code for performing the network analysis.
# ARACNE method

# Step 1. load normalized expression data.

dat<-read.table('TFandModule.csv',sep=',',as.is=T)
dim(dat)
# [1] 90 20

# Separate regulators and target genes.
regulators<-dat[1:30,]
targets<-dat[31:90,]

# To perform CLR analysis, load package parmigene
library(parmigene)

# calculate all pairwise mutual information.
allmi<-knnmi.all(dat)

# perform CLR analysis
aracne=aracne.a(allmi)

# function to make resampled matrix
randmat<-function(inmat)
{
  tmpmat<-as.vector(inmat)
  tmpmat<-sample(tmpmat)
  out<-matrix(tmpmat,ncol=ncol(inmat),nrow=nrow(inmat))
}

getpermat_ara<-function(inmat,nperm=1000)
{ # input parameters: matrix of expression
  # number of permutations
  # output: vectors of permutated average values
  inmat<-as.matrix(inmat)
  out<-matrix(0,ncol=nperm,nrow=1)
  rownames(out)<-c('aracne')
  for (i in 1:ncol(out)){
    tmp1<-randmat(inmat)
    tmp1mi<-knnmi.all(tmp1)
    tmp1mi_aracne<-aracne.a(tmp1mi)
    tmp1_aracne<-tmp1mi_aracne[upper.tri(tmp1mi_aracne)]
    out['aracne',i]<-mean(tmp1_aracne)
  }
  return(out)
}

# make permutation 
datperm<-getpermat_ara(dat,1000)

# get mean and standard deviation
datpermm<-apply(datperm,1,mean)
datpermsd<-apply(datperm,1,sd)

# calculate p values
require(reshape)
netmat<-aracne 

# convert adjancy matrix to edgelist 
elist<-melt(netmat[1:30,30:90])
elist[,1]<-as.character(elist[,1])
elist[,2]<-as.character(elist[,2])
colnames(elist)<-c('regulator','target','score')

#estimate p value with gaussian approximation.
pval<-pnorm(elist[,3],mean=datpermm,sd=datpermsd, 
            lower.tail=FALSE, log.p=TRUE)

# adjust p value with BH method. (padj function does not work with very small p values)
padj<-pval+log(length(pval))-log(rank(pval))

# prepare final output matrix
out_elist<-data.frame(elist,
                      rank=rank(-elist[,'score']),# this give the highest score rank 1
                      pval=exp(pval),
                      padj=exp(padj),
                      stringsAsFactors = FALSE)
out_elist<-out_elist[order(out_elist[,'padj']),]
thre<-0.01
elistsig<-out_elist[out_elist[,'padj']<thre,]

write.table(elistsig,'aracne_network_edges.csv',sep=',',row.names=FALSE)
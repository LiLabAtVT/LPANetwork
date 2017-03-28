# sample code for performing the network analysis.
# CLR method

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
clr=clr(allmi)

# function to make resampled matrix
randmat<-function(inmat)
{
  tmpmat<-as.vector(inmat)
  tmpmat<-sample(tmpmat)
  out<-matrix(tmpmat,ncol=ncol(inmat),nrow=nrow(inmat))
}

getpermat_clr<-function(inmat,nperm=1000)
{ # input parameters: matrix of expression
  # number of permutations
  # output: vectors of permutated average values
  inmat<-as.matrix(inmat)
  out<-matrix(0,ncol=nperm,nrow=1)
  rownames(out)<-c('clr')
  for (i in 1:ncol(out)){
    tmp1<-randmat(inmat)
    tmp1mi<-knnmi.all(tmp1)
    tmp1mi_clr<-clr(tmp1mi)
    tmp1_clr<-tmp1mi_clr[upper.tri(tmp1mi_clr)]
    out['clr',i]<-mean(tmp1_clr)
  }
  return(out)
}

# make permutation 
datperm<-getpermat_clr(dat,1000)

# get mean and standard deviation
datpermm<-apply(datperm,1,mean)
datpermsd<-apply(datperm,1,sd)

# calculate p values
require(reshape)
netmat<-clr 

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

write.table(elistsig,'CLR_network_edges.csv',sep=',',row.names=FALSE)

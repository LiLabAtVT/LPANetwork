# sample code for performing the network analysis.
# Random Forest method

# Step 1. load normalized expression data.

dat<-read.table('TFandModule.csv',sep=',',as.is=T)
dim(dat)
# [1] 90 20

# Separate regulators and target genes.
regulators<-dat[1:30,]
targets<-dat[31:90,]

source('genie3.R')
genienet <- get.weight.matrix(dat, input.idx=1:30,nb.trees = 100, trace=FALSE)

# make random matrix
randmat<-function(inmat)
{
  tmpmat<-as.vector(as.matrix(inmat))
  tmpmat<-sample(tmpmat)
  out<-matrix(tmpmat,ncol=ncol(inmat),nrow=nrow(inmat))
}

# make simulation of genie score.
geniesim<-rep(0,1000)
for(i in 1:1000)
{
  tmpmat<-randmat(dat)
  rownames(tmpmat)<-rownames(dat)
  colnames(tmpmat)<-colnames(dat)
  tmpnet<-get.weight.matrix(tmpmat,input.idx=1:30,nb.trees=100,trace=FALSE)
  geniesim[i]<-mean(tmpnet)
  if(i%%100==0){print(i)}
}

elist<-melt(genienet[1:30,30:90])
elist1<-elist[elist[,'value']>0,]
pval<-pnorm(elist1[,3],mean=mean(geniesim),sd=sd(geniesim), 
               lower.tail=FALSE, log.p=TRUE)
padj<-pval+log(length(pval))-log(rank(pval))
out_elist<-data.frame(elist1,
                      rank=rank(-elist1[,'value']),# this give the highest score rank 1
                      pval=exp(pval),
                      padj=exp(padj),
                      stringsAsFactors = FALSE)
out_elist1<-out_elist[out_elist[,'padj']<0.01,]

write.table(out_elist1,'RF_network_edges.csv',sep=',',row.names=FALSE)


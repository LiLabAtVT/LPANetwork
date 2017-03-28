# sample code for performing the network analysis.
# LARS method

# Step 1. load normalized expression data.
library(lars)
dat<-read.table('TFandModule.csv',sep=',',as.is=T)
dim(dat)
# [1] 90 20

# Separate regulators and target genes.
regulators<-dat[1:30,]
targets<-dat[31:90,]


# Make LARS prediction.
output<-NULL
tarm<-targets
regm<-regulators

ntar<-nrow(tarm) # number of targets
i=1
nreg<-ncol(dat)
nrank<-nreg-1
for(i in 1:ntar){
  y<-t(tarm[i,])
  x<-t(as.matrix(regm[,]))
  fit<-lars(x, y, type='lar',use.Gram=FALSE,max.steps = 100)
  gnames<-names(unlist(fit$actions))
  #print(gnames)
  output<-rbind(output, 
                data.frame(reg=gnames,tar=i,
                           score=as.numeric(coef(fit)[nreg,gnames]),# no score.
                           rank=c(1:nrank),
                           stringsAsFactors = FALSE)
  )
}
output1<-output
output1[,'rank']<-rank(output[,'rank'],ties.method = 'min')
write.table(output1,'LARS_network_edges.csv',sep=',',row.names=FALSE)



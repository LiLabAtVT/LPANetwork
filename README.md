## Welcome to Soybean LPA Gene Network

We applied five machine learning methods to learn gene regulatory networks from RNAseq data generated during soybean seed development. 

These machine learning methods include:

Mutual Information (ARACNE)
Random Forest (RF)
Context Likelihood Relatedness (CLR)
Least Angle Regression (LARS)
Partial Correlation (PCOR)

![Machine Learning Methods](/images/Figure1_website.png)

#### Usage
##### Input data
Input data include gene expression from a time series experiment or multiple tissue types. Data should be in matrix format. With each row represents the expression pattern of a transcription factor (TF), or a gene module. Gene module can be obtained using clustering method. 

##### Perform network analysis
Name the input data as TFandModule.csv, run the R code using Rstudio or Rscript command.

##### Output 
![Gene Regulatory Network](/images/NetworkFigureWebsite.png)

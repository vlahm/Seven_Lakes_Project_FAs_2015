setwd("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/PCA_arrangements")

#this part is not in use yet
FA_PCAer<-function(csv.name)
{  
  
  
  
  #read in the FAs and samples you want to test
  csv<-read.csv("1_allConsumers_FAoverPointFive.csv")
  
  #find number and indices of columns containing Areal data
  Nsamp<-length(grep("Area",colnames(csv)))
  samp.ind<-grep("Area",colnames(csv))
  proportional.area<-matrix(data=NA, nrow=length(csv[,1]), ncol=Nsamp)
  
  #create matrix of proportional area
  for(i in 1:Nsamp)
  {
    col<-csv[,samp.ind[i]]
    NAs<-which(is.na(col))
    col[NAs]<-0.0000001
    total.area<-sum(col)
    for(j in 1:length(col))
    {
      proportional.area[j,i]<-(col[j]*100)/total.area #find out why i and j are being created as objects
    }
  }
  proportional.area<-cbind(csv[,1:5],proportional.area)
  colnames(proportional.area)[6:length(proportional.area[1,])]<-gsub("Reten_", "", colnames(csv)[grep("Reten", colnames(csv))])
  
  #the PCA...
  pca<-prcomp(proportional.area[,c(1,6:length(proportional.area[1,]))],scale=T, scores=T)
  ordi.monte
  #Hmm.  It would seem these files are private.  Waiting to hear from Jen Lang.


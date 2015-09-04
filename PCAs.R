# setwd("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/PCA_arrangements")

####functions and packages####
library(vioplot)
library(vegan)

#FISH560 functions
pca.eigenval <-
  function(x.pca,dim=length(x.pca$sdev),digits=7){
    
    #check for dim limit
    if(dim>length(x.pca$sdev)){
      cat("Only",length(x.pca$sdev),"axes available\n")
      dim<-length(x.pca$sdev)
    }
    
    #calculate some variables
    names<-colnames(x.pca$rotation[,1:dim])
    var<-x.pca$sdev^2
    trace<-sum(var)
    prop.var<-var/trace
    
    #broken-stick distribution
    p<-length(x.pca$sdev)
    y<-rep(0,p)
    for(i in 1:p) y[i]<-1/i
    for(i in 1:p) y[i]<-sum(y[i:p])
    y<-y[1:dim]
    
    #print results    
    cat('Importance of components:\n')
    z<-rbind('Variance(eigenvalue)'=var[1:dim],
             'Proportion of Variance'=prop.var[1:dim],
             'Cumulative Proportion'=cumsum(prop.var[1:dim]),
             'Broken-stick value'=y)
    colnames(z)<-names
    z<-round(z,digits=digits)
    return(z)
  }

pca.eigenvec <-
  function(x.pca,dim=length(x.pca$sdev),
           digits=7,cutoff=0){
    
    #check for dim limit  
    if(dim>ncol(x.pca$rotation)){
      cat("Only",ncol(x.pca$rotation),"axes available\n")
      dim<-ncol(x.pca$rotation)
    }
    
    #print results    
    cat("\nEigenvectors:\n")
    z<-format(round(x.pca$rotation[,1:dim],digits=digits))
    z[abs(x.pca$rotation[,1:dim])<cutoff]<-substring('',1,nchar(z[1,1]))
    z<-as.data.frame(z)
    return(z)
  }

pca.structure <-
  function(x.pca,x,dim=length(x.pca$sdev),
           digits=3,cutoff=0){
    
    #check for dim limit
    if(dim>length(x.pca$sdev)){
      cat("Only",length(x.pca$sdev),"axes available\n")
      dim<-length(x.pca$sdev)
    }
    
    #calculate structure correlations
    z<-cor(x,x.pca$x[,1:dim])
    
    #print results 
    cat("\nStructure Correlations:\n")
    z<-round(z,digits=digits)
    z[abs(z)<cutoff]<-substring('',1,nchar(z[1,1]))
    z<-as.data.frame(z)
    return(z)
  }

ordi.monte <-
  function(x,ord,dim=length(x),perm=1000,center=TRUE,
           scale=TRUE,digits=3,plot=TRUE,col.hist='blue',col.line='red',
           lty=2,las=1,lab=c(5,5,4),...){
    
    p<-length(x)
    if(dim>p){
      cat("Only",p,"axes available\n")
      dim<-p
    }
    
    if(ord=='pca'){
      z<-prcomp(x,center=center,scale=scale) #prcomp analysis
      z.eig<-z$sdev[1:dim]^2 #compute eigenvalues
      z.teig<-t(z.eig) #transpose eigenvalues
      z.teig<-t(matrix(z.teig,length(z.teig),perm)) #fill matrix with eigenvalues
      write('',file='y.csv') #empty outfile if it exists
      for(i in 1:perm){
        y<-apply(x,2,sample) #permute data matrix
        y<-prcomp(y,center=center,scale=scale) #prcomp on permuted matrix
        y<-y$sdev[1:dim]^2 #compute eigenvalues
        y<-as.data.frame(t(y)) #coerce to data frame and transpose
        write.table(y,file='y.csv',sep=',',append=TRUE,row.names=FALSE,col.names=FALSE)
      }
      y<-read.table('y.csv',sep=',',header=FALSE) #read in permutation results
      p.value<-apply(y>z.teig,2,sum) #compute proportion of random distribution smaller than observed
      p.value<-p.value/perm #compute p-value
      names<-colnames(z$rotation[,1:dim]) #add 'PC#' names
    }
    
    else if(ord=='ca'){
      library(vegan)
      z<-cca(x) #correspondence analysis
      z.eig<-z$CA$eig[1:dim] #get eigenvalues
      z.teig<-t(z.eig) #transpose eigenvalues
      z.teig<-t(matrix(z.teig,length(z.teig),perm)) #fill matrix with eigenvalues
      write('',file='y.csv') #empty outfile if it exists
      for(i in 1:perm){
        y<-apply(x,2,sample) #permute data matrix
        y<-cca(y) #CA on permuted matrix
        y<-y$CA$eig[1:dim] #get eigenvalues
        y<-as.data.frame(t(y)) #coerce to data frame and transpose
        write.table(y,file='y.csv',sep=',',append=TRUE,row.names=FALSE,col.names=FALSE)
      }
      y<-read.table('y.csv',sep=',',header=FALSE) #read in permutation results
      p.value<-apply(y>z.teig,2,sum) #compute proportion of random distribution smaller than observed
      p.value<-p.value/perm #compute p-value
      names<-names(z$CA$eig[1:dim]) #add 'CA#' names
    }
    
    else if(ord=='dca'){
      library(vegan)
      if(dim>4){
        cat("Only",4,"axes available\n")
        dim<-4
      }
      z<-decorana(x,...) #detrended correspondence analysis
      z.eig<-z$evals[1:dim] #get eigenvalues
      z.teig<-t(z.eig) #transpose eigenvalues
      z.teig<-t(matrix(z.teig,length(z.teig),perm)) #fill matrix with eigenvalues
      write('',file='y.csv') #empty outfile if it exists
      for(i in 1:perm){
        y<-apply(x,2,sample) #permute data matrix
        y<-decorana(y,...) #DCA on permuted matrix
        y<-y$evals[1:dim] #get eigenvalues
        y<-as.data.frame(t(y)) #coerce to data frame and transpose
        write.table(y,file='y.csv',sep=',',append=TRUE,row.names=FALSE,col.names=FALSE)
      }
      y<-read.table('y.csv',sep=',',header=FALSE) #read in permutation results
      p.value<-apply(y>z.teig,2,sum) #compute proportion of random distribution smaller than observed
      p.value<-p.value/perm #compute p-value
      names<-names(z$eval[1:dim]) #add 'CA#' names
    }
    
    if(plot==TRUE){
      for(i in 1:dim){
        xmin<-min(min(y[[i]],z.eig[i]))
        xmax<-max(max(y[[i]],z.eig[i]))
        hist(y[[i]],col=col.hist,las=las,lab=lab,
             xaxs='i',yaxs='i',xlim=c(xmin,xmax),
             xlab='Eigenvalue',
             main=paste('Random Permutation Distribution of Eigenvalues for',names[i],sep=' '),...)
        abline(v=z.eig[i],col=col.line,lty=lty,lwd=2,...)
        readline("Press return for next plot ")
      }
    }  
    
    cat('Randomization Test of Eigenvalues:\n')
    z<-rbind('Eigenvalue'=z.eig,'P-value'=p.value)
    colnames(z)<-names
    z<-round(z,digits=digits)
    return(z)
  }

####FA_PCAer####

############NOTICE
#this version of FA_PCAer should be set up to calculate percent_area correctly
#(based on the overall total from each sample, rather than just the subsample used in
#the PCA.  However, the spreadsheets have not been set up to run the analysis.  
#This vein of the project was abandoned in favor of FA groupings.  

# FA_PCAer<-function(FA.subset, all.FAs, ordi.monte=F)
# {
#   #read in the FAs and samples you want to test
#   subset<-read.csv(FA.subset)
#   all<-read.csv(all.FAs) #this is necessary for computing proportional areas
#   
#   #find number and indices of columns containing Areal data
#   Nsamp<-length(grep("Area",colnames(subset)))
#   samp.ind<-grep("Area",colnames(subset))
#   proportional.area<-matrix(data=NA, nrow=length(subset[,1]), ncol=Nsamp)
#   
#   #create matrix of proportional area
#   for(i in 1:Nsamp)
#   {
#     subset.col<-subset[,samp.ind[i]]
#     all.col<-all[,samp.ind[i]]
#     total.area<-sum(na.omit(all.col))
#     for(j in 1:length(subset.col))
#     {
#       proportional.area[j,i]<-(subset.col[j]*100)/total.area #find out why i and j are being created as objects
#     }
#     NAs<-which(is.na(proportional.area[,i]))
#     proportional.area[,i][NAs]<-0.0000001
#   }
#   proportional.area<-cbind(subset[,2:5],proportional.area)
#   colnames(proportional.area)[5:length(proportional.area[1,])]<-gsub("Reten_", "", colnames(csv)[grep("Reten", colnames(csv))])
#   row.names(proportional.area)<-subset[,1]
#   
#   #the PCA...
#   transposed<-t(as.matrix(proportional.area[,-(1:4)]))
#   pca<-prcomp(transposed,scale=T, scores=T)
#   #determine eigenvalues
#   eigenvalues<-pca.eigenval(pca)
#   #see which eigenvalues are significant
#   screeplot(pca, bstick=T)
#   #do that a different way
#   if(ordi.monte==T)
#   {
#     ordimonte<-o
#   }
#   #see loadings.  square these to get percentage of variance in each original variable
#   #accounted for by each principal component
#   structure<-pca.structure(pca,transposed,dim=7,cutoff=0.2)
#   #sample.scores<-pca$x[,1:7]
#   #first plotting method
#   biplot(pca)
#   #second plotting method
#   ordiplot(pca,choices=c(1,2), type="text", display="sites")
#   arrows(0,0,pca$rotation[,1]*5, pca$rotation[,2]*5, col="red")
#   text(pca$rotation[,1]*5.2, pca$rotation[,2]*5.2, row.names(pca$rotation))
#   
#   if(ordi.monte==T)
#   {
#     details<-list(eigenvalues[,1:7], ordimonte, structure)
#   }
#   else
#   {
#     details<-list(eigenvalues[,1:7], structure)
#   }
#   return(details)
# }
# 
# FA_PCAer("1_allConsumers_FAoverPointFive.csv", )
# FA_PCAer("2_allConsumers_FAoverPointSeven.csv")
# FA_PCAer("4_caddisonly_29FAs.csv")
# FA_PCAer("5_caddisonly_20FAs.csv")
# FA_PCAer("6_calonly_43FAs.csv")
# FA_PCAer("7_cladonly_28FAs.csv")

####Grouped PCAs - Calanoid####

setwd("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/PCA_arrangements")

raw.calanoid<-read.csv("8_cal_grouped.csv")
calanoid<-t(raw.calanoid[,-(1:2)])
colnames(calanoid)<-raw.calanoid[,1]

cal.pca<-prcomp(calanoid,scale=T, scores=T)
#determine eigenvalues
pca.eigenval(cal.pca)
#see which eigenvalues are significant
screeplot(cal.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
pca.structure(cal.pca,calanoid,dim=7,cutoff=0.2)
#sample.scores<-cal.pca$x[,1:7]
#first plotting method
biplot(cal.pca, main="calanoid")
#second plotting method
ordiplot(cal.pca,choices=c(1,2), type="text", display="sites", main="Calanoid")
arrows(0,0,cal.pca$rotation[,1]*5, cal.pca$rotation[,2]*5, col="blue")
text(cal.pca$rotation[,1]*5.2, cal.pca$rotation[,2]*5.2, row.names(cal.pca$rotation), col="blue")

#comparison of distributions:
vioplot(calanoid[1,], calanoid[2,], calanoid[3,], calanoid[4,], calanoid[5,], 
        calanoid[6,], calanoid[7,], calanoid[8,], calanoid[9,], calanoid[10,], names=(rownames(calanoid)))

#Caddis####

raw.caddis<-read.csv("9_caddis_grouped.csv")
caddis<-t(raw.caddis[,-(1:2)])
colnames(caddis)<-raw.caddis[,1]

caddis.pca<-prcomp(caddis,scale=T, scores=T)
#determine eigenvalues
pca.eigenval(caddis.pca)
#see which eigenvalues are significant
screeplot(caddis.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
pca.structure(caddis.pca,caddis,dim=7,cutoff=0.2)
#sample.scores<-caddis.pca$x[,1:7]
#first plotting method
biplot(caddis.pca, main="caddis")
#second plotting method
ordiplot(caddis.pca,choices=c(1,2), type="text", display="sites", main="Caddisfly", ylim=c(-3,3))
arrows(0,0,caddis.pca$rotation[,1]*4, caddis.pca$rotation[,2]*4, col="blue")
text(caddis.pca$rotation[,1]*4.2, caddis.pca$rotation[,2]*4.2, row.names(caddis.pca$rotation), col="blue")

#comparison of distributions.  Do they seem to differ in the same way that the calanoid ones do?
vioplot(caddis[1,], caddis[2,], caddis[3,], caddis[4,], caddis[5,], 
        caddis[6,], caddis[7,], caddis[8,], caddis[9,], caddis[10,], 
        caddis[11,], names=(rownames(caddis)))

# rows<-rep(1:11, times=length(caddis[1,]))
# cols<-rep(1:7, each=length(caddis[,1]))
# counts<-round(stack(as.data.frame(caddis))[,1]*10000) #can only use integers in contingency tables, so here's
# #a "conversion" to integers that allows some retention of sigdigs.
# ctgcy<-as.data.frame(cbind(counts,rows,cols))
# ctgcy[,2]<-factor(ctgcy[,2])
# ctgcy[,3]<-factor(ctgcy[,3])


#Cladocera####

raw.clado<-read.csv("10_clado_grouped.csv")
clado<-t(raw.clado[,-(1:2)])
colnames(clado)<-raw.clado[,1]

clado.pca<-prcomp(clado,scale=T, scores=T)
#determine eigenvalues
pca.eigenval(clado.pca)
#see which eigenvalues are significant
screeplot(clado.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
pca.structure(clado.pca,clado,dim=7,cutoff=0.2)
#sample.scores<-clado.pca$x[,1:7]
#first plotting method
biplot(clado.pca, main="cladocera")
#second plotting method
ordiplot(clado.pca,choices=c(1,2), type="text", display="sites", main="Cladocera", ylim=c(-3,3))
arrows(0,0,clado.pca$rotation[,1]*5, clado.pca$rotation[,2]*5, col="blue")
text(clado.pca$rotation[,1]*5.2, clado.pca$rotation[,2]*5.2, row.names(clado.pca$rotation), col="blue")

#comparison of distributions
vioplot(clado[1,], clado[2,], clado[3,], clado[4,], clado[5,], 
        clado[6,], names=(rownames(clado)))

#Cladocera and chaoborus####

raw.cc<-read.csv("11_cladANDchaob_grouped.csv")
cc<-t(raw.cc[,-(1:2)])
colnames(cc)<-raw.cc[,1]

cc.pca<-prcomp(cc,scale=T, scores=T)
#determine eigenvalues
pca.eigenval(cc.pca)
#see which eigenvalues are significant
screeplot(cc.pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
pca.structure(cc.pca,cc,dim=7,cutoff=0.2)
#sample.scores<-cc.pca$x[,1:7]
#first plotting method
biplot(cc.pca, main="clad + chaob")
#second plotting method
ordiplot(cc.pca,choices=c(1,2), type="text", display="sites", main="Cladocera + Chaoborus", ylim=c(-3,3))
arrows(0,0,cc.pca$rotation[,1]*5, cc.pca$rotation[,2]*5, col="blue")
text(cc.pca$rotation[,1]*5.2, cc.pca$rotation[,2]*5.2, row.names(cc.pca$rotation), col="blue")

#comparison of distributions
vioplot(cc[1,], cc[2,], cc[3,], cc[4,], cc[5,], 
        cc[6,], cc[7,], cc[8,], names=(rownames(cc)))
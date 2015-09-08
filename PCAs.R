#Michael Vlah - University of Washington - Seven Lakes Project
#vlahm13@gmail.com
#created August 2015
#last edit: 9/8/15

#Notice: this script includes many trial methods that became obsolete.
#They have been retained for future reference and code-borrowing.
#To run the analysis in its final form, you'll need to run the material in 
#part 0 - functions and packages.  Then start with part 8.

# setwd("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/PCA_arrangements")

####0 -functions and packages####
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

#1 - FA_PCAer####

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

#2 - Inefficient Method (1 of 3) Grouped PCAs - probably ignore - Calanoid####

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

#3 - Caddis####

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
#Should run a k-sample K-S test here.  do this for cal and clad too? will implement later.

# rows<-rep(1:11, times=length(caddis[1,]))
# cols<-rep(1:7, each=length(caddis[,1]))
# counts<-round(stack(as.data.frame(caddis))[,1]*10000) #can only use integers in contingency tables, so here's
# #a "conversion" to integers that allows some retention of sigdigs.
# ctgcy<-as.data.frame(cbind(counts,rows,cols))
# ctgcy[,2]<-factor(ctgcy[,2])
# ctgcy[,3]<-factor(ctgcy[,3])

#4 - Cladocera####

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

#5 - Cladocera and chaoborus####

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

#6 - all together####
setwd("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/PCA_arrangements/misc spreadsheets")
all<-read.csv("4a_all_samples_unprocessed.csv")

#store the first five columns (info columns)
info.cols<-all[,1:5]
#remove spacer column
#not necessary yet

#create matrix of relative area of each FA in each sample
area.cols<-grep("Area", colnames(all))
proportional.area<-matrix(data=NA, nrow=length(all[,1]), ncol=length(area.cols))
for(i in 1:length(area.cols))
{
  working.col<-all[,area.cols[i]]
  total.area<-sum(na.omit(working.col))
  for(j in 1:length(working.col))
  {
    proportional.area[j,i]<-(working.col[j]*100)/total.area
  }
}

#replace column names with IDs (ID key here: "C:\Users\Mike\Desktop\Grad\Projects\Thesis\Seven Lakes Project 2014\Data\FA\PCA_arrangements\ID_key.csv")
#first three digits = sample number, next one is lake type, last two are organism type
IDs<-c("i001215", "i002215", "i003103", "i004103", "i006303", "i007203", "i008403", "i009403", "i010203", "011i203", "i012103", "i013120",
       "i090120", "i092303", "i108120", "i094220", "i085220", "i086220", "i087220", "i080317", "i096120", "i081318", "i082320", "i083420", "i084418", 
       "i093118", "i105320", "i106317", "i097118", "i091118", "i095118", "i017111", "i032110", "i029410", "i034210", "i041108", "i044308", "i048308",
       "i051208", "i054408", "i064108", "i071109", "i07300f", "i07400b", "i07500h", "i07600s", "i07700c", "i109107")
colnames(proportional.area)<-IDs

#identify columns with no data / 1 datum
# rowSums(is.na(proportional.area)) == 48
# rowSums(is.na(proportional.area)) == 47

#combine FAs according to biomarker classes
#switching to Excel

write.csv(proportional.area, "C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/PCA_arrangements/misc spreadsheets/4b_converted_to_rel.csv")

#7 - all together - semi-efficient Method (2 of 3)####
#added sorting columns in excel
setwd("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/PCA_arrangements")
# data<-read.csv("ultimate_spreadsheet.csv")
# IDs<-c("i001215", "i002215", "i003103", "i004103", "i006303", "i007203", "i008403", "i009403", "i010203", "011i203", "i012103", "i013120",
#        "i090120", "i092303", "i108120", "i094220", "i085220", "i086220", "i087220", "i080317", "i096120", "i081318", "i082320", "i083420", "i084418", 
#        "i093118", "i105320", "i106317", "i097118", "i091118", "i095118", "i017111", "i032110", "i029410", "i034210", "i041108", "i044308", "i048308",
#        "i051208", "i054408", "i064108", "i071109", "i00000t", "i109107")
# colnames(data)[17:60]<-IDs
# 
# #create groupings
# short_SAFA<-data[data$short_SAFA==1,c(1:5, 17:60)] #isolate rows
#   short_SAFA<-colSums(short_SAFA[,7:49], na.rm=T) #sum columns
# iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, 17:60)]
#   iSAFA_etc<-colSums(iSAFA_etc[,7:49], na.rm=T)
# other_MUFA<-data[data$other_MUFA==1,c(1:5, 17:60)]
#   other_MUFA<-colSums(other_MUFA[,7:49], na.rm=T)
# LIN<-data[data$LIN==1,c(1:5, 17:60)]
#   LIN<-colSums(LIN[,7:49], na.rm=T)
# other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, 17:60)]
#   other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,7:49], na.rm=T)
# OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, 17:60)]
#   OA_16.4n3<-colSums(OA_16.4n3[,7:49], na.rm=T)
# ALA_SDA<-data[data$ALA_SDA==1,c(1:5, 17:60)]
#   ALA_SDA<-colSums(ALA_SDA[,7:49], na.rm=T)
# long_SAFA<-data[data$long_SAFA==1,c(1:5, 17:60)]
#   long_SAFA<-colSums(long_SAFA[,7:49], na.rm=T)
# ARA<-data[data$ARA==1,c(1:5, 17:60)]
#   ARA<-colSums(ARA[,7:49], na.rm=T)
# EPA_DHA<-data[data$EPA_DHA==1,c(1:5, 17:60)]
#   EPA_DHA<-colSums(EPA_DHA[,7:49], na.rm=T)
# other_PUFA<-data[data$other_PUFA==1,c(1:5, 17:60)]
#   other_PUFA<-colSums(other_PUFA[,7:49], na.rm=T)
# 
# grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
#       ALA_SDA, long_SAFA, ARA, EPA_DHA, other_PUFA)
# 
# #need to combine samples for periphyton and POM.  too many zero values in those 
# #categories.  Is it legal to combine, or are the distributions very dissimilar?
# #comparing probability distributions of periphyton
# qqplot(data[,49], data[,50]); abline(0,1)
# qqplot(data[,49], data[,51]); abline(0,1) 
# qqplot(data[,51], data[,50]); abline(0,1)
# #probably not cool to combine periphyton.  Comparing POM.
# qqplot(data[,52], data[,53]); abline(0,1)
# qqplot(data[,52], data[,54]); abline(0,1)
# qqplot(data[,52], data[,55]); abline(0,1) #similar
# qqplot(data[,52], data[,56]); abline(0,1) #similar
# qqplot(data[,52], data[,57]); abline(0,1) #only two values
# qqplot(data[,53], data[,54]); abline(0,1)
# qqplot(data[,53], data[,55]); abline(0,1)
# qqplot(data[,53], data[,56]); abline(0,1)
# qqplot(data[,53], data[,57]); abline(0,1) #only two values
# qqplot(data[,54], data[,55]); abline(0,1)
# qqplot(data[,54], data[,56]); abline(0,1)
# qqplot(data[,54], data[,57]); abline(0,1) #only two values
# qqplot(data[,55], data[,56]); abline(0,1) #similar
# qqplot(data[,55], data[,57]); abline(0,1) #only two values
# qqplot(data[,56], data[,57]); abline(0,1) #only two values
# #probably not cool to combine.  Gonna do it anyway.

#combined POM and peri samples in Excel
data<-read.csv("ultimate_spreadsheet_POMandPeriCombinedLakes.csv")
IDs<-c("i001215", "i002215", "i003103", "i004103", "i006303", "i007203", "i008403", "i009403", "i010203", "011i203", "i012103", "i013120",
       "i090120", "i092303", "i108120", "i094220", "i085220", "i086220", "i087220", "i080317", "i096120", "i081301", "i082320", "i083420", "i084401", 
       "i093101", "i105320", "i106317", "i097101", "i091101", "i095101", "i017111", "i000504", "i000508", "i071109", "i000514", "i109107")
colnames(data)[16:52]<-IDs

#create groupings
short_SAFA<-data[data$short_SAFA==1,c(1:5, 16:52)] #isolate rows
short_SAFA<-colSums(short_SAFA[,6:42], na.rm=T) #sum columns
iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, 16:52)]
iSAFA_etc<-colSums(iSAFA_etc[,6:42], na.rm=T)
other_MUFA<-data[data$other_MUFA==1,c(1:5, 16:52)]
other_MUFA<-colSums(other_MUFA[,6:42], na.rm=T)
LIN<-data[data$LIN==1,c(1:5, 16:52)]
LIN<-colSums(LIN[,6:42], na.rm=T)
other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, 16:52)]
other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,6:42], na.rm=T)
OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, 16:52)]
OA_16.4n3<-colSums(OA_16.4n3[,6:42], na.rm=T)
ALA_SDA<-data[data$ALA_SDA==1,c(1:5, 16:52)]
ALA_SDA<-colSums(ALA_SDA[,6:42], na.rm=T)
long_SAFA<-data[data$long_SAFA==1,c(1:5, 16:52)]
long_SAFA<-colSums(long_SAFA[,6:42], na.rm=T)
EPA_DHA<-data[data$EPA_DHA==1,c(1:5, 16:52)]
EPA_DHA<-colSums(EPA_DHA[,6:42], na.rm=T)
other_PUFA<-data[data$other_PUFA==1,c(1:5, 16:52)]
other_PUFA<-colSums(other_PUFA[,6:42], na.rm=T)

#grouping 1
grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
               ALA_SDA, long_SAFA, EPA_DHA, other_PUFA)
#grouping 2 (removed "other_MUFA" and "Other PUFA"
# grouped<-rbind(short_SAFA, iSAFA_etc, LIN, other_n3_n6_PUFA, OA_16.4n3,
#                ALA_SDA, long_SAFA, EPA_DHA)

#replace zeros with small values
for (i in 1:length(grouped[1,]))
{
  for (j in 1:length(grouped[,1]))
  {
    if (grouped[j,i]==0)
    {
    grouped[j,i]<-0.000001
    }
  }
}

#PCA, ignoring "other_PUFA" (which has many zeros and is not very informative;
#note that it could be combined with "other_n3_n6_PUFA", which also has a lot of zeros.)

tgrouped<-t(as.matrix(grouped))
pca<-prcomp(tgrouped,scale=T, scores=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(pca)
#see which eigenvalues are significant
screeplot(pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(pca,tgrouped,dim=7,cutoff=0.2)
#sample.scores<-pca$x[,1:7]
#first plotting method
biplot(pca)
#second plotting method
fig1<-ordiplot(pca, choices=c(1,2), type="none", xlim=c(-4,4), ylim=c(-4,4))
points(fig1, "sites", pch=as.numeric(substr(row.names(tgrouped),6,7)), col=as.numeric(substr(row.names(tgrouped),5,5)))  #good example of argument matching for later use (match.arg)
arrows(0,0,pca$rotation[,1]*5, pca$rotation[,2]*5, col="black")
text(pca$rotation[,1]*5.8, pca$rotation[,2]*5.8, row.names(pca$rotation))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly", 
                           "fish", "POM", "periphyton", "filamentous", "moss", 
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))

####8 - efficient method (3 of 3)####

#combined POM and peri samples in Excel
data<-read.csv("ultimate_spreadsheet_POMandPeriCombinedLakes.csv")
IDs<-c("i001215", "i002215", "i003103", "i004103", "i006303", "i007203", "i008403", "i009403", "i010203", "011i203", "i012103", "i013120",
       "i090120", "i092303", "i108120", "i094220", "i085220", "i086220", "i087220", "i080317", "i096120", "i081301", "i082320", "i083420", "i084401", 
       "i093101", "i105320", "i106317", "i097101", "i091101", "i095101", "i017111", "i000504", "i000508", "i071109", "i000514", "i109107")
colnames(data)[16:52]<-IDs

cal.cols<-grep("20", substr(colnames(data),6,7))
clad.cols<-grep("01", substr(colnames(data),6,7))
caddis.cols<-grep("03", substr(colnames(data),6,7))
chaob.cols<-grep("17", substr(colnames(data),6,7))
fish.cols<-grep("15", substr(colnames(data),6,7))
POM.cols<-grep("08", substr(colnames(data),6,7))
peri.cols<-grep("04", substr(colnames(data),6,7))
fil.cols<-grep("11", substr(colnames(data),6,7))
moss.cols<-grep("07", substr(colnames(data),6,7))
soil.cols<-grep("09", substr(colnames(data),6,7))
plant.cols<-grep("14", substr(colnames(data),6,7))

####9 - Calanoid####
#create groupings
short_SAFA<-data[data$short_SAFA==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
short_SAFA<-colSums(short_SAFA[,6:22], na.rm=T)
iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
iSAFA_etc<-colSums(iSAFA_etc[,6:22], na.rm=T)
other_MUFA<-data[data$other_MUFA==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_MUFA<-colSums(other_MUFA[,6:22], na.rm=T)
LIN<-data[data$LIN==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
LIN<-colSums(LIN[,6:22], na.rm=T)
other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,6:22], na.rm=T)
OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
OA_16.4n3<-colSums(OA_16.4n3[,6:22], na.rm=T)
ALA_SDA<-data[data$ALA_SDA==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
ALA_SDA<-colSums(ALA_SDA[,6:22], na.rm=T)
long_SAFA<-data[data$long_SAFA==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
long_SAFA<-colSums(long_SAFA[,6:22], na.rm=T)
EPA_DHA<-data[data$EPA_DHA==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
EPA_DHA<-colSums(EPA_DHA[,6:22], na.rm=T)
other_PUFA<-data[data$other_PUFA==1,c(1:5, cal.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_PUFA<-colSums(other_PUFA[,6:22], na.rm=T)

#grouping 1
grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
               ALA_SDA, long_SAFA, EPA_DHA, other_PUFA)
#grouping 2 (removed "other_MUFA" and "Other PUFA"
# grouped<-rbind(short_SAFA, iSAFA_etc, LIN, other_n3_n6_PUFA, OA_16.4n3,
#                ALA_SDA, long_SAFA, EPA_DHA)

#replace zeros with small values
for (i in 1:length(grouped[1,]))
{
  for (j in 1:length(grouped[,1]))
  {
    if (grouped[j,i]==0)
    {
      grouped[j,i]<-0.000001
    }
  }
}

#PCA, ignoring "other_PUFA" (which has many zeros and is not very informative;
#note that it could be combined with "other_n3_n6_PUFA", which also has a lot of zeros.)

tgrouped<-t(as.matrix(grouped))
pca<-prcomp(tgrouped,scale=T, scores=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(pca)
#see which eigenvalues are significant
screeplot(pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(pca,tgrouped,dim=7,cutoff=0.2)
#sample.scores<-pca$x[,1:7]
#first and second principal components
fig1<-ordiplot(pca, choices=c(1,2), type="none", xlim=c(-4,5), ylim=c(-4,4))
points(fig1, "sites", pch=as.numeric(substr(row.names(tgrouped),6,7)), col=as.numeric(substr(row.names(tgrouped),5,5)))
arrows(0,0,pca$rotation[,1]*5, pca$rotation[,2]*5, col="black")
text(pca$rotation[,1]*5.8, pca$rotation[,2]*5.8, row.names(pca$rotation))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly", 
                           "fish", "POM", "periphyton", "filamentous", "moss", 
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))


#plotting 3rd and 4th principal components
fig1<-ordiplot(pca, choices=c(3,4), type="none", xlim=c(-4,5), ylim=c(-4,4))
points(fig1, "sites", pch=as.numeric(substr(row.names(tgrouped),6,7)), col=as.numeric(substr(row.names(tgrouped),5,5)))
arrows(0,0,pca$rotation[,3]*5, pca$rotation[,4]*5, col="black")
text(pca$rotation[,3]*5.8, pca$rotation[,4]*5.8, row.names(pca$rotation))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly", 
                           "fish", "POM", "periphyton", "filamentous", "moss", 
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))

####10 - Caddisfly####
#create groupings
short_SAFA<-data[data$short_SAFA==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
short_SAFA<-colSums(short_SAFA[,6:21], na.rm=T)
iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
iSAFA_etc<-colSums(iSAFA_etc[,6:21], na.rm=T)
other_MUFA<-data[data$other_MUFA==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_MUFA<-colSums(other_MUFA[,6:21], na.rm=T)
LIN<-data[data$LIN==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
LIN<-colSums(LIN[,6:21], na.rm=T)
other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,6:21], na.rm=T)
OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
OA_16.4n3<-colSums(OA_16.4n3[,6:21], na.rm=T)
ALA_SDA<-data[data$ALA_SDA==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
ALA_SDA<-colSums(ALA_SDA[,6:21], na.rm=T)
long_SAFA<-data[data$long_SAFA==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
long_SAFA<-colSums(long_SAFA[,6:21], na.rm=T)
EPA_DHA<-data[data$EPA_DHA==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
EPA_DHA<-colSums(EPA_DHA[,6:21], na.rm=T)
other_PUFA<-data[data$other_PUFA==1,c(1:5, caddis.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_PUFA<-colSums(other_PUFA[,6:21], na.rm=T)

#grouping 1
grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
               ALA_SDA, long_SAFA, EPA_DHA, other_PUFA)
#grouping 2 (removed "other_MUFA" and "Other PUFA"
# grouped<-rbind(short_SAFA, iSAFA_etc, LIN, other_n3_n6_PUFA, OA_16.4n3,
#                ALA_SDA, long_SAFA, EPA_DHA)

#replace zeros with small values
for (i in 1:length(grouped[1,]))
{
  for (j in 1:length(grouped[,1]))
  {
    if (grouped[j,i]==0)
    {
      grouped[j,i]<-0.000001
    }
  }
}

#PCA, ignoring "other_PUFA" (which has many zeros and is not very informative;
#note that it could be combined with "other_n3_n6_PUFA", which also has a lot of zeros.)

tgrouped<-t(as.matrix(grouped))
pca<-prcomp(tgrouped,scale=T, scores=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(pca)
#see which eigenvalues are significant
screeplot(pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(pca,tgrouped,dim=7,cutoff=0.2)
#sample.scores<-pca$x[,1:7]
#first and second principal components
fig1<-ordiplot(pca, choices=c(1,2), type="none", xlim=c(-4,5), ylim=c(-4,4))
points(fig1, "sites", pch=as.numeric(substr(row.names(tgrouped),6,7)), col=as.numeric(substr(row.names(tgrouped),5,5)))
arrows(0,0,pca$rotation[,1]*5, pca$rotation[,2]*5, col="black")
text(pca$rotation[,1]*5.8, pca$rotation[,2]*5.8, row.names(pca$rotation))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly", 
                           "fish", "POM", "periphyton", "filamentous", "moss", 
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))


#plotting 3rd and 4th principal components
fig1<-ordiplot(pca, choices=c(3,4), type="none", xlim=c(-4,5), ylim=c(-4,4))
points(fig1, "sites", pch=as.numeric(substr(row.names(tgrouped),6,7)), col=as.numeric(substr(row.names(tgrouped),5,5)))
arrows(0,0,pca$rotation[,3]*5, pca$rotation[,4]*5, col="black")
text(pca$rotation[,3]*5.8, pca$rotation[,4]*5.8, row.names(pca$rotation))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly", 
                           "fish", "POM", "periphyton", "filamentous", "moss", 
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))

####11 - Cladocera####
#create groupings
short_SAFA<-data[data$short_SAFA==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
short_SAFA<-colSums(short_SAFA[,6:17], na.rm=T)
iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
iSAFA_etc<-colSums(iSAFA_etc[,6:17], na.rm=T)
other_MUFA<-data[data$other_MUFA==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_MUFA<-colSums(other_MUFA[,6:17], na.rm=T)
LIN<-data[data$LIN==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
LIN<-colSums(LIN[,6:17], na.rm=T)
other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,6:17], na.rm=T)
OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
OA_16.4n3<-colSums(OA_16.4n3[,6:17], na.rm=T)
ALA_SDA<-data[data$ALA_SDA==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
ALA_SDA<-colSums(ALA_SDA[,6:17], na.rm=T)
long_SAFA<-data[data$long_SAFA==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
long_SAFA<-colSums(long_SAFA[,6:17], na.rm=T)
EPA_DHA<-data[data$EPA_DHA==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
EPA_DHA<-colSums(EPA_DHA[,6:17], na.rm=T)
other_PUFA<-data[data$other_PUFA==1,c(1:5, clad.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_PUFA<-colSums(other_PUFA[,6:17], na.rm=T)

#grouping 1
grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
               ALA_SDA, long_SAFA, EPA_DHA, other_PUFA)
#grouping 2 (removed "other_MUFA" and "Other PUFA"
# grouped<-rbind(short_SAFA, iSAFA_etc, LIN, other_n3_n6_PUFA, OA_16.4n3,
#                ALA_SDA, long_SAFA, EPA_DHA)

#replace zeros with small values
for (i in 1:length(grouped[1,]))
{
  for (j in 1:length(grouped[,1]))
  {
    if (grouped[j,i]==0)
    {
      grouped[j,i]<-0.000001
    }
  }
}

#PCA, ignoring "other_PUFA" (which has many zeros and is not very informative;
#note that it could be combined with "other_n3_n6_PUFA", which also has a lot of zeros.)

tgrouped<-t(as.matrix(grouped))
pca<-prcomp(tgrouped,scale=T, scores=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(pca)
#see which eigenvalues are significant
screeplot(pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(pca,tgrouped,dim=7,cutoff=0.2)
#sample.scores<-pca$x[,1:7]
#first and second principal components
fig1<-ordiplot(pca, choices=c(1,2), type="none", xlim=c(-4,5), ylim=c(-4,5))
points(fig1, "sites", pch=as.numeric(substr(row.names(tgrouped),6,7)), col=as.numeric(substr(row.names(tgrouped),5,5)))
arrows(0,0,pca$rotation[,1]*5, pca$rotation[,2]*5, col="black")
text(pca$rotation[,1]*5.8, pca$rotation[,2]*5.8, row.names(pca$rotation))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly", 
                           "fish", "POM", "periphyton", "filamentous", "moss", 
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))


#plotting 3rd and 4th principal components
fig1<-ordiplot(pca, choices=c(3,4), type="none", xlim=c(-4,5), ylim=c(-4,4))
points(fig1, "sites", pch=as.numeric(substr(row.names(tgrouped),6,7)), col=as.numeric(substr(row.names(tgrouped),5,5)))
arrows(0,0,pca$rotation[,3]*5, pca$rotation[,4]*5, col="black")
text(pca$rotation[,3]*5.8, pca$rotation[,4]*5.8, row.names(pca$rotation))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly", 
                           "fish", "POM", "periphyton", "filamentous", "moss", 
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))

####12 - All####

short_SAFA<-data[data$short_SAFA==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
short_SAFA<-colSums(short_SAFA[,6:42], na.rm=T)
iSAFA_etc<-data[data$iSAFA_etc==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
iSAFA_etc<-colSums(iSAFA_etc[,6:42], na.rm=T)
other_MUFA<-data[data$other_MUFA==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_MUFA<-colSums(other_MUFA[,6:42], na.rm=T)
LIN<-data[data$LIN==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
LIN<-colSums(LIN[,6:42], na.rm=T)
other_n3_n6_PUFA<-data[data$other_n3_n6_PUFA==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_n3_n6_PUFA<-colSums(other_n3_n6_PUFA[,6:42], na.rm=T)
OA_16.4n3<-data[data$OA_16.4n3==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
OA_16.4n3<-colSums(OA_16.4n3[,6:42], na.rm=T)
ALA_SDA<-data[data$ALA_SDA==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
ALA_SDA<-colSums(ALA_SDA[,6:42], na.rm=T)
long_SAFA<-data[data$long_SAFA==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
long_SAFA<-colSums(long_SAFA[,6:42], na.rm=T)
EPA_DHA<-data[data$EPA_DHA==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
EPA_DHA<-colSums(EPA_DHA[,6:42], na.rm=T)
other_PUFA<-data[data$other_PUFA==1,c(1:5, cal.cols, fish.cols, caddis.cols, clad.cols, chaob.cols, POM.cols, peri.cols, fil.cols, moss.cols, soil.cols, plant.cols)]
other_PUFA<-colSums(other_PUFA[,6:42], na.rm=T)

#grouping 1
grouped<-rbind(short_SAFA, iSAFA_etc, other_MUFA, LIN, other_n3_n6_PUFA, OA_16.4n3,
               ALA_SDA, long_SAFA, EPA_DHA, other_PUFA)
#grouping 2 (removed "other_MUFA" and "Other PUFA"
# grouped<-rbind(short_SAFA, iSAFA_etc, LIN, other_n3_n6_PUFA, OA_16.4n3,
#                ALA_SDA, long_SAFA, EPA_DHA)

#replace zeros with small values
for (i in 1:length(grouped[1,]))
{
  for (j in 1:length(grouped[,1]))
  {
    if (grouped[j,i]==0)
    {
      grouped[j,i]<-0.000001
    }
  }
}

#PCA, ignoring "other_PUFA" (which has many zeros and is not very informative;
#note that it could be combined with "other_n3_n6_PUFA", which also has a lot of zeros.)

tgrouped<-t(as.matrix(grouped))
pca<-prcomp(tgrouped,scale=T, scores=T)
#determine eigenvalues
eigenvalues<-pca.eigenval(pca)
#see which eigenvalues are significant
screeplot(pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
structure<-pca.structure(pca,tgrouped,dim=7,cutoff=0.2)
#sample.scores<-pca$x[,1:7]
#first and second principal components
fig1<-ordiplot(pca, choices=c(1,2), type="none", xlim=c(-4,5), ylim=c(-4,4))
points(fig1, "sites", pch=as.numeric(substr(row.names(tgrouped),6,7)), col=as.numeric(substr(row.names(tgrouped),5,5)))
arrows(0,0,pca$rotation[,1]*5, pca$rotation[,2]*5, col="black")
text(pca$rotation[,1]*5.8, pca$rotation[,2]*5.8, row.names(pca$rotation))
legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly", 
                           "fish", "POM", "periphyton", "filamentous", "moss", 
                           "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
                                                            11, 7, 9, 14))
legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
       fill=c(1, 2, 3, 4, 5))


# #plotting 3rd and 4th principal components
# fig1<-ordiplot(pca, choices=c(3,4), type="none", xlim=c(-4,5), ylim=c(-4,4))
# points(fig1, "sites", pch=as.numeric(substr(row.names(tgrouped),6,7)), col=as.numeric(substr(row.names(tgrouped),5,5)))
# arrows(0,0,pca$rotation[,3]*5, pca$rotation[,4]*5, col="black")
# text(pca$rotation[,3]*5.8, pca$rotation[,4]*5.8, row.names(pca$rotation))
# legend("topleft", legend=c("calanoida", "cladocera", "chaoborus", "caddisfly", 
#                            "fish", "POM", "periphyton", "filamentous", "moss", 
#                            "soil", "terrest. plant"), pch=c(20, 1, 17, 3, 15, 8, 4,
#                                                             11, 7, 9, 14))
# legend("bottomleft", legend=c("puddles", "high", "milk+clear", "big", "combined"),
#        fill=c(1, 2, 3, 4, 5))
# 

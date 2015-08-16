dirs<-c("1_mirror_fish","2_mirror_fish","3_C_caddis","4_C_caddis","6_milk_caddis","7_mirror_caddis","8_morg_caddis","9_morg_caddis","10_noname_caddis","11_y015_caddis","12_y025_caddis","13_Z_caddis","90_C_cal","92_O_cal","108_clear_caddis")

FAmerger<-function(dirs.list)
{
  #list that will contain each data frame
  reports<-list()
  
  #vector X will be used to merge and sort each data frame
  X<-as.data.frame(1:98)
  colnames(X)<-"X"
  
  #loop through directories and perform data filtering/extraction
  for(i in 1:length(dirs.list))
  {
    data<-read.csv(paste("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/GC readouts/New folder (2)/",dirs.list[i],"/REPORT01.CSV",sep=""))
    data<-data[!is.na(data$X),]
    data<-merge(X,data, all.y=T)
    data<-data[,c(1,3,6)]
    colnames(data)<-c("Peak", paste("Reten",strsplit(substr(dirs.list[i],1,3), "[^0-9]+"), sep=""), paste("Area",strsplit(substr(dirs.list[i],1,3), "[^0-9]+"),sep=""))
    reports[[i]]<-data
  }
  
  #bind all reports together
  merged<-Reduce(function (...) {merge(..., all=T)}, reports)
  
  #make vector of total appearances of each FA
  FA.frequency<-vector(length=length(merged[,1]))
  for(i in 1:length(merged[,1]))
  {
    FA.frequency[i]<-sum(sapply(na.omit(as.numeric(as.vector(merged[i,2:31]))),is.numeric))/2
  }
  
  #bind vector to "merged"
  merged[[32]]<-FA.frequency
  colnames(merged)[32]<-"total_appearances"
  
  return(merged)
}

x<-FAmerger(dirs)
write.csv(x, file="C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/consumer_FA_data.csv")

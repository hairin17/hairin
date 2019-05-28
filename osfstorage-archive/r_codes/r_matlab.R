#read matlab data & I made a structure like our "class sample dataset" 
#STUDY1 (n = 28, trials = 80)
rm(list=ls())
setwd("/Users/hairin/Library/Mobile Documents/com~apple~CloudDocs/exp_seminar_2019/HWs/osfstorage-archive/Data")

library(R.matlab)
matlabFile  <- readMat('DATA1.mat')
datList     <- matlabFile$Data.Study1
n_subj <-length(matlabFile$Data.Study1)
x<-c("trial","cue","cuetype","heat","ex_pain","pain_rate","sub_id")
dat<- as.data.frame(datList[[1]][[1]]) 
dat["id"] <- NA
dat$id<-1
colnames(dat)<-x

for (i in 2:n_subj) {
  tmp_dat<-as.data.frame(datList[[i]][[1]]) #sub1 dat
  tmp_dat["id"] <- NA
  tmp_dat$id<-i
  colnames(tmp_dat)<-x
  dat<-rbind(dat,tmp_dat)
}
# Write out data
write.table(dat, file = "study1_data.txt", row.names = F, col.names = T, sep = "\t")


#STUDY2 (n = 33, trials = 70)
matlabFile  <- readMat('DATA2.mat')
datList     <- matlabFile$Data.Study2
n_subj <-length(matlabFile$Data.Study2)
x<-c("trial","cue","cuetype","heat","ex_pain","pain_rate","nps","sub_id")
dat_2<- as.data.frame(datList[[1]][[1]]) 
dat_2["id"] <- NA
dat_2$id<-1
colnames(dat_2)<-x

for (i in 2:n_subj) {
  tmp_dat<-as.data.frame(datList[[i]][[1]]) #sub1 dat
  tmp_dat["id"] <- NA
  tmp_dat$id<-i
  colnames(tmp_dat)<-x
  dat_2<-rbind(dat_2,tmp_dat)
}

# Write out data
write.table(dat_2, file = "study2_data.txt", row.names = F, col.names = T, sep = "\t")

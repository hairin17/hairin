rm(list=ls())

setwd("/Users/hairin/Library/Mobile Documents/com~apple~CloudDocs/exp_seminar_2019/HWs/osfstorage-archive/Data")
data_raw <- read.table('study1_data.txt', header = TRUE, sep = '\t')

list_subj <- unique(data_raw$sub_id)
n_subj <- length(list_subj)
t_subj <- sapply(list_subj, function(id) { sum(data_raw$sub_id == id) })
t_max <- max(t_subj)
nCues = 8

painRating   <- array(0,  c(t_max, n_subj)) #[nTrial,nSubj] 
expectPain   <- array(0,  c(t_max, n_subj))
noxInput   <- array(0,  c(t_max, n_subj))
cueType   <- array(0,  c(t_max, n_subj))

for (i in 1:n_subj) {
  tmp_subj <- data_raw[data_raw$sub_id == list_subj[i], ]
  
  painRating[1:t_subj[i],i]   <- tmp_subj$pain_rate[1:t_subj[i]]
  expectPain[1:t_subj[i],i]   <- tmp_subj$ex_pain[1:t_subj[i]]
  noxInput[1:t_subj[i],i]   <- tmp_subj$heat[1:t_subj[i]]
  cueType[1:t_subj[i],i] <-tmp_subj$cuetype[1:t_subj[i]]
}

#initExpect mean according to nCue (1-8)
initExpect <-array(0, c(nCues,n_subj))
for (s in 1:n_subj) {
  muCue <- array(0, c(nCues,1))
      for (k in 1:nCues) {
      tmpCue <- subset(tmp_subj,cue==k)
      tmpMuCue<-tmpCue$ex_pain
      muCue[k,] <-mean(tmpMuCue)
      }
initExpect[,s]<-muCue
}

###need to be specified
#uniformPriorParameter
#Mymulti?
data_ra <- list(nCues = nCues, nSubject = n_subj, nTrials = t_max, painRating = painRating,
                expectPain = expectPain, noxInput = noxInput, cueType = cueType, myMulti = 1,initExpect = initExpect)


#Stan data structure
#data
{
int<lower=1> nCues; // number of cues
int<lower=1> nSubject; // number of subjects
int<lower=1> nTrials; // number of subjects
matrix[nTrials, nSubject] painRating; // observations
matrix[nTrials, nSubject] expectPain; // observations
matrix[nTrials, nSubject] noxInput; // observations
int<lower=1> myMulti;
int<lower=1> uniformPriorParameter;
int cueIndexVec[(nTrials-nCues), nSubject]; // indexing next trial for spec cue
int firstTimeCue[nCues, nSubject];
matrix[nTrials, nSubject] cueType; // 
  matrix[nCues,  nSubject] initExpect; //  initial expectation mu = expectation rating
}
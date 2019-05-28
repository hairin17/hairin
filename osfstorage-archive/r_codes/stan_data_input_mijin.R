rm(list=ls())

# Reading Data

setwd("/Users/mijinkwon/OneDrive/2019 Spring/courses/computational modeling/final project/compmodeling_finalproj")
setwd("/Users/hairin/Library/Mobile Documents/com~apple~CloudDocs/exp_seminar_2019/HWs/osfstorage-archive") #for hairin

data_raw <- read.table('DATA/study1_data.txt', header = TRUE, sep = '\t')


list_subj <- unique(data_raw$sub_id)
n_subj <- length(list_subj)
t_subj <- sapply(list_subj, function(id) { sum(data_raw$sub_id == id) })
t_max <- max(t_subj)
nCues = 8

# pain rating, expected pain, noxious input, cue type

painRating   <- array(0,  c(t_max, n_subj)) # [nTrial,nSubj] 
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

# initExpect mean according to nCue (1-8)

initExpect <- array(0, c(nCues,n_subj))
muCue <- array(0, c(nCues,1))

for (s in 1:n_subj) {
  tmp_subj <- data_raw[data_raw$sub_id == list_subj[s], ]
  for (k in 1:nCues) {
    tmpCue <- subset(tmp_subj, cue == k)
    tmpMuCue <- tmpCue$ex_pain
    muCue[k, ] <- mean(tmpMuCue)
    initExpect[k, s] <- muCue[k, ]
  }
}


# issues to solve ( --> temp solution):
### need to be specified

## temp: uniformPriorParameter --> uniformPriorParameter <- 0.125 (unform distribution for 8 cues)
uniformPriorParameter <- 1

### temp: myMulti? --> myMulti = 1
myMulti = 1


### temp: first-time cue
firstTimeCue = array(1, c(nCues, n_subj))

### cue index (vector, indexing next trial for spec cue)
cueIndexVec = array(1, c((t_max-nCues), n_subj))

# datalist

data_ra <- list(nCues = nCues, 
                nSubject = n_subj, 
                nTrials = t_max, 
                painRating = painRating,
                expectPain = expectPain, 
                noxInput = noxInput, 
                cueType = cueType, 
                myMulti = myMulti,
                uniformPriorParameter = uniformPriorParameter,
                firstTimeCue = firstTimeCue,
                cueIndexVec = cueIndexVec,
                initExpect = initExpect)

library(rstan)
test_01 = stan(file.choose(), data = data_ra,
              pars = c('alphaCoefC', 'painError', 'expectError', 'alphaC', 'betaC', 'expectScale', 'painScale'), init = 'random',
              chains = 1, cores = 1, iter = 100, warmup = 20)

#Stan data structure
#data
{
#int<lower=1> nCues; // number of cues
#int<lower=1> nSubject; // number of subjects
#int<lower=1> nTrials; // number of subjects
#matrix[nTrials, nSubject] painRating; // observations
#matrix[nTrials, nSubject] expectPain; // observations
#matrix[nTrials, nSubject] noxInput; // observations
#int<lower=1> myMulti;
#int<lower=1> uniformPriorParameter;
#int cueIndexVec[(nTrials-nCues), nSubject]; // indexing next trial for spec cue
#int firstTimeCue[nCues, nSubject];
#matrix[nTrials, nSubject] cueType; // 
#matrix[nCues,  nSubject] initExpect; //  initial expectation mu = expectation rating
}
``

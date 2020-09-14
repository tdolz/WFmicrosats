#### Running STRUCTURE ######

##created Sept 1, 2020
## from script skreportJAN_28_2020MATTITUCK.R

library('plyr')
library("dplyr")
library("poppr")
library("tidyr")
library('purrr')
library("ggplot2")
library("adegenet")
library("pegas")
library("cowplot")

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")

##### Formating the dataset #####
# We are going to use the doubl0 version. 
wfpop <- read.genalex("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_20_sept20204genalex_doubl0.csv")
wfpop4df <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_20_sept2020_doubl0.csv", header = TRUE) #csv version 

splitStrata(wfpop) <-~Ocean/Bay/Con/Year
setPop(wfpop) <-~Bay

#take out the loci and samples. 
wfpopLD <-genclone2genind(wfpop)
all_loci <- locNames(wfpopLD)# create vector of all loci
removeloc <- c("WF27","WF06","WF32")# create vector containing loci to remove
keeploc <- setdiff(all_loci, removeloc)# create list of loci to keep
wfpopLD <- wfpopLD[loc = keeploc]# filter loci in genind object
wfpopLD <-wfpopLD %>% missingno("loci", cutoff=0.30) # removes loci where overall missingness > 30%, so WF01
wfpopLD <- wfpopLD %>% missingno("geno", cutoff=0.25) # remove samples that are > 25% missing
length(locNames(wfpopLD))# check number of loci in genind obj

#convert to df
df <- genind2df(wfpopLD, usepop=TRUE, oneColPerAll = TRUE)
df$Ind <-rownames(df)
df <-df[,c(34,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33)]
df[df=="NA"] <- 0 # missing data must be 0
df[df==0]<- -9  #just kidding, actually make it -9

write.csv(df, "wfpopLD4STRUCTURE.csv") #to save it. 

#make the strata.scheme
wf.strata <-dplyr::select(df, Ind, pop)
strata.schemes <-wf.strata[,c("Ind","pop")] 
row.names(strata.schemes) <-wf.strata$Ind

#convert to gtypes
wf.gtype <-df2gtypes(df,ploidy=2, id.col=1, loc.col =3, schemes=strata.schemes)
wf.gtype

# Use only 2-5 populations and low burn in because otherwise it takes too long. 

#NO ADMIXTURE RUN This model is appropriate for studying fully discrete populations 
#and is often more powerful than the admixture model at detecting subtle structure.
sr <-structureRun(wf.gtype, k.range=1:20,num.k.rep=10, noadmix=TRUE, seed=10, burnin=15000,numreps=350000)

# Calculate Evanno metrics
evno <- evanno(sr)
evno
ggsave('evno_noadmix.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)


#q.mat1 <-clumpp(sr,k=1)
q.mat2 <-clumpp(sr,k=2)
q.mat3 <-clumpp(sr,k=3)
q.mat4 <-clumpp(sr,k=4)
q.mat5 <-clumpp(sr,k=5)
q.mat6 <-clumpp(sr,k=6)
q.mat7 <-clumpp(sr,k=7)
q.mat8 <-clumpp(sr,k=8)
q.mat9 <-clumpp(sr,k=9)
q.mat10 <-clumpp(sr,k=10)
q.mat11 <-clumpp(sr,k=11)
q.mat12 <-clumpp(sr,k=12)
q.mat13 <-clumpp(sr,k=13)
q.mat14 <-clumpp(sr,k=14)
q.mat15 <-clumpp(sr,k=15)
q.mat16 <-clumpp(sr,k=16)
q.mat17 <-clumpp(sr,k=17)
q.mat18 <-clumpp(sr,k=18)
q.mat19 <-clumpp(sr,k=19)
q.mat20 <-clumpp(sr,k=20)

#plot qmat.
#structurePlot(q.mat1,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
structurePlot(q.mat2,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot2_noadmix.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
C:\Users\tara\Desktop
structurePlot(q.mat3,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot3_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat4,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot4_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat5,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot5_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat6,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot6_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat7,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot7_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat8,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot8_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat9,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot9_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat10,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot10_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat11,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot11_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat12,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot12_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat13,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot13_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat14,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot14_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat15,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot15_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat16,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot16_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat17,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot17_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)



####################################################################################
#ADMIXTURE, INDEPENDENT ALLELE FREQUENCIES. 
set.seed(10)
sr <-structureRun(wf.gtype, k.range=1:20,num.k.rep=10, noadmix=FALSE, seed=10, burnin=15000,numreps=350000)

# Calculate Evanno metrics
evno <- evanno(sr)
evno
ggsave('evno_IndAllelesadmix.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

#q.mat1 <-clumpp(sr,k=1)
q.mat2 <-clumpp(sr,k=2)
q.mat3 <-clumpp(sr,k=3)
q.mat4 <-clumpp(sr,k=4)
q.mat5 <-clumpp(sr,k=5)
q.mat6 <-clumpp(sr,k=6)
q.mat7 <-clumpp(sr,k=7)
q.mat8 <-clumpp(sr,k=8)
q.mat9 <-clumpp(sr,k=9)
q.mat10 <-clumpp(sr,k=10)
q.mat11 <-clumpp(sr,k=11)
q.mat12 <-clumpp(sr,k=12)
q.mat13 <-clumpp(sr,k=13)
q.mat14 <-clumpp(sr,k=14)
q.mat15 <-clumpp(sr,k=15)
q.mat16 <-clumpp(sr,k=16)
q.mat17 <-clumpp(sr,k=17)
q.mat18 <-clumpp(sr,k=18)
q.mat19 <-clumpp(sr,k=19)
q.mat20 <-clumpp(sr,k=20)

#plot qmat.
#structurePlot(q.mat1,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
structurePlot(q.mat2,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot2_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat3,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot3_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat4,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot4_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat5,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot5_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat6,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot6_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat7,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot7_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat8,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot8_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat9,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot9_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat10,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot10_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat11,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot11_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat12,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot12_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat13,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot13_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat14,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot14_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat15,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot15_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat16,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot16_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
structurePlot(q.mat17,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
ggsave('splot17_noadmix', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)
















#ADMIXTURE, CORRELATED ALLELE FREQUENCIES.
set.seed(10)
sr <-structureRun(wf.gtype, k.range=1:20,num.k.rep=10, noadmix=FALSE, seed=10,freqscorr=TRUE, burnin=15000,numreps=350000)

#ADMIXTURE, INDEPENDENT ALLELE FREQUENCIES, LOCATION PRIOR. 
#going to increase the reps & burnin, 
set.seed(10)
sr <-structureRun(wf.gtype, k.range=1:20,num.k.rep=10, pop.prior="locprior",noadmix=FALSE, seed=10, burnin=15000,numreps=350000)
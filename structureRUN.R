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
wfpop <- read.genalex("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrected12satsAug1920204genalex_doubl0.csv")
wfpop4df <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrected_12_satsAug2020_doubl0.csv", header = TRUE) #csv version 

splitStrata(wfpop) <-~Ocean/Bay/Con/Year
setPop(wfpop) <-~Bay

#create a second dataset that removes one loci in the WF27/WF33 pair 
#go back to wfpop for this. 
wfpopLD <-genclone2genind(wfpop)
all_loci <- locNames(wfpopLD)# create vector of all loci
removeloc <- c("WF27")# create vector containing loci to remove
keeploc <- setdiff(all_loci, removeloc)# create list of loci to keep
wfpopLD <- wfpopLD[loc = keeploc]# filter loci in genind object
length(locNames(wfpopLD))# check number of loci in genind obj
#now re-remove the missing data. 
info_table(wfpopLD, plot = TRUE, scaled =FALSE)
wfpopLD <- wfpopLD %>% missingno("geno", cutoff=0.20) # remove samples that are > 20% missing



# save everything first. 

#convert to df
df <- genind2df(wfpopLD, usepop=TRUE, oneColPerAll = TRUE)
df$Ind <-rownames(df)
df <-df[,c(24,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,23)]
df[df=="NA"] <- 0
df[df==0]<- -9  #just kidding, actually make it -9

write.csv(df, "wfpopLD4STRUCTURE.csv")

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
sr <-structureRun(wf.gtype, k.range=1:7,num.k.rep=10, noadmix=TRUE, seed=10, burnin=1000,numreps=3500)

# Calculate Evanno metrics
evno <- evanno(sr)
evno

#q.mat1 <-clumpp(sr,k=1)
q.mat2 <-clumpp(sr,k=2)
q.mat3 <-clumpp(sr,k=3)
q.mat4 <-clumpp(sr,k=4)
q.mat5 <-clumpp(sr,k=5)
q.mat6 <-clumpp(sr,k=6)
q.mat7 <-clumpp(sr,k=7)

#plot qmat.
#structurePlot(q.mat1,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
structurePlot(q.mat2,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
structurePlot(q.mat3,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
structurePlot(q.mat4,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
structurePlot(q.mat5,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
structurePlot(q.mat6,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)
structurePlot(q.mat7,pop.col=3,prob.col=4,sort.probs=TRUE,label.pops=TRUE, horiz=FALSE)

#ADMIXTURE, INDEPENDENT ALLELE FREQUENCIES. 
sr <-structureRun(wf.gtype, k.range=1:7,num.k.rep=10, noadmix=FALSE, seed=10, burnin=1000,numreps=3500)

#ADMIXTURE, CORRELATED ALLELE FREQUENCIES.
set.seed(10)
sr <-structureRun(wf.gtype, k.range=1:7,num.k.rep=10, noadmix=FALSE, seed=10,freqscorr=TRUE, burnin=1000,numreps=3500)

#ADMIXTURE, INDEPENDENT ALLELE FREQUENCIES, LOCATION PRIOR. 
#going to increase the reps & burnin, 
set.seed(10)
sr <-structureRun(wf.gtype, k.range=1:7,num.k.rep=10, pop.prior="locprior",noadmix=FALSE, seed=10, burnin=15000,numreps=350000)

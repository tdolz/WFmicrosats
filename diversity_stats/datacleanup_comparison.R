#Data cleanup rarefaction check
# 5/30/23
# here we check to see if it makes a difference when we do data cleanup
# with only 2016 YOY vs. all and with or without rarefaction of the 2016 YOY. 
# I think it will make a difference but it might be sort of too complicated to do them all separately? 
# We will also do the Mt and Shin cleanups here. 
#
##install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)
library("tidyverse")
library('plyr')
library("dplyr")
library("poppr")
library("tidyr")
library('purrr')
library("ggplot2")
library("adegenet")
library("pegas")
library("lattice")
library("related")
library("cowplot")
library('readr')
library("coin")
library("reshape2")
library("strataG")
library("hierfstat")
library("viridis")
library("forcats")
library("data.table")


#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")

############## THE ORIGINAL WAY ###########################

##### Formating the dataset #####
# We are going to use the doubl0 version. 
##### Formating the dataset #####
wfpop <- read.genalex("./data/popcorrect_17_sept20204genalex_doubl0ABC.csv")
wfpop4df <-read.csv("./data/popcorrect_17_sept2020_doubl0ABC.csv", header = TRUE) #csv

splitStrata(wfpop) <-~Ocean/Bay/Con/Year
setPop(wfpop) <-~Bay

#look at missing data
#Where missing data is greater than 10% for that locus...the locus is not informative for that bay.
info_table(wfpop, plot = TRUE, scaled =FALSE)
ggsave("rawinfotable20.png", path="./diversity_stats/diversity_figs/bay")
#dev.off()

#Create the LD dataset to check for Linkage Disequilibreum 
#Linkage disequilibreum  in a NONE MISSING dataset- have to do this. 
wfpop2CLEAN <- wfpop %>% missingno("geno", cutoff=0.0)
setPop(wfpop2CLEAN) <-~Bay
wfia.pair <-wfpop2CLEAN %>% clonecorrect(strata= ~Bay) %>% pair.ia(quiet=FALSE)
wfia.pair <- seppop(wfpop2CLEAN) %>% lapply(pair.ia) #by bay!
ggsave("rawLD20.png", path="./diversity_stats/diversity_figs/bay")
dev.off()

#create a second dataset that removes one loci in the WF27/WF33 pair 
#also remove WF06 & WF32 due to null alllels, see null_checker script and "microsat record Sept2020.ppt" for details. 
#go back to wfpop for this. 
wfpopLD <-genclone2genind(wfpop)
all_loci <- locNames(wfpopLD)# create vector of all loci
removeloc <- c("WF06","WF27")# create vector containing loci to remove
keeploc <- setdiff(all_loci, removeloc)# create list of loci to keep
wfpopLD <- wfpopLD[loc = keeploc]# filter loci in genind object
length(locNames(wfpopLD))# check number of loci in genind obj

#now re-remove the missing data. 
setPop(wfpop2CLEAN) <-~Bay
info_table(wfpopLD, plot = TRUE, scaled =FALSE)
wfpopLD <-wfpopLD %>% missingno("loci", cutoff=0.30) # removes loci where overall missingness > 30%, so WF01
info_table(wfpopLD, plot = TRUE, scaled =FALSE)
wfpopLD <- wfpopLD %>% missingno("geno", cutoff=0.25) # remove samples that are > 25% missing
info_table(wfpopLD, plot = TRUE, scaled =FALSE)
ggsave("25percutinfotable16LD.png", path="./diversity_stats/diversity_figs/bay")
dev.off()

#HWE Heatmap#
setPop(wfpopLD) <-~Ocean
hw.test(wfpopLD, B=1000) #permutation based
hw.test(wfpop2CLEAN, B=1000)
hw.test(wfpopLD, B=0) #analytical p value
wfhwe.pop <- seppop(wfpopLD) %>% lapply(hw.test)
write.csv(wfhwe.pop, file="./diversity_stats/diversity_output_files/bay/wfhwepop16.csv")
(wfhwe.mat <- sapply(wfhwe.pop, "[", i = TRUE, j = 3)) # Take the third column with all rows ---> output this for supplementary tables.
wfhw.mc <-sapply(wfhwe.pop, "[", i = TRUE, j = 4) #the PR exact based on the Monte carlo test! ----> p.values on the hw.test
wfhw.mc  # this is just the p values. 
alpha  <- 0.05
newmat <- wfhwe.mat
newmat[newmat > alpha] <- 1 #where the p value on the chi square is greater than 0.05, give it a 1.
# so pink means zero, which means p < 0.05, which means out of HWE. 
levelplot(t(newmat),scales=list(x=list(rot=90)))
ggsave("HWEtest16.png", path="./diversity_stats/diversity_figs/bay")
dev.off()



############## YOY16 ONLY NO RARIFACTION ###########################

setPop(wfpop) <-~Bay/Year
#remove all individuals that aren't 2016 YOY
yoy16 <-popsub(wfpop, exclude=c("Mt_2015","Mt_adults","Shin_2017"))
setPop(yoy16) <-~Bay

#look at missing data
#Where missing data is greater than 10% for that locus...the locus is not informative for that bay.
info_table(yoy16, plot = TRUE, scaled =FALSE)
ggsave("rawinfotable_yoy16.png", path="./diversity_stats/diversity_figs/YOY16")
#dev.off()

#Create the LD dataset to check for Linkage Disequilibreum 
#Linkage disequilibreum  in a NONE MISSING dataset- have to do this. 
wfpop2CLEAN <- yoy16 %>% missingno("geno", cutoff=0.0)
setPop(wfpop2CLEAN) <-~Bay
wfia.pair <-wfpop2CLEAN %>% clonecorrect(strata= ~Bay) %>% pair.ia(quiet=FALSE)
wfia.pair <- seppop(wfpop2CLEAN) %>% lapply(pair.ia) #by bay!
ggsave("rawLD_yoy16.png", path="./diversity_stats/diversity_figs/YOY16")
dev.off()

#In this case it's still the WF27/WF33 pair!
#create a second dataset that removes one loci in the WF27/WF33 pair 
#also remove WF06 & WF32 due to null alllels, see null_checker script and "microsat record Sept2020.ppt" for details. 
#go back to wfpop for this. 
yoy16<-genclone2genind(yoy16)
all_loci <- locNames(yoy16)# create vector of all loci
removeloc <- c("WF06","WF27")# create vector containing loci to remove
keeploc <- setdiff(all_loci, removeloc)# create list of loci to keep
yoy16 <- yoy16[loc = keeploc]# filter loci in genind object
length(locNames(yoy16))# check number of loci in genind obj

#now re-remove the missing data. A larger percentage will be removed than previous. 
setPop(yoy16) <-~Bay
info_table(yoy16, plot = TRUE, scaled =FALSE)
yoy16 <-yoy16 %>% missingno("loci", cutoff=0.30) # removes loci where overall missingness > 30%
info_table(yoy16, plot = TRUE, scaled =FALSE)
yoy16 <- yoy16 %>% missingno("geno", cutoff=0.25) # remove samples that are > 25% missing
info_table(yoy16, plot = TRUE, scaled =FALSE)
ggsave("25percutinfotable16LD.png", path="./diversity_stats/diversity_figs/YOY16")
dev.off()

#WF32, Pam21, A441 look really bad here. But WF06 is not so bad. 
#So WF06 is getting removed "due to null alleles" unless we re-run nullchecker. 
#WF27 removed due to LD with WF33
#I think we should seriously consider removing WF32 & A441 because they have over 50% missing in some bays. 
#especially due to imbalance in missing data. 
#but is it really the overall NUMBER of missing individuals or is it the percent because NAP and MOR are so low. 

#let's remove WF32 and A441?
all_loci <- locNames(yoy16)# create vector of all loci
removeloc <- c("WF32","A441")# create vector containing loci to remove
keeploc <- setdiff(all_loci, removeloc)# create list of loci to keep
yoy16 <- yoy16[loc = keeploc]# filter loci in genind object
length(locNames(yoy16))# check number of loci in genind obj
setPop(yoy16) <-~Bay
info_table(yoy16, plot = TRUE, scaled =FALSE)
ggsave("cut_WF32_A441.png", path="./diversity_stats/diversity_figs/YOY16")
#now we are down to 15 loci which I am not thrilled about. 

#HWE Heatmap#
setPop(yoy16) <-~Ocean
hw.test(yoy16, B=1000) #permutation based
hw.test(wfpop2CLEAN, B=1000)
hw.test(yoy16, B=0) #analytical p value
wfhwe.pop <- seppop(yoy16) %>% lapply(hw.test)
write.csv(wfhwe.pop, file="./diversity_stats/diversity_output_files/YOY16/wfhwepop16.csv")
(wfhwe.mat <- sapply(wfhwe.pop, "[", i = TRUE, j = 3)) # Take the third column with all rows ---> output this for supplementary tables.
wfhw.mc <-sapply(wfhwe.pop, "[", i = TRUE, j = 4) #the PR exact based on the Monte carlo test! ----> p.values on the hw.test
wfhw.mc  # this is just the p values. 
alpha  <- 0.05
newmat <- wfhwe.mat
newmat[newmat > alpha] <- 1 #where the p value on the chi square is greater than 0.05, give it a 1.
# so pink means zero, which means p < 0.05, which means out of HWE. 
levelplot(t(newmat),scales=list(x=list(rot=90)))
ggsave("HWEtest16.png", path="./diversity_stats/diversity_figs/YOY16")
dev.off()
# we have a few loci which are out of HWE. So really we have 12?


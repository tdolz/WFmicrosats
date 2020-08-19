### Genetic analysis Bays ###
## created 8-19-2020

## NOTES ##
## includes only 12 loci for now.
## no half 0 genotypes.
## population corrected for new cohort assignment
## based on the script "skreportdnaJAN_28_2020MATTITUCK.R"
## Only including analysis that is for the publication. No exploratory analysis. 


#install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)
library("dplyr")
library("poppr")
library("tidyr")
library('purrr')
library("ggplot2")
library("adegenet")
library("pegas")
library("lattice")
library("related")

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")

##### Formating the dataset #####
wfpop <- read.genalex("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrected12satsAug192020.csv")
wfpop4df <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrected_12_satsAug2020.csv", header = TRUE) #csv version 

splitStrata(wfpop) <-~Ocean/Bay/Con/Year
setPop(wfpop) <-~Bay

#look at missing data
#Where missing data is greater than 10% for that locus...the locus is not informative for that bay.
info_table(wfpop, plot = TRUE, scaled =FALSE)
wfpop2 <- wfpop %>% missingno("geno", cutoff=0.20) # remove samples that are > 20% missing
info_table(wfpop2, plot = TRUE, scaled =FALSE)

#Create the LD dataset to check for Linkage Disequilibreum 
#Linkage disequilibreum  in a NONE MISSING dataset- have to do this. 
wfpop2CLEAN <- wfpop2 %>% missingno("geno", cutoff=0.0)
setPop(wfpop2CLEAN) <-~Bay
wfia.pair <-wfpop2CLEAN %>% clonecorrect(strata= ~Bay) %>% pair.ia(quiet=FALSE)
#wfia.pair <- seppop(wfpop2CLEAN) %>% lapply(pair.ia) #by bay!

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

#HWE Heatmap#
setPop(wfpopLD) <-~Bay
wfhwe.pop <- seppop(wfpopLD) %>% lapply(hw.test)
(wfhwe.mat <- sapply(wfhwe.pop, "[", i = TRUE, j = 3)) # Take the third column with all rows ---> output this for supplementary tables.
wfhw.mc <-sapply(wfhwe.pop, "[", i = TRUE, j = 4) #the PR exact based on the Monte carlo test! ----> p.values on the hw.test
wfhw.mc
alpha  <- 0.05
newmat <- wfhwe.mat
newmat[newmat > alpha] <- 1
levelplot(t(newmat),scales=list(x=list(rot=90)))
# a few loci are out of HWE, but not bad. 
# This is a supplementary figure

#Test HWE over Con/Year
#Test HWE over Year population
setPop(wfpopLD) <-~Bay/Con/Year
wfhwe.pop <- seppop(wfpopLD) %>% lapply(hw.test)
(wfhwe.mat <- sapply(wfhwe.pop, "[", i = TRUE, j = 3)) # Take the third column with all rows
wfhw.mc <-sapply(wfhwe.pop, "[", i = TRUE, j = 4) #the PR exact based on the Monte carlo test! 
wfhw.mc
alpha  <- 0.05
newmat <- wfhwe.mat
newmat[newmat > alpha] <- 1
library("lattice")
levelplot(t(newmat),scales=list(x=list(rot=90)))#slightly different if you do it with genind object (which it calls for) instead of genclone.

######## Pairwise Relatedness ######

#There are many different estimators for relatedness. In previous scripts we have tried more robust ways. 
# for now, Shannon's way. https://gist.github.com/sjoleary/3efd4a7d56b115fad319781298765a31

#first convert to the right data format
setPop(wfpopLD) <-~Bay
df <-genind2df(wfpopLD, usepop = FALSE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(23,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)]
df[df=="NA"] <- 0 # missing data must be 0
write.table(df, "scratch",row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) #write it.
genotypedata <- readgenotypedata("scratch")# import input file as list (gdata, nloci, nalleles, ninds, freqs)

# pairwise relatedness (Lynch & Ritland 1999)
relatedness_lynchrd <- coancestry(genotypedata$gdata,
                                  lynchrd = 1)

# write relatedness to file
relatedn <- relatedness_lynchrd$relatedness %>%
  dplyr::select(pair.no, ind1.id, ind2.id, lynchrd)

#not sure what's going on here. 
library("readr")
write_delim(relatedn, "pairwise_relatedness")
# write inbreeding to file
inbreed <- relatedness_lynchrd$inbreeding %>%
  dplyr::select(ind.id, L3, LH) %>%
  dplyr::rename(INDV = ind.id)
write_delim(inbreed, "inbreeding")

#the simulation
#change bootstrap back to 100 because 1000 just takes too long. 
cosim <-compareestimators(relatedness_lynchrd, ninds=100)
cosim

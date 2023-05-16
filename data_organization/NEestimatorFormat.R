#Format data for NEestimator. 

## September 1, 2020
## originally from the script "skreportJAN_28_MATTITUCK

#install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)
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

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")

##### Formating the dataset ##### use the doubl0 version.
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
#ggsave("20percutinfotable11LD.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()


#Try another way to format the input file for NEestimator
#could do genind to genpop? 
setPop(wfpopLD) <- ~Bay/Con/Year
yoy2016 <- popsub(wfpopLD, blacklist = c("Shin_1_2017", "Shin_2_2017","Mt_1_2015", "Mt_2_2015", "Mt_3_2015", "Mt_5_2015", "Mt_3_2016","Mt_4_2016", "Mt_5_2016")) #remove Mattituck adults and remove non- 2016 fish. 
setPop(yoy2016) <-~Bay

yoy2016gen <- genind2genpop(yoy2016, pop=yoy2016@pop)

#library(devtools)
#install_github("romunov/zvau")
#install.packages("remotes")
#remotes::install_github("romunov/zvau")

library("zvau")
#writeGenPop(yoy2016, "yoy2016MT.gen", "NEEstimatorMT")

#this works, but you have to go to the file and manually fix the number of zeroes. 
# remember to structure it like yoy2016.txt
#NAP MOR JAM SHIN

#use the file yoy2016.gen

#make a no population file

setPop(wfpopLD) <- ~Bay/Con/Year
yoy2016 <- popsub(wfpopLD, blacklist = c("Shin_1_2017", "Shin_2_2017","Mt_1_2015", "Mt_2_2015", "Mt_3_2015", "Mt_5_2015", "Mt_3_2016","Mt_4_2016", "Mt_5_2016")) #remove Mattituck adults and remove non- 2016 fish. 
setPop(yoy2016) <-~Ocean
#writeGenPop(yoy2016, "yoy2016OceanMT.gen", "NEEstimator")

#lets compare Mattituck Ne from 2015 data vs. from 2016 data to see if they agree. 
setPop(wfpopLD) <- ~Bay/Con/Year
yoyMT2015 <- popsub(wfpopLD, sublist = c("Mt_1_2015", "Mt_2_2015")) #only include 2015 mattituck YOY
setPop(yoyMT2015) <-~Bay
#writeGenPop(yoyMT2015, "yoy2015MT.gen", "NEEstimator")

#compare Shinnecock Ne from 2016 vs 2017
setPop(wfpopLD) <- ~Bay/Con/Year
yoyShi2017 <- popsub(wfpopLD, sublist = c("Shin_1_2017", "Shin_2_2017")) #only include 2015 mattituck YOY
setPop(yoyShi2017) <-~Bay
#writeGenPop(yoyShi2017, "yoy2017Shi.gen", "NEEstimator")

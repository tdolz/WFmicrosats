####FST comparisons ####### 
#Nov 1, 2020, updated 5/29/23

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
library('readr')
library("reshape2")
library("strataG")

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")

##### Formating the dataset #####
wfpop <- read.genalex("./data/popcorrect_17_sept20204genalex_doubl0ABC.csv")
wfpop4df <-read.csv("./data/popcorrect_17_sept2020_doubl0ABC.csv", header = TRUE) #csv

splitStrata(wfpop) <-~Ocean/Bay/Con/Year
setPop(wfpop) <-~Bay

#clean dataset. 
wfpopLD <-genclone2genind(wfpop)
all_loci <- locNames(wfpopLD)# create vector of all loci
removeloc <- c("WF06","WF27")# create vector containing loci to remove
keeploc <- setdiff(all_loci, removeloc)# create list of loci to keep
wfpopLD <- wfpopLD[loc = keeploc]
wfpopLD <-wfpopLD %>% missingno("loci", cutoff=0.30) # removes loci where overall missingness > 30%, so WF01
wfpopLD <- wfpopLD %>% missingno("geno", cutoff=0.25) # remove samples that are > 25% missing

## Subset databases
#shinnecock only database
setPop(wfpopLD) <-~Bay
shinco <-popsub(wfpopLD, sublist=c("Shin"))
setPop(shinco) <-~Bay/Con/Year

#mattituck only database, exclude MT_5
setPop(wfpopLD) <-~Bay/Con
mtco <-popsub(wfpopLD, sublist=c("Mt_1","Mt_2","Mt_3","Mt_4"))
setPop(mtco) <-~Bay/Con/Year

####HWE - New tests with Benjamini-Hochberg.####
#bh adjust function
pradj <-function(df){
  df <-as.data.frame(df) %>%tibble::rownames_to_column() %>%
    mutate(Pr.adj=p.adjust(Pr.exact, method="BH"))%>%mutate(sig=ifelse(Pr.adj < 0.05, "sig","not sig"))}

###### FST ########
#Weir and Cockheram FST global & pairwise. - probably your best bet. 

#bays.
setPop(wfpopLD) <-~Bay
wf.bay <- genind2gtypes(wfpopLD)
popStruct.b <- popStructTest(wf.bay, nrep = 1000, quietly = TRUE)
popStruct.b
#BH correct the Fst pvals
mtres <- as.data.frame(popStruct.b$pairwise$result)
mtres <-dplyr::select(mtres, strata.1,strata.2, Fst, Fst.p.val) %>% dplyr::rename(Pr.exact=Fst.p.val)
mtres <-pradj(mtres)
write.csv(mtres, file="./diversity_stats/diversity_output_files/bay/bay_fst_notrarified.csv")

#BAY - CON Mattituck
#explore structure within Mattituck
setPop(mtco) <-~Bay/Con
wf.M2 <- genind2gtypes(mtco)
popStruct.M <- popStructTest(wf.M2, nrep = 1000, quietly = TRUE)
popStruct.M
#BH correct the Fst pvals
mtres <- as.data.frame(popStruct.M$pairwise$result)
mtres <-dplyr::select(mtres, strata.1,strata.2, Fst, Fst.p.val) %>% dplyr::rename(Pr.exact=Fst.p.val)
mtres <-pradj(mtres)
write.csv(mtres, file="./diversity_stats/diversity_output_files/Mattituck/MT_fst_BayCon_notrarified.csv")

#BAY/YEAR/CON Mattituck. 
setPop(mtco) <-~Bay/Con/Year
wf.M2 <- genind2gtypes(mtco)
popStruct.M <- popStructTest(wf.M2, nrep = 1000, quietly = TRUE)
popStruct.M
#BH correct the Fst pvals
mtres <- as.data.frame(popStruct.M$pairwise$result)
mtres <-dplyr::select(mtres, strata.1,strata.2, Fst, Fst.p.val) %>% dplyr::rename(Pr.exact=Fst.p.val)
mtres <-pradj(mtres)
write.csv(mtres, file="./diversity_stats/diversity_output_files/Mattituck/MT_fst_BayConYEAR_notrarified.csv")

#Shinnecock. 
setPop(shinco) <-~Bay/Con/Year
wf.s2 <- genind2gtypes(shinco)
#if for some annoying reason, that doesn't work, try converting from df
df <-genind2df(shinco, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(36,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)]
df[df=="NA"] <- 0 # missing data must be 0
wf.g <-df2gtypes(df,ploidy=2)
#back to business 
popStruct.S <- popStructTest(wf.g, nrep = 1000, quietly = TRUE)
popStruct.S
mtres <- as.data.frame(popStruct.S$pairwise$result)
mtres <-dplyr::select(mtres, strata.1,strata.2, Fst, Fst.p.val) %>% dplyr::rename(Pr.exact=Fst.p.val)
mtres <-pradj(mtres)
write.csv(mtres, file="./diversity_stats/diversity_output_files/Mattituck/SHIN_fst_BayConYEAR_notrarified.csv")

#### The bays comparison where we do a rareifaction########## 
setPop(wfpopLD) <-~Bay/Year
#remove all individuals that aren't 2016 YOY
yoy16 <-popsub(wfpopLD, exclude=c("Mt_2015","Mt_adults","Shin_2017"))
setPop(yoy16) <-~Bay
#convert to a df. 
df <-genind2df(yoy16, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(36,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)]
df[df=="NA"] <- 0 # missing data must be 0

#randomly sample from within shinnecock and mattituck 100x
rare.out <-list()
for(i in 1:2){
  df.split <-split(df, df$pop)
  new.shin <-sample_n(df.split$Shin,30)
  new.mt <-sample_n(df.split$Mt,30)
  rare.bays <-bind_rows(new.shin,new.mt,df.split$Jam,df.split$Mor,df.split$Nap)
  wf.rare <-df2gtypes(rare.bays,ploidy=2)
  popStruct.b <- popStructTest(wf.rare, nrep = 1000, quietly = TRUE) 
  mtres <- as.data.frame(popStruct.b$pairwise$result)
  mtres <-dplyr::select(mtres, strata.1,strata.2, Fst, Fst.p.val) %>% dplyr::rename(Pr.exact=Fst.p.val)%>%
    mutate(count=i)
  #mtres <-pradj(mtres) #don't need to do this yet.
  rare.out[[i]]<-mtres
}

mean.fst <-bind_rows(rare.out)%>%group_by(strata.1,strata.2)%>%
  summarize(mean.fst=mean(Fst), sd.fst=sd(Fst), mean.pval=mean(Pr.exact), sd.pval=sd(Pr.exact),.groups="keep")
  # i don't actually think you can combine the pval this way, but let's calculate them anyway. 
write.csv(mean.fst, file="./diversity_stats/diversity_output_files/YOY16_bay/bay16_fst_rarified.csv")

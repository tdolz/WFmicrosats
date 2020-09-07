### checking for null alleles #####

#install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)
library('plyr')
library("dplyr")
library("poppr")
library("tidyr")
library("ggplot2")
library("adegenet")
library("PopGenReport")
library('readr')

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")

##### Formating the dataset #####
# We are going to use the doubl0 version. 
wfpop <- read.genalex("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_20_sept20204genalex_doubl0.csv")
wfpop4df <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_20_sept2020_doubl0.csv", header = TRUE) #csv version 

splitStrata(wfpop) <-~Ocean/Bay/Con/Year
setPop(wfpop) <-~Bay
###### The null alleles analysis ######
# this was mainly done offline. The full explanation can be found in the powerpoint "microsat record Sept2020.ppt"
# eventually we decided to exclude WF06 and WF32 due to null alleles. a further two candidates would include WF517 and WF3.
# We will keep a separate script for checking for nulls

setPop(wfpop) <-~Bay
wf.gen <-genclone2genind(wfpop) 
popgenreport(wf.gen,mk.counts=TRUE,mk.locihz = TRUE, mk.fst=TRUE, mk.allele.dist=TRUE, mk.null.all=TRUE,mk.allel.rich = TRUE,mk.differ.stats = TRUE,path.pgr=getwd(),mk.Rcode=TRUE,mk.pdf=TRUE )

#can also just check for nulls
wf.nulls <- null.all(wf.gen)
#If the 95% confidence interval includes zero, it indicates that the frequency of null alleles at a locus does not significantly differ from zero.

#Check each population separately for null alleles, as suggested. 
popNames(wfpop)
wf.shin <-popsub(wf.gen, sublist=c("Shin"))
#wf.shin.nulls <-null.all(wf.shin)
popgenreport(wf.shin,mk.counts=TRUE,mk.locihz = TRUE, mk.fst=TRUE, mk.allele.dist=TRUE, mk.null.all=TRUE,mk.allel.rich = TRUE,mk.differ.stats = TRUE,path.pgr="/Users/tdolan/Documents/WIP research/microsats",mk.Rcode=TRUE,mk.pdf=TRUE )

#Check each population separately for null alleles, as suggested. 
popNames(wfpop)
wf.shin <-popsub(wf.gen, sublist=c("Shin"))
#wf.shin.nulls <-null.all(wf.shin)
popgenreport(wf.shin,mk.counts=TRUE,mk.locihz = TRUE, mk.fst=TRUE, mk.allele.dist=TRUE, mk.null.all=TRUE,mk.allel.rich = TRUE,mk.differ.stats = TRUE,path.pgr="/Users/tdolan/Documents/WIP research/microsats",mk.Rcode=TRUE,mk.pdf=TRUE )

#hacky way to convert to genepop file for offline checking. 
wf.gp <-wfpop4df %>% unite(J42, J42.1,J42.2,sep='',remove=TRUE) %>%unite(WF22,WF22.1,WF22.2,sep='',remove=TRUE)%>%
  unite(PAM27, PAM27.1,PAM27.2,sep='',remove=TRUE) %>%unite(PAM79,PAM79.1,PAM79.2,sep='',remove=TRUE)%>%
  unite(WF27, WF27.1,WF27.2,sep='',remove=TRUE) %>%unite(WF16,WF16.1,WF16.2,sep='',remove=TRUE)%>%
  unite(WF33, WF33.1,WF33.2,sep='',remove=TRUE) %>%unite(WF3,WF3.1,WF3.2,sep='',remove=TRUE)%>%
  unite(WF517, WF517.1,WF517.2,sep='',remove=TRUE) %>%unite(WF196,WF196.1,WF196.2,sep='',remove=TRUE)%>%
  unite(WF223, WF223.1,WF223.2,sep='',remove=TRUE) %>%unite(WF421,WF421.1,WF421.2,sep='',remove=TRUE)%>%
  unite(A441, A441.1,A441.2,sep='',remove=TRUE) %>%unite(PAM21,PAM21.1,PAM21.2,sep='',remove=TRUE)%>%
  unite(Psy087, Psy087.1,Psy087.2,sep='',remove=TRUE) %>%unite(Psy022,Psy022.1,Psy022.2,sep='',remove=TRUE)%>%
  unite(WF06, WF06.1,WF06.2,sep='',remove=TRUE) %>%unite(WF12,WF12.1,WF12.2,sep='',remove=TRUE)%>%
  unite(WF01, WF01.1,WF01.2,sep='',remove=TRUE) %>%unite(WF32,WF32.1,WF32.2,sep='',remove=TRUE)
wf.gp[wf.gp == "00"] <- "000000"
wf.gp[is.na(wf.gp)]<-"000000"
wf.gp <-dplyr::select(wf.gp, -Pop)
#write_delim(wf.gp, "wfgp")

#try checking for nulls with genepop
library("genepop")
#wf.genpop <-genind2genpop(wf.gen)
#genepop file had to be manually reformatted to check for 00 instead of 000 and also the pop line. 
nulls(inputFile="/Users//tdolan/Documents//R-Github//WFmicrosats/wfgp", outputFile= "/Users//tdolan/Documents//R-Github//WFmicrosats/wfnulls", verbose=TRUE)


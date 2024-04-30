#loading plots 
#3/19/23

library("poppr")
library("tidyr")
library('purrr')
library("ggplot2")
library("adegenet")
library("pegas")
library("cowplot")
library("vcfR")
library("adegenet")
library("forcats")
detach(package:plyr)    
library(dplyr)


###############functions##################
#which quantile function. 
which.quantile <- function (x, probs, na.rm = FALSE){
  if (! na.rm & any (is.na (x)))
    return (rep (NA_integer_, length (probs)))
  o <- order (x)
  n <- sum (! is.na (x))
  o <- o [seq_len (n)]
  nppm <- n * probs - 0.5
  j <- floor(nppm)
  h <- ifelse((nppm == j) & ((j%%2L) == 0L), 0, 1)
  j <- j + h
  j [j == 0] <- 1
  o[j]
}

###################
#scaffold for data frame of summary results. 
r_names <- c("LD1","locus max sum loading","loading","mean loading","std dev loading","locus min sum loading","loading","allele max loading","loading",
             "LD2","locus max sum loading","loading","mean sum loading","std dev sum loading","locus min sum loading","loading","allele max loading","loading")


#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")


##### Formating the dataset #######################################
# We are going to use the doubl0 version. 
wfpop <- read.genalex("/Users//tdolan/Documents//R-Github//WFmicrosats/data/popcorrect_17_sept20204genalex_doubl0ABC.csv")
wfpop4df <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/data/popcorrect_17_sept2020_doubl0ABC.csv", header = TRUE) #csv version 

splitStrata(wfpop) <-~Ocean/Bay/Con/Year
setPop(wfpop) <-~Bay

#clean the dataset
wfpopLD <-genclone2genind(wfpop)
all_loci <- locNames(wfpopLD)# create vector of all loci
removeloc <- c("WF06")# create vector containing loci to remove
keeploc <- setdiff(all_loci, removeloc)# create list of loci to keep
wfpopLD <- wfpopLD[loc = keeploc]# filter loci in genind object
wfpopLD <-wfpopLD %>% missingno("loci", cutoff=0.30) # removes loci where overall missingness > 30%, so WF01
wfpopLD <- wfpopLD %>% missingno("geno", cutoff=0.25) # remove samples that are > 25% missing
length(locNames(wfpopLD))# check number of loci in genind obj

################################# 2016 YOY ONLY ########################################
wfpopLD3 <-wfpopLD
setPop(wfpopLD3) <-~Bay/Con/Year
#exclude everything but yoy 2016
wfbays16 <-popsub(wfpopLD3, exclude =c("Mt_1_2015", "Mt_2_2015", "Mt_3_adults", "Mt_4_adults",
                                       "Mt_5_2015","Mt_5_2016","Shin_1_2017","Shin_2_2017"))
#set back to bay
setPop(wfbays16) <-~Bay

#################################################################################################################
################## LI BAYS 2016 A PRIORI ###################################################################################################
ldvar <-read.csv(file="./DAPC/DAPC_figs/bay16_loadings_apriori.csv",header=T)
rownames(ldvar)<-ldvar[,1]
ldvar <-ldvar[,-1]

#summarize loadings
lds <- ldvar %>% as.data.frame()%>%rownames_to_column(var="Allele")%>%
  tidyr::separate(col=Allele, into=c("locus","allele"), remove=FALSE)%>%mutate(locus=as.factor(locus))

lds.summary <-lds %>%dplyr::group_by(locus)%>%summarize(sum.ld1 = sum(LD1), sum.ld2=sum(LD2))

#save these
lds_baysAP <-lds
ldsum_baysAP <-lds.summary

#ummary results
btwn_bay_apriori <-c(NA,paste(lds.summary[which.max(lds.summary$sum.ld1),][[1]]),paste(lds.summary[which.max(lds.summary$sum.ld1),][[2]]),paste(mean(lds$LD1)),paste(sd(lds$LD1)),
             paste(lds.summary[which.min(lds.summary$sum.ld1),][[1]]),paste(lds.summary[which.min(lds.summary$sum.ld1),][[2]]),paste(lds[which.max(lds$LD1),][[1]]),paste(lds[which.max(lds$LD1),][[4]]),
             NA, paste(lds.summary[which.max(lds.summary$sum.ld2),][[1]]),paste(lds.summary[which.max(lds.summary$sum.ld2),][[3]]),paste(mean(lds$LD2)),paste(sd(lds$LD2)),
             paste(lds.summary[which.min(lds.summary$sum.ld2),][[1]]),paste(lds.summary[which.min(lds.summary$sum.ld2),][[3]]),paste(lds[which.max(lds$LD2),][[1]]),paste(lds[which.max(lds$LD2),][[5]]))


#LOADING PLOTS - diverging alleles between a priori groupings. 
#we want a loading plot where the top 1% of loadings are shown.

#Axis1
thresh1 = lds[which.quantile(lds$LD1,0.99),][[4]]
#what alleles are above the thresh?
b <-filter(lds, LD1 >= thresh1);b
althresh1 <-unique(b$Allele)

#Axis2 
thresh2 = lds[which.quantile(lds$LD2,0.99),][[5]]
#what alleles are above the thresh?
b <-filter(lds, LD2 >= thresh2);b
althresh2 <-unique(b$Allele)

#LOADING PLOT - diverging alleles between a priori groupings. 
set.seed(4)
par(mfrow=c(1,2))
#AXIS 1
contrib1 <- loadingplot(ldvar, axis=1,
                        thres=thresh1, lab.jitter=0, main="Axis 1", xlab="Alleles", ylab="DAPC Loadings")
#AXIS 2
contrib2 <- loadingplot(ldvar, axis=2,
                        thres=thresh2, lab.jitter=4, main="Axis 2", xlab="Alleles", ylab="DAPC Loadings")
#save 900x500

althresh1
althresh2

#Axis 1
#co<-round(length(althresh1)/2)
#par(mfrow=c(2,co))
par(mfrow=c(4,3))
for (i in 1:length(althresh1)){
  p <- tab(genind2genpop(wfbays16[loc=i]),freq=TRUE)
  p.2 <-p[,i]
  matplot(p.2,pch=19, type="b",
          xlab="pop",ylab="allele frequency", xaxt="n",
          cex=1, main=althresh1[i])
  axis(side=1, at=1:5, lab=c("Jam","Mor","Mt","Nap","Shin"))
}
#dev.off()

#Axis 2
#co<-round(length(althresh1)/2)
#par(mfrow=c(2,co))
#par(mfrow=c(4,3))
for (i in 1:length(althresh2)){
  p <- tab(genind2genpop(wfbays16[loc=i]),freq=TRUE)
  p.2 <-p[,i]
  matplot(p.2,pch=19, type="b",
          xlab="pop",ylab="allele frequency", xaxt="n",
          cex=1, main=althresh2[i])
  axis(side=1, at=1:5, lab=c("Jam","Mor","Mt","Nap","Shin"))
}
#900 by 900
#save as 
dev.off()

############################################################################################################################################
################## LI BAYS k=2 ###################################################################################################
ldvar <-read.csv(file="./DAPC/DAPC_figs/bay16_loadings_k2.csv",header=T)
rownames(ldvar)<-ldvar[,1]

#summarize loadings
lds <- ldvar %>% as.data.frame()%>%rownames_to_column(var="Allele")%>%
  tidyr::separate(col=Allele, into=c("locus","allele"), remove=FALSE)%>%mutate(locus=as.factor(locus))%>%select(-X)

lds.summary <-lds %>%dplyr::group_by(locus)%>%summarize(sum.ld1 = sum(LD1))

lds_baysk2 <-lds
ldsum_baysk2 <-lds.summary

#ummary results
btwn_bay_k2 <-c(NA,paste(lds.summary[which.max(lds.summary$sum.ld1),][[1]]),paste(lds.summary[which.max(lds.summary$sum.ld1),][[2]]),paste(mean(lds$LD1)),paste(sd(lds$LD1)),
                     paste(lds.summary[which.min(lds.summary$sum.ld1),][[1]]),paste(lds.summary[which.min(lds.summary$sum.ld1),][[2]]),paste(lds[which.max(lds$LD1),][[1]]),paste(lds[which.max(lds$LD1),][[4]]),
                     NA, NA,NA,NA,NA,NA,NA,NA,NA)

#LOADING PLOTS - diverging alleles between a priori groupings. 
#we want a loading plot where the top 1% of loadings are shown.

#Axis1
thresh1 = lds[which.quantile(lds$LD1,0.99),][[4]]
#what alleles are above the thresh?
b <-filter(lds, LD1 >= thresh1);b
althresh1 <-unique(b$Allele)


#LOADING PLOT - diverging alleles between a priori groupings. 
ldvar <-as.data.frame(ldvar)%>%select(-X)
set.seed(4)
par(mfrow=c(1,2))
#AXIS 1
contrib1 <- loadingplot(ldvar, axis=1,
                        thres=thresh1, lab.jitter=4, main="Axis 1", xlab="Alleles", ylab="DAPC Loadings")

#save 900x500

dev.off()

#Axis 1
#co<-round(length(althresh1)/2)
#par(mfrow=c(2,co))
par(mfrow=c(4,3))
for (i in 1:length(althresh1)){
  p <- tab(genind2genpop(wfbays16[loc=i]),freq=TRUE)
  p.2 <-p[,i]
  matplot(p.2,pch=19, type="b",
          xlab="pop",ylab="allele frequency", xaxt="n",
          cex=1, main=althresh1[i])
  axis(side=1, at=1:5, lab=c("Jam","Mor","Mt","Nap","Shin"))
}


dev.off()

############################################################################################################################################
################## MATTITUCK WITH 2015 APRIORI ###################################################################################################
wfpopLD4 <-wfpopLD
setPop(wfpopLD4) <-~Bay
mt <-popsub(wfpopLD4, sublist =c("Mt"))
setPop(mt) <-~Bay/Con/Year
#exclude unidentified adults but we will include the 2015 YOY for now. 
mt <-popsub(mt,exclude=c("Mt_5_2015", "Mt_5_2016"))

### read in data
ldvar <-read.csv(file="./DAPC/DAPC_figs/MT15_loadings_apriori.csv",header=T)
rownames(ldvar)<-ldvar[,1]
ldvar <-ldvar[,-1]

#summarize loadings
lds <- ldvar %>% as.data.frame()%>%rownames_to_column(var="Allele")%>%
  tidyr::separate(col=Allele, into=c("locus","allele"), remove=FALSE)%>%mutate(locus=as.factor(locus))

lds.summary <-lds %>%dplyr::group_by(locus)%>%summarize(sum.ld1 = sum(LD1), sum.ld2=sum(LD2))

#save these
lds_mtAP <-lds
ldsum_mtAP <-lds.summary


#ummary results
Mt15_apriori <-c(NA,paste(lds.summary[which.max(lds.summary$sum.ld1),][[1]]),paste(lds.summary[which.max(lds.summary$sum.ld1),][[2]]),paste(mean(lds$LD1)),paste(sd(lds$LD1)),
                     paste(lds.summary[which.min(lds.summary$sum.ld1),][[1]]),paste(lds.summary[which.min(lds.summary$sum.ld1),][[2]]),paste(lds[which.max(lds$LD1),][[1]]),paste(lds[which.max(lds$LD1),][[4]]),
                     NA, paste(lds.summary[which.max(lds.summary$sum.ld2),][[1]]),paste(lds.summary[which.max(lds.summary$sum.ld2),][[3]]),paste(mean(lds$LD2)),paste(sd(lds$LD2)),
                     paste(lds.summary[which.min(lds.summary$sum.ld2),][[1]]),paste(lds.summary[which.min(lds.summary$sum.ld2),][[3]]),paste(lds[which.max(lds$LD2),][[1]]),paste(lds[which.max(lds$LD2),][[5]]))


#LOADING PLOTS - diverging alleles between a priori groupings. 
#we want a loading plot where the top 1% of loadings are shown.

#Axis1
thresh1 = lds[which.quantile(lds$LD1,0.99),][[4]]
#what alleles are above the thresh?
b <-filter(lds, LD1 >= thresh1);b
althresh1 <-unique(b$Allele)

#Axis2 
thresh2 = lds[which.quantile(lds$LD2,0.99),][[5]]
#what alleles are above the thresh?
b <-filter(lds, LD2 >= thresh2);b
althresh2 <-unique(b$Allele)

par(mfrow=c(1,2))
#LOADING PLOT - diverging alleles between a priori groupings. 
set.seed(4)
#AXIS 1
contrib1 <- loadingplot(ldvar, axis=1,
                        thres=thresh1, lab.jitter=4, main="Axis 1", xlab="Alleles", ylab="DAPC Loadings")
#AXIS 2
contrib2 <- loadingplot(ldvar, axis=2,
                        thres=thresh2, lab.jitter=4, main="Axis 2", xlab="Alleles", ylab="DAPC Loadings")
#save 900x500

althresh1
althresh2

#Axis 1
#co<-round(length(althresh1)/2)
#par(mfrow=c(2,co))
par(mfrow=c(4,3))
for (i in 1:length(althresh1)){
  p <- tab(genind2genpop(mt[loc=i]),freq=TRUE)
  p.2 <-p[,i]
  matplot(p.2,pch=19, type="b",
          xlab="pop",ylab="allele frequency", xaxt="n",
          cex=1, main=althresh1[i])
  axis(side=1, at=1:6, lab=c("Early 2015","Early 2016","Late 2015","Late 2016","Migrants","Residents"))}
#dev.off()

#Axis 2
#co<-round(length(althresh1)/2)
#par(mfrow=c(2,co))
#par(mfrow=c(4,3))
for (i in 1:length(althresh2)){
  p <- tab(genind2genpop(mt[loc=i]),freq=TRUE)
  p.2 <-p[,i]
  matplot(p.2,pch=19, type="b",
          xlab="pop",ylab="allele frequency", xaxt="n",
          cex=1, main=althresh2[i])
  axis(side=1, at=1:6, lab=c("Early 2015","Early 2016","Late 2015","Late 2016","Migrants","Residents"))
}

dev.off()

############################################################################################################################################
################## MATTITUCK WITH 2015 k2 ###################################################################################################
ldvar <-read.csv(file="./DAPC/DAPC_figs/MT15_loadings_k2.csv",header=T)
rownames(ldvar)<-ldvar[,1]


#summarize loadings
lds <- ldvar %>% as.data.frame()%>%rownames_to_column(var="Allele")%>%
  tidyr::separate(col=Allele, into=c("locus","allele"), remove=FALSE)%>%mutate(locus=as.factor(locus))%>%select(-X)

lds.summary <-lds %>%dplyr::group_by(locus)%>%summarize(sum.ld1 = sum(LD1))

#save these
lds_mtk2 <-lds
ldsum_mtk2 <-lds.summary


#ummary results
Mt15_k2 <-c(NA,paste(lds.summary[which.max(lds.summary$sum.ld1),][[1]]),paste(lds.summary[which.max(lds.summary$sum.ld1),][[2]]),paste(mean(lds$LD1)),paste(sd(lds$LD1)),
            paste(lds.summary[which.min(lds.summary$sum.ld1),][[1]]),paste(lds.summary[which.min(lds.summary$sum.ld1),][[2]]),paste(lds[which.max(lds$LD1),][[1]]),paste(lds[which.max(lds$LD1),][[4]]),
            NA, NA, NA, NA,NA,NA,NA,NA,NA)


#LOADING PLOTS - diverging alleles between a priori groupings. 
#we want a loading plot where the top 1% of loadings are shown.

#Axis1
thresh = lds[which.quantile(lds$LD1,0.99),][[4]]
#what alleles are above the thresh?
b <-filter(lds, LD1 >= thresh);b
althresh1 <-unique(b$Allele)

set.seed(4)
par(mfrow=c(1,2))
ldvar<-select(ldvar,-X)
#AXIS 1
contrib1 <- loadingplot(ldvar, axis=1,
                        thres=thresh, lab.jitter=4, main="Axis 1", xlab="Alleles", ylab="DAPC Loadings")

#save 900x500

althresh1


#Axis 1
#co<-round(length(althresh1)/2)
#par(mfrow=c(2,co))
par(mfrow=c(4,3))
for (i in 1:length(althresh1)){
  p <- tab(genind2genpop(mt[loc=i]),freq=TRUE)
  p.2 <-p[,i]
  matplot(p.2,pch=19, type="b",
          xlab="pop",ylab="allele frequency", xaxt="n",
          cex=1, main=althresh1[i])
  axis(side=1, at=1:6, lab=c("Early 2015","Early 2016","Late 2015","Late 2016","Migrants","Residents"))
}


  dev.off()
  

############################################################################################################################################
################## SHINNECOCK APRIORI ###################################################################################################
wfpopLD5 <-wfpopLD
setPop(wfpopLD5) <-~Bay
sh <-popsub(wfpopLD5, sublist =c("Shin"))
setPop(sh) <-~Bay/Con/Year
unique(sh@pop)


ldvar <-read.csv(file="./DAPC/DAPC_figs/shin_loadings_apriori.csv",header=T)
rownames(ldvar)<-ldvar[,1]
ldvar <-ldvar[,-1]

#summarize loadings
lds <- ldvar %>% as.data.frame()%>%rownames_to_column(var="Allele")%>%
  tidyr::separate(col=Allele, into=c("locus","allele"), remove=FALSE)%>%mutate(locus=as.factor(locus))

lds.summary <-lds %>%dplyr::group_by(locus)%>%summarize(sum.ld1 = sum(LD1), sum.ld2=sum(LD2))

lds_shinAP <-lds
ldsum_ShinAP <-lds.summary


#ummary results
shin_apriori  <-c(NA,paste(lds.summary[which.max(lds.summary$sum.ld1),][[1]]),paste(lds.summary[which.max(lds.summary$sum.ld1),][[2]]),paste(mean(lds$LD1)),paste(sd(lds$LD1)),
                                  paste(lds.summary[which.min(lds.summary$sum.ld1),][[1]]),paste(lds.summary[which.min(lds.summary$sum.ld1),][[2]]),paste(lds[which.max(lds$LD1),][[1]]),paste(lds[which.max(lds$LD1),][[4]]),
                                  NA, paste(lds.summary[which.max(lds.summary$sum.ld2),][[1]]),paste(lds.summary[which.max(lds.summary$sum.ld2),][[3]]),paste(mean(lds$LD2)),paste(sd(lds$LD2)),
                                  paste(lds.summary[which.min(lds.summary$sum.ld2),][[1]]),paste(lds.summary[which.min(lds.summary$sum.ld2),][[3]]),paste(lds[which.max(lds$LD2),][[1]]),paste(lds[which.max(lds$LD2),][[5]]))



#LOADING PLOTS - diverging alleles between a priori groupings. 
#we want a loading plot where the top 1% of loadings are shown.

#Axis1
thresh1 = lds[which.quantile(lds$LD1,0.99),][[4]]
#what alleles are above the thresh?
b <-filter(lds, LD1 >= thresh1);b
althresh1 <-unique(b$Allele)
#Axis2 
thresh2 = lds[which.quantile(lds$LD2,0.99),][[5]]
#what alleles are above the thresh?
b <-filter(lds, LD2 >= thresh2);b
althresh2 <-unique(b$Allele)

set.seed(4)
par(mfrow=c(1,2))
#AXIS 1
contrib1 <- loadingplot(ldvar, axis=1,
                        thres=thresh1, lab.jitter=4, main="Axis 1", xlab="Alleles", ylab="DAPC Loadings")

#AXIS 2
contrib2 <- loadingplot(ldvar, axis=2,
                        thres=thresh2, lab.jitter=4, main="Axis 2", xlab="Alleles", ylab="DAPC Loadings")
#save 900x500

althresh1
althresh2

#Axis 1
#co<-round(length(althresh1)/2)
#par(mfrow=c(2,co))
par(mfrow=c(4,3))
for (i in 1:length(althresh1)){
  p <- tab(genind2genpop(sh[loc=i]),freq=TRUE)
  p.2 <-p[,i]
  matplot(p.2,pch=19, type="b",
          xlab="pop",ylab="allele frequency", xaxt="n",
          cex=1, main=althresh1[i])
  axis(side=1, at=1:4, lab=c("Early 2016","Early 2017","Late 2016","Late 2017"))}
  #dev.off()
  
  #Axis 2
  #co<-round(length(althresh1)/2)
  #par(mfrow=c(2,co))
  #par(mfrow=c(4,3))
  for (i in 1:length(althresh2)){
    p <- tab(genind2genpop(sh[loc=i]),freq=TRUE)
    p.2 <-p[,i]
    matplot(p.2,pch=19, type="b",
            xlab="pop",ylab="allele frequency", xaxt="n",
            cex=1, main=althresh2[i])
    axis(side=1, at=1:4, lab=c("Early 2016","Early 2017","Late 2016","Late 2017"))
  }
  
  
  dev.off()
  
  
  
############################################################################################################################################
################## SHINNECOCK k2 ###################################################################################################
ldvar <-read.csv(file="./DAPC/DAPC_figs/shin_loadings_k2.csv",header=T)
  rownames(ldvar)<-ldvar[,1]
  
  #summarize loadings
  lds <- ldvar %>% as.data.frame()%>%rownames_to_column(var="Allele")%>%
    tidyr::separate(col=Allele, into=c("locus","allele"), remove=FALSE)%>%mutate(locus=as.factor(locus))%>%select(-X)
  
  lds.summary <-lds %>%dplyr::group_by(locus)%>%summarize(sum.ld1 = sum(LD1))

  lds_shink2 <-lds
  ldsum_shink2 <-lds.summary
  

  #ummary results
  shin_k2 <-c(NA,paste(lds.summary[which.max(lds.summary$sum.ld1),][[1]]),paste(lds.summary[which.max(lds.summary$sum.ld1),][[2]]),paste(mean(lds$LD1)),paste(sd(lds$LD1)),
              paste(lds.summary[which.min(lds.summary$sum.ld1),][[1]]),paste(lds.summary[which.min(lds.summary$sum.ld1),][[2]]),paste(lds[which.max(lds$LD1),][[1]]),paste(lds[which.max(lds$LD1),][[4]]),
              NA, NA,NA,NA,NA,NA,NA,NA,NA)
  
  
  #LOADING PLOTS - diverging alleles between a priori groupings. 
  #we want a loading plot where the top 1% of loadings are shown.
  
  #Axis1
  thresh1 = lds[which.quantile(lds$LD1,0.99),][[4]]
  #what alleles are above the thresh?
  b <-filter(lds, LD1 >= thresh1);b
  althresh1 <-unique(b$Allele)
  
  set.seed(4)
  par(mfrow=c(1,2))
  ldvar<-as.data.frame(ldvar)%>%select(-X)
  #AXIS 1
  contrib1 <- loadingplot(ldvar, axis=1,
                          thres=thresh1, lab.jitter=4, main="Axis 1", xlab="Alleles", ylab="DAPC Loadings")

  
  #Axis 1
  #co<-round(length(althresh1)/2)
  #par(mfrow=c(2,co))
  par(mfrow=c(4,3))
  for (i in 1:length(althresh1)){
    p <- tab(genind2genpop(sh[loc=i]),freq=TRUE)
    p.2 <-p[,i]
    matplot(p.2,pch=19, type="b",
            xlab="pop",ylab="allele frequency", xaxt="n",
            cex=1, main=althresh1[i])
    axis(side=1, at=1:4, lab=c("Early 2016","Early 2017","Late 2016","Late 2017"))}
  
  
  dev.off()

######################################### combine all the summary stats. ######################################################################
 
lstats <-bind_cols(r_names,btwn_bay_apriori,btwn_bay_k2,Mt15_apriori,Mt15_k2,shin_apriori,shin_k2)  
names(lstats) <-c("stat","btwn_bay_apriori","btwn_bay_k2","Mt15_apriori","Mt15_k2","shin_apriori","shin_k2")  

write.csv(lstats,file="./DAPC/DAPC_figs/loading_stats.csv")


################ LDS VS ALLELIC DIVERSITY ANALYSIS ##################################################

#"Summing the variance contributions from all alleles for each locus revealed a positive and roughly linear 
#relationship between loading values and the total number of alleles for all DAPC analyses"

lds_baysAP <- lds_baysAP %>%mutate(analysis="LI Bays",clust="a priori",LDtotal=LD1+LD2)
lds_mtAP <- lds_mtAP %>%mutate(analysis="Mattituck",clust="a priori",LDtotal=LD1+LD2)
lds_shinAP <- lds_shinAP %>%mutate(analysis="Shinnecock",clust="a priori",LDtotal=LD1+LD2)

lds_baysk2<-lds_baysk2 %>%mutate(LD2=NA,analysis="LI Bays",clust="k-means",LDtotal=LD1)
lds_mtk2<-lds_mtk2 %>%mutate(LD2=NA,analysis="Mattituck",clust="k-means",LDtotal=LD1)
lds_shink2<-lds_shink2 %>%mutate(LD2=NA,analysis="Shinnecock",clust="k-means",LDtotal=LD1)

lds_all <-bind_rows(lds_baysAP,lds_mtAP,lds_shinAP,lds_baysk2,lds_mtk2,lds_shink2)

lds_all.sum <-lds_all %>%group_by(analysis,clust,locus)%>%summarize(sum.ld1 = sum(LD1), 
              n.alleles=n_distinct(Allele), sum.ld2=sum(LD2,na.rm=T), sum.bothLD=sum(LDtotal), .groups="keep")

ap.plot <-lds_all.sum%>%
  ggplot(aes(n.alleles,sum.bothLD, label=locus))+
  geom_point()+
  geom_label_repel(aes(label = locus),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  geom_smooth(method = "lm", color="black", size=0.5,)+
  #stat_cor(label.y = 10)+ #this means at 35th unit in the y axis, the r squared and p value will be shown
  #stat_regline_equation(label.y = 10)+ #this means at 30th unit regresion line equation will be shown
  #stat_poly_line() +  #this would be great but its from ggpmisc so i cant use it because that package is blocked. 
  #stat_poly_eq(use_label(c("eq", "R2"))) +#this would be great but its from ggpmisc so i cant use it because that package is blocked. 
  xlab("Allelic richness")+ylab("DAPC loadings")+
  facet_grid(analysis~clust)+
  theme_classic()+
  theme(
    #strip.background =element_rect(fill="transparent"),
    panel.spacing = unit(1.5, "mm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size=10), axis.text.x = element_text(size=10),
    strip.text.x=element_text(size=10),strip.text.y=element_text(size=10),
    #panel.background = element_rect(fill = "white",colour = "dark grey"),
    plot.background = element_rect(fill = NA,colour = "black", size=0.5),
    panel.border = element_rect(colour="black",fill=NA, size=0.5)
  );ap.plot
ggsave("loadingsVrichness.png", width=12, height=10, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

summary(lm(sum.bothLD~n.alleles,data=filter(lds_all.sum,clust=="a priori" & analysis=="LI Bays")))
summary(lm(sum.bothLD~n.alleles,data=filter(lds_all.sum,clust=="a priori" & analysis=="Mattituck")))
summary(lm(sum.bothLD~n.alleles,data=filter(lds_all.sum,clust=="a priori" & analysis=="Shinnecock")))

summary(lm(sum.bothLD~n.alleles,data=filter(lds_all.sum,clust=="k-means" & analysis=="LI Bays")))
summary(lm(sum.bothLD~n.alleles,data=filter(lds_all.sum,clust=="k-means" & analysis=="Mattituck")))
summary(lm(sum.bothLD~n.alleles,data=filter(lds_all.sum,clust=="k-means" & analysis=="Shinnecock")))





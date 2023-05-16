#Created March 17, 2023 from script ClusterDAPCnew.R for Shinnecock and Mattituck only. 

library('plyr')
library("dplyr")
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


#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")


##### Formating the dataset #####
# We are going to use the doubl0 version. 
wfpop <- read.genalex("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_17_sept20204genalex_doubl0ABC.csv")
wfpop4df <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_17_sept2020_doubl0ABC.csv", header = TRUE) #csv version 

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


####################################### ALL GROUPS SHINNECOCK ONLY ####################################
setPop(wfpopLD) <-~Bay/Con/Year
#exclude everything but yoy 2016
shinpop <-popsub(wfpopLD, sublist=c("Shin_1_2016", "Shin_1_2017", "Shin_2_2016", "Shin_2_2017"))

#first set to all groups. 
setPop(shinpop) <-~Bay/Con/Year

shincolors <-c("#29bf12", "#abff4f", "#3f8efc", "#3b28cc")
tolcolors <-c("#332288","#117733","#88CCEE","#DDCC77","#CC6677","#AA4499","#882255","#44AA99")

##################################### A priori assignment ###########################################
setPop(shinpop) <-~Bay/Con/Year
shin_dapc <- dapc(shinpop, n.pca = 150, n.da = 2) 
as <- optim.a.score(shin_dapc)#how many pc? 42

### Apriori with cross validation, by bay ###
df <- genind2df(shinpop, usepop=TRUE, oneColPerAll = TRUE)
df$Ind <-rownames(df)
df <-df[,c(38,1:37)]
df[df=="NA"] <- 0
#make the strata.scheme
wf.strata <-dplyr::select(df, Ind, pop)
strata.schemes <-wf.strata[,c("Ind","pop")] 
row.names(strata.schemes) <-wf.strata$Ind
GRP <- strata.schemes %>% dplyr::select(Ind, pop)

set.seed(1213)
#xval the assignment. 
X <- scaleGen(shinpop, NA.method = "mean")
xval <- xvalDapc(X, GRP$pop, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 50, xval.plot = TRUE)
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`); retain 
#10pc

#re-do plots with 10 pc based on xval
shin_dapc <- dapc(shinpop, n.pca = 10, n.da = 2) 
scatter(shin_dapc, posi.da="topleft",scree.pca=TRUE,posi.pca="topright",col=shincolors,legend=FALSE,solid=1.0) #adjust the solidity value for better plots. 

#re-do plots with 40 pc based on 2nd best xval
shin_dapc <- dapc(shinpop, n.pca = 40, n.da = 2) 
scatter(shin_dapc, posi.da="topleft",scree.pca=TRUE,posi.pca="topright",col=shincolors,legend=FALSE,solid=1.0) #adjust the solidity value for better plots. 
#let's go with 40. 

## Compare geographic regions of origin and assigned group memberships
grp_membership <- as.data.frame(shin_dapc$posterior) %>% tibble::rownames_to_column("Ind") %>%
  gather(key = GRP, value = MEMBSHIP, 2:5) %>% arrange(MEMBSHIP)%>%
  left_join(GRP)
#NEW compopolot - fewer PC - It looks worse.... 
memprob <-ggplot(grp_membership, aes(x = Ind, y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = shincolors) +theme_cowplot() ;memprob
ggsave("shinBYC_apriori_membship_40pc.png", width=10, height=5, path="/Users/tdolan/documents/WIP research/microsats/microsat_figs/figs2023")

#to see how many were assigned to their population of origin, compare how many individuals have a GRP that matches POP
boo <-filter(grp_membership, GRP==pop) %>%n_distinct()
toto <-n_distinct(grp_membership$Ind)
boo/toto
boo<-filter(grp_membership, GRP==pop)
mean(boo$MEMBSHIP) #the percent out of total that are assigned to the correct group
sem(boo$MEMBSHIP) #standard error of the mean percentage


###########################################naive k means clustering #####################################
############ let algo determine number of clusters ###########
grp_BIC <- find.clusters.genind(shinpop, stat = "BIC",choose.n.clust = TRUE,max.n.clust = 40)
bay_kstat <-grp_BIC$Kstat
bay_kstat

grp_BIC <- find.clusters.genind(shinpop, stat = "BIC",n.clust = 2, n.pca=40)
bay_kstat <-grp_BIC$Kstat

# # perform stratified cross validation
X <- scaleGen(shinpop, NA.method = "mean")
xval <- xvalDapc(X, grp_BIC$grp, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$`Mean Successful Assignment by Number of PCs of PCA`
# identify optimum number of PCs to retain  --> 30
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)
retain #50 

#make scatterplot
# perform DAPC using k-mean clusters as groups
dapck <-dapc(shinpop, grp_BIC$grp, n.pca = 50, n.da=2) #not going to specifiy the number of da, retain 4.
scatter(dapck,grp=grp_BIC$grp, posi.da="topright", bg="white", pch=17:22, scree.pca=TRUE, posi.pca="topleft", col=tolcolors)

table(pop(wfpopLD), grp_BIC$grp)

#cluster assignments relative to original group assignments
kgrp <- strata.schemes %>% dplyr::select(Ind, pop)
kgrp$grp <-grp_BIC$grp
# show cluster distribution across populations

## Compare geographic regions of origin and assigned group memberships
kgrps <-kgrp %>% group_by(pop,grp)%>%summarize(group=n_distinct(Ind),.groups="drop")
kgrps %>%
  ggplot(aes(x=pop, y=group, fill=grp)) +
  scale_fill_manual(values=tolcolors)+ylab("number of individuals")+
  geom_bar(stat="identity",alpha=0.7, position=position_dodge())+
  theme(axis.text.x = element_text(angle = 90),axis.text = element_text(size = 14, face="bold"),axis.title = element_text(size = 14,face="bold"),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),panel.spacing = unit(1, "lines"))

#now create the compoplot
kgrp <- as.data.frame(dapck$grp) %>% tibble::rownames_to_column(var="Ind") #we will use this later. 
GRP <- strata.schemes %>% dplyr::select(Ind, pop)

kgrp <-left_join(kgrp,GRP, by="Ind")
kgrps <-kgrp %>% group_by(pop,`dapck$grp`)%>%summarize(group=n_distinct(Ind),.groups="drop")

grp_membership <- as.data.frame(dapck$posterior) %>% tibble::rownames_to_column("Ind") 
grp_membership <- grp_membership %>%
  gather(key = GRP, value = MEMBSHIP, 2:ncol(grp_membership)) %>% arrange(MEMBSHIP)%>%
  left_join(GRP)%>%mutate(groupnum=as.integer(GRP)) %>%arrange(Ind,MEMBSHIP)

#order the bars by their predominant association
mjgrp <-grp_membership%>%group_by(Ind)%>%summarize(mxmem =max(MEMBSHIP))
mjgrp<-left_join(grp_membership,mjgrp)%>%
  mutate(major.group=ifelse(mxmem==MEMBSHIP,groupnum,NA))%>%
  select(Ind,major.group)
mjgrp <-mjgrp[!is.na(mjgrp$major.group), ] 
grp_membership <-left_join(grp_membership, mjgrp)

#structure plot
memprob <-ggplot(grp_membership, aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = tolcolors)+
  theme_cowplot();memprob
ggsave("shinBYC_kmeans_membship_50pc_2clust.png", width=10, height=5, path="/Users/tdolan/documents/WIP research/microsats/microsat_figs/figs2023")
#dev.off()
#dim 510x345 gets rid of the white lines.

########### specify number of clusters ##########################
############ let algo determine number of clusters ###########
grp_BIC <- find.clusters.genind(shinpop, stat = "BIC",n.clust = 4, n.pca=50)
bay_kstat <-grp_BIC$Kstat

# # perform stratified cross validation
X <- scaleGen(shinpop, NA.method = "mean")
xval <- xvalDapc(X, grp_BIC$grp, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$`Mean Successful Assignment by Number of PCs of PCA`
# identify optimum number of PCs to retain  --> 30
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)
retain #50 

#make scatterplot
# perform DAPC using k-mean clusters as groups
dapck <-dapc(shinpop, grp_BIC$grp, n.pca = 50, n.da=2) #not going to specifiy the number of da, retain 4.
scatter(dapck,grp=grp_BIC$grp, posi.da="topright", bg="white", pch=17:22, scree.pca=TRUE, posi.pca="topleft", col=tolcolors)

#cluster assignments relative to original group assignments
kgrp <- strata.schemes %>% dplyr::select(Ind, pop)
kgrp$grp <-grp_BIC$grp
# show cluster distribution across populations

## Compare geographic regions of origin and assigned group memberships
kgrps <-kgrp %>% group_by(pop,grp)%>%summarize(group=n_distinct(Ind),.groups="drop")
kgrps %>%
  ggplot(aes(x=pop, y=group, fill=grp)) +
  scale_fill_manual(values=tolcolors)+ylab("number of individuals")+
  geom_bar(stat="identity",alpha=0.7, position=position_dodge())+
  theme(axis.text.x = element_text(angle = 90),axis.text = element_text(size = 14, face="bold"),axis.title = element_text(size = 14,face="bold"),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),panel.spacing = unit(1, "lines"))

#now create the compoplot
kgrp <- as.data.frame(dapck$grp) %>% tibble::rownames_to_column(var="Ind") #we will use this later. 
GRP <- strata.schemes %>% dplyr::select(Ind, pop)

kgrp <-left_join(kgrp,GRP, by="Ind")
kgrps <-kgrp %>% group_by(pop,`dapck$grp`)%>%summarize(group=n_distinct(Ind),.groups="drop")

grp_membership <- as.data.frame(dapck$posterior) %>% tibble::rownames_to_column("Ind") 
grp_membership <- grp_membership %>%
  gather(key = GRP, value = MEMBSHIP, 2:ncol(grp_membership)) %>% arrange(MEMBSHIP)%>%
  left_join(GRP)%>%mutate(groupnum=as.integer(GRP)) %>%arrange(Ind,MEMBSHIP)

#order the bars by their predominant association
mjgrp <-grp_membership%>%group_by(Ind)%>%summarize(mxmem =max(MEMBSHIP))
mjgrp<-left_join(grp_membership,mjgrp)%>%
  mutate(major.group=ifelse(mxmem==MEMBSHIP,groupnum,NA))%>%
  select(Ind,major.group)
mjgrp <-mjgrp[!is.na(mjgrp$major.group), ] 
grp_membership <-left_join(grp_membership, mjgrp)

#structure plot
memprob <-ggplot(grp_membership, aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = tolcolors)+
  theme_cowplot();memprob
ggsave("shinBYC_kmeans_membship_50pc_4clust.png", width=10, height=5, path="/Users/tdolan/documents/WIP research/microsats/microsat_figs/figs2023")
#dev.off()
#dim 510x345 gets rid of the white lines.



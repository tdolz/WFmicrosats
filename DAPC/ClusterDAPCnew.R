#Created February 9, 2023 from script clusteringDAPC.R in order to figure out what we did last time. 
#updated May 15, 2023


### Here, in this sheet we are trying to do the kmeans tuning. 
### in a second sheet we willcompare bays based on 2016 YOY alone
### in a third sheet, we will compare within Mattituck and Within Shinnecock. 

library('plyr')
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
library("dplyr")


#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")


##### Formating the dataset #####
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


########################################## BAYS DAPC ALL #################################################################

tolcolors <-c("#332288","#117733","#88CCEE","#DDCC77","#CC6677","#AA4499","#882255","#44AA99")
drabcolors <-c("#d0d1e6","#a6bddb", "#67a9cf", "#1c9099", "#016450")

### A priori assignment ###
setPop(wfpopLD) <-~Bay
bay_dapc <- dapc(wfpopLD, n.pca = 150, n.da = 2) 
as <- optim.a.score(bay_dapc)#how many pc? 150

#PC plots 150pc
scatter(bay_dapc, posi.da="topleft",scree.pca=TRUE,posi.pca="topright",col=drabcolors,legend=FALSE,solid=1.0) #adjust the solidity value for better plots. 
scatter(bay_dapc,1,1, bg="white", scree.da=FALSE, legend=TRUE, solid=.4,col=drabcolors,)

### Apriori with cross validation, by bay ###
setPop(wfpopLD)<-~Bay
df <- genind2df(wfpopLD, usepop=TRUE, oneColPerAll = TRUE)
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
X <- scaleGen(wfpopLD, NA.method = "mean")
xval <- xvalDapc(X, GRP$pop, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 50, xval.plot = TRUE)
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`); retain 
#100pc

#re-do plots with 100 pc based on xval
bay_dapc <- dapc(wfpopLD, n.pca = 100, n.da = 2) 
scatter(bay_dapc, posi.da="topleft",scree.pca=TRUE,posi.pca="topright",col=drabcolors,legend=FALSE,solid=1.0) #adjust the solidity value for better plots. 

## Compare geographic regions of origin and assigned group memberships
grp_membership <- as.data.frame(bay_dapc$posterior) %>% tibble::rownames_to_column("Ind") %>%
  gather(key = GRP, value = MEMBSHIP, 2:6) %>% arrange(MEMBSHIP)%>%
  left_join(GRP)
#NEW compopolot - fewer PC - It looks worse.... 
memprob <-ggplot(grp_membership, aes(x = Ind, y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = drabcolors) +theme_cowplot() ;memprob
#to see how many were assigned to their population of origin, compare how many individuals have a GRP that matches POP
boo <-filter(grp_membership, GRP==pop) %>%n_distinct()
toto <-n_distinct(grp_membership$Ind)
boo/toto
boo<-filter(grp_membership, GRP==pop)
mean(boo$MEMBSHIP) #the percent out of total that are assigned to the correct group
#44.5%
sem(boo$MEMBSHIP) #standard error of the mean percentage
 #0.011 

################################ NAIVE K MEANS CLUSTERING BAY ######################################################
setPop(wfpopLD)<-~Bay
set.seed(1713)

############ let algo determine number of clusters ###########
#first let the algo choose
grp_BIC <- find.clusters.genind(wfpopLD, stat = "BIC",choose.n.clust = FALSE,max.n.clust = 40)
plot(grp_BIC$Kstat)


grp_BIC <- find.clusters.genind(wfpopLD, stat = "BIC",choose.n.clust = TRUE,max.n.clust = 40)
bay_kstat <-grp_BIC$Kstat
bay_kstat
# naive number of groups is 3 BIC=586.1856, with 100PC
# naive number of groups is 2 BIC=625.8268, with 150PC
# overall best number of groups is 3, I am going to say because BIC is lower
grp_BIC <- find.clusters.genind(wfpopLD, stat = "BIC",
                                n.clust = 4, n.pca=150)
bay_kstat <-grp_BIC$Kstat

# # perform stratified cross validation
X <- scaleGen(wfpopLD, NA.method = "mean")
xval <- xvalDapc(X, grp_BIC$grp, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$`Mean Successful Assignment by Number of PCs of PCA`
# identify optimum number of PCs to retain  --> 100, mean successful assignment 74.6%
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)
retain #100 or 150PC is best. 

#make scatterplot
# perform DAPC using k-mean clusters as groups
dapck <-dapc(wfpopLD, grp_BIC$grp, n.pca = retain, n.da=2) #not going to specifiy the number of da, retain 4.
scatter(dapck,grp=grp_BIC$grp, posi.da="topright", bg="white", pch=17:22, scree.pca=TRUE, posi.pca="topleft", col=tolcolors)

table(pop(wfpopLD), grp_BIC$grp)

#cluster assignments relative to original group assignments
kgrp <- strata.schemes %>% dplyr::select(Ind, pop)
kgrp$grp <-grp_BIC$grp
# show cluster distribution across populations

## Compare geographic regions of origin and assigned group memberships
kgrps <-kgrp %>% group_by(pop,grp)%>%dplyr::summarize(group=n_distinct(Ind),.groups="drop")
kgrps %>%
  ggplot(aes(x=pop, y=group, fill=grp)) +
  #scale_fill_manual(values=tolcolors)+ylab("number of individuals")+
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
mjgrp <-grp_membership%>%group_by(Ind)%>%dplyr::summarize(mxmem =max(MEMBSHIP))
mjgrp<-left_join(grp_membership,mjgrp)%>%
mutate(major.group=ifelse(mxmem==MEMBSHIP,groupnum,NA))%>%
  dplyr::select(Ind,major.group)
mjgrp <-mjgrp[!is.na(mjgrp$major.group), ] 
grp_membership <-left_join(grp_membership, mjgrp)

#structure plot
memprob <-ggplot(grp_membership, aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  #scale_fill_manual(values = tolcolors)+
  theme_cowplot();memprob
#ggsave("bays_kmeans_membship_150pc_3clust.png", width=10, height=5, path="/Users/tdolan/documents/WIP research/microsats/microsat_figs/figs2023")
#dev.off()
#dim 510x345 gets rid of the white lines.

#### Now, structure plot within bay strata ################
## Create new strata schemes 
wfpopLD2 <-wfpopLD
setPop(wfpopLD2)<-~Bay/Con/Year
df <- genind2df(wfpopLD2, usepop=TRUE, oneColPerAll = TRUE)
df$Ind <-rownames(df)
df <-df[,c(38,1:37)]
df[df=="NA"] <- 0
#make the strata.scheme
wf.strata <-dplyr::select(df, Ind, pop)
strata.schemes <-wf.strata[,c("Ind","pop")] 
row.names(strata.schemes) <-wf.strata$Ind
GRP <- strata.schemes %>% dplyr::select(Ind, pop)

#repeat the process for joining DAPC results to strata.schemes
kgrp <- as.data.frame(dapck$grp) %>% tibble::rownames_to_column(var="Ind") #we will use this later. 
kgrp <-left_join(kgrp,GRP, by="Ind")
kgrps <-kgrp %>% group_by(pop,`dapck$grp`)%>%summarize(group=n_distinct(Ind),.groups="drop")
grp_membership <- as.data.frame(dapck$posterior) %>% tibble::rownames_to_column("Ind") 
grp_membership <- grp_membership %>%
  gather(key = GRP, value = MEMBSHIP, 2:ncol(grp_membership)) %>% arrange(MEMBSHIP)%>%
  left_join(GRP)%>%mutate(groupnum=as.integer(GRP)) %>%arrange(Ind,MEMBSHIP)

#order the bars by their predominant association
mjgrp <-grp_membership%>%group_by(Ind)%>%dplyr::summarize(mxmem =max(MEMBSHIP))
mjgrp<-left_join(grp_membership,mjgrp)%>%
  mutate(major.group=ifelse(mxmem==MEMBSHIP,groupnum,NA))%>%
  dplyr::select(Ind,major.group)
mjgrp <-mjgrp[!is.na(mjgrp$major.group), ] 
grp_membership <-left_join(grp_membership, mjgrp)

#structure plot
memprob <-ggplot(grp_membership, aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  #scale_fill_manual(values = tolcolors)+
  theme_cowplot();memprob


#############################################################################################
################################# TRY THIS AGAIN BUT 2016 YOY ONLY ########################################
wfpopLD3 <-wfpopLD
setPop(wfpopLD3) <-~Bay/Con/Year
#exclude everything but yoy 2016
wfbays16 <-popsub(wfpopLD3, exclude =c("Mt_1_2015", "Mt_2_2015", "Mt_3_adults", "Mt_4_adults",
                                            "Mt_5_2015","Mt_5_2016","Shin_1_2017","Shin_2_2017"))
#set back to bay
setPop(wfbays16) <-~Bay

## A priori assignment ###
drabcolors <-c("#d0d1e6","#a6bddb", "#67a9cf", "#1c9099", "#016450")

### A priori assignment ###
bay_dapc <- dapc(wfbays16, n.pca = 150, n.da = 2) 
as <- optim.a.score(bay_dapc)#optimum number of pc is 61
bay_dapc <- dapc(wfbays16, n.pca = 61, n.da = 2) 

#PC plots 61pc
scatter(bay_dapc, posi.da="bottomleft",scree.pca=TRUE,posi.pca="bottomright",col=drabcolors,legend=FALSE,solid=1.0) #adjust the solidity value for better plots. 
scatter(bay_dapc,1,1, bg="white", scree.da=FALSE, legend=TRUE, solid=.4,col=drabcolors,)

### Apriori with cross validation, by bay ###
setPop(wfbays16)<-~Bay
df <- genind2df(wfbays16, usepop=TRUE, oneColPerAll = TRUE)
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
X <- scaleGen(wfbays16, NA.method = "mean")
xval <- xvalDapc(X, GRP$pop, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 50, xval.plot = TRUE)
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`); retain 
#100pc

#re-do plots with 100 pc based on xval
bay_dapc <- dapc(wfbays16, n.pca = 100, n.da = 2) 
scatter(bay_dapc, posi.da="topleft",scree.pca=TRUE,posi.pca="topright",col=drabcolors,legend=FALSE,solid=1.0) #adjust the solidity value for better plots. 

## Compare geographic regions of origin and assigned group memberships
grp_membership <- as.data.frame(bay_dapc$posterior) %>% tibble::rownames_to_column("Ind") %>%
  gather(key = GRP, value = MEMBSHIP, 2:6) %>% arrange(MEMBSHIP)%>%
  left_join(GRP)
#NEW compopolot - fewer PC - It looks worse.... 
memprob <-ggplot(grp_membership, aes(x = Ind, y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = drabcolors) +theme_cowplot() ;memprob
ggsave("bays16_apriori_membship_100pc.png", width=10, height=5, path="/Users/tdolan/documents/WIP research/microsats/microsat_figs/figs2023")

############ kmeans ###########
as <- optim.a.score(bay_dapc)#optimum number of pc is 61
## Let the algo choose k
grp_BIC <- find.clusters.genind(wfbays16, stat = "BIC",choose.n.clust = FALSE,max.n.clust = 40)
plot(grp_BIC$Kstat)
#now look at optimal pc
grp_BIC <- find.clusters.genind(wfbays16, stat = "BIC",choose.n.clust = TRUE,max.n.clust = 40)
#now you choose k
grp_BIC <- find.clusters.genind(wfbays16, stat = "BIC",
                                n.clust = 4, n.pca=150)
bay_kstat <-grp_BIC$Kstat

# # perform stratified cross validation
X <- scaleGen(wfbays16, NA.method = "mean")
xval <- xvalDapc(X, grp_BIC$grp, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$`Mean Successful Assignment by Number of PCs of PCA`
# identify optimum number of PCs to retain  --> 100, mean successful assignment 74.6%
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)
retain #100 or 150PC is best. 

#make scatterplot
# perform DAPC using k-mean clusters as groups
dapck <-dapc(wfbays16, grp_BIC$grp, n.pca = retain, n.da=2) #not going to specifiy the number of da, retain 4.
scatter(dapck,grp=grp_BIC$grp, posi.da="topright", bg="white", pch=17:22, scree.pca=TRUE, posi.pca="topleft", col=tolcolors)

table(pop(wfbays16), grp_BIC$grp)

#cluster assignments relative to original group assignments
kgrp <- strata.schemes %>% dplyr::select(Ind, pop)
kgrp$grp <-grp_BIC$grp
# show cluster distribution across populations

## Compare geographic regions of origin and assigned group memberships
kgrps <-kgrp %>% group_by(pop,grp)%>%dplyr::summarize(group=n_distinct(Ind),.groups="drop")
kgrps %>%
  ggplot(aes(x=pop, y=group, fill=grp)) +
  #scale_fill_manual(values=tolcolors)+ylab("number of individuals")+
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
mjgrp <-grp_membership%>%group_by(Ind)%>%dplyr::summarize(mxmem =max(MEMBSHIP))
mjgrp<-left_join(grp_membership,mjgrp)%>%
  mutate(major.group=ifelse(mxmem==MEMBSHIP,groupnum,NA))%>%
  dplyr::select(Ind,major.group)
mjgrp <-mjgrp[!is.na(mjgrp$major.group), ] 
grp_membership <-left_join(grp_membership, mjgrp)

#structure plot
memprob <-ggplot(grp_membership, aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  #scale_fill_manual(values = tolcolors)+
  theme_cowplot();memprob
#ggsave("bays16_kmeans_membship_150pc_3clust.png", width=10, height=5, path="/Users/tdolan/documents/WIP research/microsats/microsat_figs/figs2023")
#dev.off()
#dim 510x345 gets rid of the white lines.


####################### WITHIN MATTITUCK ##########################################
##################################################################################

wfpopLD4 <-wfpopLD
setPop(wfpopLD4) <-~Bay
mt <-popsub(wfpopLD4, sublist =c("Mt"))
setPop(mt) <-~Bay/Con/Year

### A priori assignment ###
bay_dapc <- dapc(mt, n.pca = 150, n.da = 2) 
as <- optim.a.score(bay_dapc)#optimum number of pc is 53
bay_dapc <- dapc(mt, n.pca = 53, n.da = 2) 

#PC plots 61pc
scatter(bay_dapc, posi.da="bottomleft",scree.pca=TRUE,posi.pca="bottomright",col=drabcolors,legend=FALSE,solid=1.0) #adjust the solidity value for better plots. 
scatter(bay_dapc,1,1, bg="white", scree.da=FALSE, legend=TRUE, solid=.4,col=drabcolors,)

### Apriori with cross validation, by bay ###
setPop(mt)<-~Bay/Con/Year
df <- genind2df(mt, usepop=TRUE, oneColPerAll = TRUE)
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
X <- scaleGen(mt, NA.method = "mean")
xval <- xvalDapc(X, GRP$pop, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 50, xval.plot = TRUE)
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`); retain 
#redo scatter
bay_dapc <- dapc(mt, n.pca = retain, n.da = 2) 
scatter(bay_dapc, posi.da="topleft",scree.pca=TRUE,posi.pca="topright",legend=FALSE,solid=1.0) #adjust the solidity value for better plots. 

## Compare geographic regions of origin and assigned group memberships
grp_membership <- as.data.frame(bay_dapc$posterior) %>% tibble::rownames_to_column("Ind") %>%
  gather(key = GRP, value = MEMBSHIP, 2:9) %>% arrange(MEMBSHIP)%>%
  left_join(GRP)
#compopolot 
memprob <-ggplot(grp_membership, aes(x = Ind, y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  theme_cowplot() ;memprob
ggsave("mt_apriori_membship_80pc.png", width=10, height=5, path="/Users/tdolan/documents/WIP research/microsats/microsat_figs/figs2023")

###################### kmeans ###################################
## Let the algo choose k
grp_BIC <- find.clusters.genind(mt, stat = "BIC",choose.n.clust = FALSE,max.n.clust = 40)

grp_BIC <- find.clusters.genind(mt, stat = "BIC",choose.n.clust = TRUE,max.n.clust = 40)
# the above was really just to look at optimal pc/
# now we specify clusters
grp_BIC <- find.clusters.genind(mt, stat = "BIC",
                                n.clust = 9, n.pca=100)
bay_kstat <-grp_BIC$Kstat

# perform stratified cross validation
X <- scaleGen(mt, NA.method = "mean")
xval <- xvalDapc(X, grp_BIC$grp, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$`Mean Successful Assignment by Number of PCs of PCA`
# identify optimum number of PCs to retain
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)
retain  

#make scatterplot
# perform DAPC using k-mean clusters as groups
dapck <-dapc(mt, grp_BIC$grp, n.pca = retain, n.da=2) #not going to specifiy the number of da, retain 4.
scatter(dapck,grp=grp_BIC$grp, posi.da="topright", bg="white", pch=17:22, scree.pca=TRUE, posi.pca="topleft", col=tolcolors)

table(pop(mt), grp_BIC$grp)

#cluster assignments relative to original group assignments
kgrp <- strata.schemes %>% dplyr::select(Ind, pop)
kgrp$grp <-grp_BIC$grp
# show cluster distribution across populations
kgrps <-kgrp %>% group_by(pop,grp)%>%dplyr::summarize(group=n_distinct(Ind),.groups="drop")

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
mjgrp <-grp_membership%>%group_by(Ind)%>%dplyr::summarize(mxmem =max(MEMBSHIP))
mjgrp<-left_join(grp_membership,mjgrp)%>%
  mutate(major.group=ifelse(mxmem==MEMBSHIP,groupnum,NA))%>%
  dplyr::select(Ind,major.group)
mjgrp <-mjgrp[!is.na(mjgrp$major.group), ] 
grp_membership <-left_join(grp_membership, mjgrp)

#structure plot
memprob <-ggplot(grp_membership, aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  #scale_fill_manual(values = tolcolors)+
  theme_cowplot();memprob
#ggsave("bays16_kmeans_membship_150pc_3clust.png", width=10, height=5, path="/Users/tdolan/documents/WIP research/microsats/microsat_figs/figs2023")
#dev.off()
#dim 510x345 gets rid of the white lines.



####################### WITHIN SHINNECOCK ##########################################
##################################################################################

wfpopLD5 <-wfpopLD
setPop(wfpopLD5) <-~Bay
sh <-popsub(wfpopLD5, sublist =c("Shin"))
setPop(sh) <-~Bay/Con/Year

### A priori assignment ###
bay_dapc <- dapc(sh, n.pca = 150, n.da = 2) 
as <- optim.a.score(bay_dapc)#optimum number of pc is 38
bay_dapc <- dapc(sh, n.pca = 38, n.da = 2) 

#PC plots 38pc
scatter(bay_dapc, posi.da="bottomleft",scree.pca=TRUE,posi.pca="bottomright",col=drabcolors,legend=FALSE,solid=1.0) #adjust the solidity value for better plots. 
scatter(bay_dapc,1,1, bg="white", scree.da=FALSE, legend=TRUE, solid=.4,col=drabcolors,)

### Apriori with cross validation, by bay ###
setPop(sh)<-~Bay/Con/Year
df <- genind2df(sh, usepop=TRUE, oneColPerAll = TRUE)
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
X <- scaleGen(sh, NA.method = "mean")
xval <- xvalDapc(X, GRP$pop, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 50, xval.plot = TRUE)
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`); retain 
#redo scatter
bay_dapc <- dapc(sh, n.pca = retain, n.da = 2) 
scatter(bay_dapc, posi.da="topleft",scree.pca=TRUE,posi.pca="topright",legend=FALSE,solid=1.0) #adjust the solidity value for better plots. 

## Compare geographic regions of origin and assigned group memberships
grp_membership <- as.data.frame(bay_dapc$posterior) %>% tibble::rownames_to_column("Ind") %>%
  gather(key = GRP, value = MEMBSHIP, 2:5) %>% arrange(MEMBSHIP)%>%
  left_join(GRP)
#compopolot 
memprob <-ggplot(grp_membership, aes(x = Ind, y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  theme_cowplot() ;memprob
#ggsave("mt_apriori_membship_80pc.png", width=10, height=5, path="/Users/tdolan/documents/WIP research/microsats/microsat_figs/figs2023")

###################### kmeans ###################################
## Let the algo choose k
grp_BIC <- find.clusters.genind(sh, stat = "BIC",choose.n.clust = FALSE,max.n.clust = 40)
plot(grp_BIC$Kstat)

grp_BIC <- find.clusters.genind(sh, stat = "BIC",choose.n.clust = TRUE,max.n.clust = 40)
# the above was really just to look at optimal pc/
# now we specify clusters
grp_BIC <- find.clusters.genind(sh, stat = "BIC",
                                n.clust = 5, n.pca=80)
bay_kstat <-grp_BIC$Kstat

# perform stratified cross validation
X <- scaleGen(sh, NA.method = "mean")
xval <- xvalDapc(X, grp_BIC$grp, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$`Mean Successful Assignment by Number of PCs of PCA`
# identify optimum number of PCs to retain
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)
retain  

#make scatterplot
# perform DAPC using k-mean clusters as groups
dapck <-dapc(sh, grp_BIC$grp, n.pca = retain, n.da=2) #not going to specifiy the number of da, retain 4.
scatter(dapck,grp=grp_BIC$grp, posi.da="topright", bg="white", pch=17:22, scree.pca=TRUE, posi.pca="topleft", col=tolcolors)

table(pop(sh), grp_BIC$grp)

#cluster assignments relative to original group assignments
kgrp <- strata.schemes %>% dplyr::select(Ind, pop)
kgrp$grp <-grp_BIC$grp
# show cluster distribution across populations
kgrps <-kgrp %>% group_by(pop,grp)%>%dplyr::summarize(group=n_distinct(Ind),.groups="drop")

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
mjgrp <-grp_membership%>%group_by(Ind)%>%dplyr::summarize(mxmem =max(MEMBSHIP))
mjgrp<-left_join(grp_membership,mjgrp)%>%
  mutate(major.group=ifelse(mxmem==MEMBSHIP,groupnum,NA))%>%
  dplyr::select(Ind,major.group)
mjgrp <-mjgrp[!is.na(mjgrp$major.group), ] 
grp_membership <-left_join(grp_membership, mjgrp)

#structure plot
memprob <-ggplot(grp_membership, aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  #scale_fill_manual(values = tolcolors)+
  theme_cowplot();memprob
#ggsave("bays16_kmeans_membship_150pc_3clust.png", width=10, height=5, path="/Users/tdolan/documents/WIP research/microsats/microsat_figs/figs2023")
#dev.off()
#dim 510x345 gets rid of the white lines.












######################### WITHIN SHINNECOCK ########################################
######################### ############################################################

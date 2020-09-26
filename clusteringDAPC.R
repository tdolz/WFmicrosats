##### Clustering 

### Created Sept 1, 2020
### From script skreportdnaJAN_28_2020MATTITUCK.R

#install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)
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

#create a dataset that is just 2016 YOY only. 
setPop(wfpopLD) <-~Ocean/Bay/Con/Year
wfyoy16 <-popsub(wfpopLD, blacklist=c("Atl_Mt_3_adults","Atl_Mt_4_adults","Atl_Mt_5_2015","Atl_Mt_5_2016","Atl_Shin_1_2017","Atl_Shin_2_2017","Atl_Mt_1_2015","Atl_Mt_2_2015"))

#shinnecock only database
setPop(wfpopLD) <-~Bay
shinco <-popsub(wfpopLD, sublist=c("Shin"))
setPop(shinco) <-~Bay/Con/Year

#mattituck only database, exclude MT_5
setPop(wfpopLD) <-~Bay/Con
mtco <-popsub(wfpopLD, sublist=c("Mt_1","Mt_2","Mt_3","Mt_4"))
setPop(mtco) <-~Bay/Con/Year

#colorschemes
drabcolors <-c("#d0d1e6","#a6bddb", "#67a9cf", "#1c9099", "#016450")
shincolors <-c("#143601","#68b0ab","#538d22","#aad576")
Mtcolors <-c("#f4d35e","#ee964b","#f95738","#ee4266","#15616d","#0d3b66")
lindsaycolors <-c("#ff5d8f","#ff90b3","#ce4257","#8a2846","#84bcda","#ffd166")
admixcols <-c("#8ecae6","#219ebc","#023047","#ffb703","#fb8500")
sunnycolors <-c("#ec4847","#e6df44","#6600ff", "#061283", "#5bd0c8")
boldcols <-c("firebrick","orange","black","blue","#016450")

library("wesanderson")
royal <-wes_palette(5,name="Darjeeling1",type="continuous")
fox <-wes_palette(5,name="FantasticFox1",type="continuous")
#### Clustering ####
####DAPC USE THIS ONE####
  #http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
  #this is a NEW DAPC, which I am pretty sure is NOT using prior group assignment. 

#optim.a.score <-choosing the right number of PC to retain. 
  # overfitting is very possible and will result in even randomize groups being 100% assigned. 
  # use the A score to optimize the number of pc to choose. 

# Assignplot <- check how well the DAPC is consistent with the original clusters
  #only do a subset of individuals (say 40 at a time) or the plot will be unreadable.
# compoplot <- make a structure-like plot
  #(dapc1, posi="bottomright",
  # a blue cross is the original population, a red shade means high probability according to assingment. 
  # so a blue cross on top of a red box is a good thing. 

#####BAY#####
setPop(wfpopLD) <-~Bay
bay_dapc <- dapc(wfpopLD, n.pca = 150, n.da = 2) 
as <- optim.a.score(bay_dapc)#how many clusters
#ggsave("bayclust_optim17.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs"); dev.off()
#PC plots
scatter(bay_dapc, posi.da="topleft",scree.pca=TRUE,posi.pca="topright",col=boldcols,legend=FALSE,solid=1.0) #adjust the solidity value for better plots. 
#ggsave("bayclust17_150pC_2.png", path="/Users//tdolan//Documents//WIP research//microsats//microsat_figs"); dev.off()
scatter(bay_dapc,1,1, bg="white", scree.da=FALSE, legend=TRUE, solid=.4,col=drabcolors,)
#ggsave("/Users/tdolan/Documents/WIP research/microsats/microsat_figs/bayclust_150pC_1.png");dev.off()
assignplot(bay_dapc)#check assignments
#ggsave("bayclust_assignments.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", height=20, width =5); dev.off()
compoplot(bay_dapc, posi="topright", txt.leg=paste("Cluster", 1:4),ncol=1, xlab="individuals", col=admixcols, show.lab=TRUE)#make a structure-like plot
#ggsave("bayclust_strucplot.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()
temp <- which(apply(bay_dapc$posterior,1, function(e) all(e<0.50)));temp #most admixed individuals? (no more than 50% probability of membership to any group)
compoplot(bay_dapc, subset=temp, show.lab=TRUE, posi="topright", txt.leg=paste("Cluster", 1:4),ncol=2, col=admixcols)
#ggsave("bayclust_admixindplot.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()


######yoy 2016 only#######
# seems to think that the optimal number of pc is 1, which makes no sense. 
setPop(wfyoy16) <-~Bay
yoy_dapc <- dapc(wfyoy16, n.pca = 150, n.da = 2) 
as <- optim.a.score(yoy_dapc)
yoy_dapc <- dapc(wfyoy16, n.pca = 1, n.da = 2) 
#ggsave("yoyclust_optim.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()
scatter(yoy_dapc, posi.da="bottomright",scree.pca=TRUE,posi.pca="bottomleft", possub="bottomleft", legend=TRUE, col=drabcolors)
#ggsave("yoy16clust_150pC_2dc.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()
scatter(yoy_dapc,1,1, bg="white", scree.da=FALSE, legend=TRUE, solid=.4,col=drabcolors,)
#ggsave("yoy16clust_150pC_1.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()
assignplot(yoy_dapc, subset=161:197)##check assignments
#ggsave("yoy16clust_assignments.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", height=20, width =5); dev.off()
compoplot(yoy_dapc, posi="bottomright", txt.leg=paste("Cluster", 1:4), lab="",ncol=1, xlab="individuals", col=drabcolors)#make a structure-like plot
#ggsave("yoy16clust_strucplot.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()
temp <- which(apply(yoy_dapc$posterior,1, function(e) all(e<0.50)));temp #most admixed individuals? (no more than 50% probability of membership to any group)
compoplot(yoy_dapc, subset=temp, show.lab=TRUE, posi="bottomright", txt.leg=paste("Cluster", 1:4),ncol=2, col=drabcolors)
#ggsave("yoy16clust_admixindplot.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()


#######Shin cohorts#######
setPop(shinco) <-~Bay/Con/Year
shin_dapc <- dapc(shinco, n.pca = 70, n.da = 2) 
as <- optim.a.score(shin_dapc)#how many clusters
shin_dapc <- dapc(shinco, n.pca = 43, n.da = 2) 
#ggsave("shinclust_optim.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()
scatter(shin_dapc, posi.da="bottomright",scree.pca=TRUE,posi.pca="bottomleft",solid=1.0, posi.leg="topleft",possub="bottomleft", legend=TRUE, col=shincolors)
#ggsave("shincoClust_70pC_2dc.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()
scatter(shin_dapc,1,1, bg="white", scree.da=FALSE, legend=TRUE, solid=.4,col=shincolors)
#ggsave("shinclust_150pC_1.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()
assignplot(shin_dapc)##check assignments
#ggsave("shinclust_assignments.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", height=20, width =5); dev.off()
compoplot(shin_dapc, posi="bottomright", txt.leg=paste("Cluster", 1:4), lab="",ncol=1, xlab="individuals", col=drabcolors)#make a structure-like plot
#ggsave("shinclust_strucplot.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()
temp <- which(apply(shin_dapc$posterior,1, function(e) all(e<0.50)));temp #most admixed individuals? (no more than 50% probability of membership to any group)
compoplot(shin_dapc, subset=temp, show.lab=TRUE, posi="bottomright", txt.leg=paste("Cluster", 1:4),ncol=2, col=drabcolors)
#ggsave("shinclust_admixindplot.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()


#Mt cohorts - bay con year
setPop(mtco) <-~Bay/Con/Year
mt_dapc <- dapc(mtco, n.pca = 100, n.da = 2) 
as <- optim.a.score(mt_dapc)#how many clusters
mt_dapc <- dapc(mtco, n.pca = 57, n.da = 2) 
#ggsave("mtclust_optim.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()
scatter(mt_dapc, posi.da="topright",scree.pca=TRUE,posi.pca="bottomleft",solid=1.0, posi.leg="topleft",possub="bottomleft", legend=TRUE, col=lindsaycolors)
#ggsave("mtClust_100pC_2dc.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()
scatter(mt_dapc,1,1, bg="white", scree.da=FALSE, legend=TRUE, solid=.6,col=lindsaycolors,posi.pca="topleft")
#ggsave("Mtclust_100pC_1.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()
assignplot(mt_dapc)#check assignments
#ggsave("Mtclust_assignments.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", height=20, width =5); dev.off()
compoplot(mt_dapc, posi="bottomright", txt.leg=paste("Cluster", 1:4), lab="",ncol=1, xlab="individuals", col=lindsaycolors)#make a structure-like plot
#ggsave("Mtclust_strucplot.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()
temp <- which(apply(mt_dapc$posterior,1, function(e) all(e<0.50)));temp #most admixed individuals? (no more than 50% probability of membership to any group)
compoplot(mt_dapc, subset=temp, show.lab=TRUE, posi="bottomright", txt.leg=paste("Cluster", 1:4),ncol=2, col=lindsaycolors)
#ggsave("Mtclust_admixindplot.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()





#need to K fold cross-validate the DAPC. but having a really hard time with this one.... 
wfclust <-find.clusters(wfpopLD, choose.n.clust=TRUE) #kmeans clustering

# You need to edit the original df so that it's identical to wfpopLD, which means manually removing specific individuals and loci. 
wfpop4df <-dplyr::select(wfpop4df,-"WF01.1",-"WF01.2",-"WF06.1",-"WF06.2")
wfpop4df <-filter(wfpop4df, !Ind %in% c(83, 41, 50, 66, 189, 190, 191, 197, 209, 238, 173, 174, 175, 179, 185, 186, 216, 8, 19, 24, 30, 32, 123, 126, 128, 133,
                                        158, 159, 160, 161, 165, 122, 138, 141, 142, 150, "S34", "S39","S74", "S87", "S94"))

#try with the original df? 
#doing this works
x <- dplyr::select(wfpop4df, -Pop)
x$Ind <-rownames(x)
x <-dplyr::select(x, -Ind)
x[is.na(x)] <-"0"
x <-mutate_all(x,as.numeric)

pcs <- seq(2,150,1) #making our vector of pcs to retain
grp <- pop(wfpopLD) #group is assigned population
kgrp <-wfclust$grp # group is done through k-means clustering & 100 pc
dgrp <-bay_dapc$assign # group is done through posterior group assignment after dapc

mat <- base::as.matrix(x)
xval <- xvalDapc.matrix(x, kgrp, training.set = 0.9,n.pca.max = 150,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
#looks overfit to me. 
#how do we get the xlim to be at 0-150, like shannons. 



###### DAPC Shannon's Code #####

# load supplementary functions (https://gist.github.com/sjoleary/88330fec3bf6b7088bc381441be5c7b3)
source("PCA_DAPCfunctions.R")

# genind object of genotypes
# assume strata loaded to group individuals (LIB_ID or SAMPLE_ID) - in this example hierarchical grouping of 
# individuals grouped into populations (POP)
# populations grouped by geographic region (REGION)
# regions grouped by ocean basin (OCEAN)

gen <- wfpopLD #shinco, mtco, wfyoy16
setPop(gen) <- ~Bay

# dataframe with sample information (population designations, hierarchical strata, sex, ...)
# must have column with sample IDs in common to be able to perform joins
setPop(wfpopLD)  <- ~Bay   
df <- genind2df(gen, usepop=TRUE, oneColPerAll = TRUE)
df$Ind <-rownames(df)
df <-df[,c(34,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33)]
df[df=="NA"] <- 0
#make the strata.scheme
wf.strata <-dplyr::select(df, Ind, pop)
strata.schemes <-wf.strata[,c("Ind","pop")] 
row.names(strata.schemes) <-wf.strata$Ind
strata.schemes

# NO ASSUMED POP GEN MODEL/UNSUPERVISED CLUSTERING: K-MEANS CLUSTERING OF INDIVIDUALS + DAPC ----
# k-means-clustering ====
# see documentation to choose K interactively
grp_BIC <- find.clusters.genind(gen, n.pca = 500,
                                stat = "BIC", 
                                #choose.n.clust = FALSE, criterion = "min",
                                choose.n.clust = TRUE,
                                max.n.clust = 40)

# plot BIC per K; chosen value for K highlighted in red
plot.Kstat(grp_BIC, 2)

# determine optimum number of PCs to retain ====
# scale allele frq/missing data replaced with mean
X <- scaleGen(gen, NA.method = "mean")

# perform stratified cross validation
# 90% of individuals used as training individuals, remaining 10% are assigned
xval <- xvalDapc(X, grp_BIC$grp, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)

# identify optimum number of PCs to retain  --> 150
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)

# perform DAPC ====
# perform DAPC using k-mean clusters as groups
dapck <-dapc(gen, grp_BIC$grp, n.pca = retain, n.da = 10)

#compoplot(dapck, posi="bottomright",txt.leg=paste("Cluster", 1:2), lab="", ncol=2, xlab="individuals", col=funky(3))


# evaluate results ====
# coordinates of indv/groups used in scatterplot
DAPC_Ind <- as.data.frame(dapck$ind.coord) %>%
  tibble::rownames_to_column("Ind") %>%
  left_join(strata.schemes)

# add information alloted genetic cluster
Clust <- as.data.frame(grp_BIC$grp) %>%
  tibble::rownames_to_column("Ind") %>%
  dplyr::rename(CLUSTER = `grp_BIC$grp`) %>%
  left_join(DAPC_Ind)

# scatterplot individuals (groups)
scatter.indv <- ggplot(Clust, aes(x = Ind, y = LD1, fill = pop)) +
  geom_point(shape = 21, size = 2) +
  geom_hline(yintercept = 0, color = "darkred", linetype = "dashed", size = 1) +
  labs(x = "Individual", y = "DF 1")+
  facet_grid(. ~ pop, scales = "free", space = "free") +
  theme_cowplot() +
  theme(axis.text.x = element_blank())
#ggsave(scatter.indv, file="scatterindvdapckmeans.png", width = 10, height = 7)

#DF2 by DF1  - there is no LD2 if it's a 2 pc thing, i think. 
df12<-ggplot(DAPC_Ind, aes(x = LD1, y = LD2, fill = pop)) +
  geom_point(size = 2, shape = 21, color = "black", alpha = 1) +
  labs(x = "Discriminant Function 1", y = "Discriminant Function2") +
  scale_fill_manual(values = c("Jam", "Mor", "Mt", "Nap", "Shin")) +
  theme_cowplot()
#ggsave(df12, file="df12kmeans.png", width = 10, height = 7)

#loadingplot
set.seed(4)
contrib <- loadingplot(dapck$var.contr, axis=2,
                       thres=.07, lab.jitter=1)

#DF3 by DF2
df23 <-ggplot(Clust, aes(x = LD2, y = LD3, fill = pop)) +
  geom_point(size = 2, shape = 21, color = "black", alpha = 1) +
  labs(x = "Discriminant Function 2", y = "Discriminant Function3") +
  #scale_fill_manual(values = c("Nap", "Mor", "Jam", "Shin")) +
  theme_cowplot()
#ggsave(df23, file="df23kmeans.png", width = 10, height = 7)

# density plots discriminant functions
densplot <-ggplot(Clust, aes(x = LD1, stat = "position", fill = pop)) +
  geom_density(alpha = .9) +
  labs(x = "Discriminant Function 1", y = "Density") +
  scale_fill_manual(values = drabcolors) +
  theme_cowplot()
#ggsave(densplot, file="densplotkmeans1.png", width = 10, height = 7)

# density plot on LD3
densplot3 <-ggplot(DAPC_Ind, aes(x = LD2, stat = "position", fill = pop)) +
  geom_density(alpha = .9) +
  labs(x = "Discriminant Function 2", y = "Density") +
  #scale_fill_manual(values = col_regs) +
  theme_cowplot()
#ggsave(densplot3, file="densplotkmeans3.png", width = 10, height = 7)

## Compare geographic regions of origin and assigned group memberships
grp_membership <- as.data.frame(dapck$posterior) %>%
  tibble::rownames_to_column("Ind") %>%
  gather(key = GRP, value = MEMBSHIP, 1:2) %>%
  left_join(GRP)
  
  memprob <-ggplot(grp_membership, aes(x = Ind, y = MEMBSHIP, fill = pop)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  theme_cowplot() +
  theme(axis.text.x = element_blank())
#ggsave(memprob, file="membprobpopprior.png", width = 10, height = 7)



#BIPLOTS
scatter(dapck, posi.da="bottomright", bg="white",
        pch=17:22, cstar=0, scree.pca=TRUE,
        posi.pca="bottomleft")



######## A PRIORI GROUPING OF INDIVIDUALS + DAPC ----##############

# use group by region ====
GRP <- strata.schemes %>% dplyr::select(Ind, pop)

# run intial DAPC ====
dapcPop <-dapc(gen, GRP$pop, n.pca = 500, n.da = 6)

temp <- optim.a.score(dapcPop);temp #very flat. #111 pc

#https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html

# cross validation ====
X <- scaleGen(gen, NA.method = "mean")

xval <- xvalDapc(X, GRP$pop, 
                 n.pca.max = 500, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)

retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`); retain #says 200

# run dapc ====
dapcPop <-dapc(gen, GRP$pop, n.pca = retain, n.da = 6)

# evaluate results ====
# coordinates of indv/groups used in scatterplot
DAPC_Ind <- as.data.frame(dapcPop$ind.coord) %>%
  tibble::rownames_to_column("Ind") %>%
  left_join(strata.schemes) %>%
  mutate(pop = ordered(pop))

scatter.indv <- ggplot(DAPC_Ind, aes(x = Ind, y = LD1, fill = pop)) +
  geom_point(shape = 21, size = 2) +
  geom_hline(yintercept = 0, color = "darkred", linetype = "dashed", size = 1) +
  labs(x = "Individual", y = "DF 1")+
  facet_grid(. ~ pop, scales = "free", space = "free") +
  theme_cowplot() +
  theme(axis.text.x = element_blank())
scatter.indv
#ggsave(scatter.indv, file="scatterindvdapcpopprio.png", width = 10, height = 7)


#DF2 by DF1
df12<-ggplot(DAPC_Ind, aes(x = LD1, y = LD2, fill = pop)) +
  geom_point(size = 2, shape = 21, color = "black", alpha = 1) +
  labs(x = "Discriminant Function 1", y = "Discriminant Function2") +
  #scale_fill_manual(values = drabcolors) +
  theme_cowplot(); df12
#ggsave(df12, file="df12popprior.png", width = 10, height = 7)

#DF3 by DF2
df23 <-ggplot(DAPC_Ind, aes(x = LD2, y = LD3, fill = pop)) +
  geom_point(size = 2, shape = 21, color = "black", alpha = 1) +
  labs(x = "Discriminant Function 2", y = "Discriminant Function3") +
  #scale_fill_manual(values = drabcolors) +
  theme_cowplot();df23
#ggsave(df23, file="df23popprior.png", width = 10, height = 7)

# density plots discriminant functions
densplot <-ggplot(DAPC_Ind, aes(x = LD1, stat = "position", fill = pop)) +
  geom_density(alpha = .6) +
  labs(x = "Discriminant Function 1", y = "Density") +
  scale_fill_manual(values = drabcolors) +
  theme_cowplot();densplot
ggsave(densplot, file="densplotpopprior1.png", width = 10, height = 7)

# density plot on LD3
densplot3 <-ggplot(DAPC_Ind, aes(x = LD2, stat = "position", fill = pop)) +
  geom_density(alpha = .6) +
  labs(x = "Discriminant Function 2", y = "Density") +
  #scale_fill_manual(values = col_regs) +
  theme_cowplot();densplot3
ggsave(densplot3, file="densplotpopprior3.png", width = 10, height = 7)

## Compare geographic regions of origin and assigned group memberships
grp_membership <- as.data.frame(dapcPop$posterior) %>%
  tibble::rownames_to_column("Ind") %>%
  gather(key = GRP, value = MEMBSHIP, 2:4) %>%
  left_join(GRP)

memprob <-ggplot(grp_membership, aes(x = Ind, y = MEMBSHIP, fill = pop)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = drabcolors) +
  theme_cowplot() +
  theme(axis.text.x = element_blank());memprob
#ggsave(memprob, file="membprobpopprior.png", width = 10, height = 7)


#compare memembership posterior probabilities
#summary
summary(dapcPop)
round(head(dapcPop$posterior),3)

#assignment plot
assignplot(dapcPop)

compoplot(dapcPop, posi="bottomright",
          txt.leg=paste("Cluster", 1:5), lab="",
          ncol=2, xlab="individuals", col=inferno(5))

#This graph comes out good when you use 200 pc. 
scatter(dapcPop, posi.da="topright", bg="white",
        pch=17:22, cstar=0,  scree.pca=TRUE,
        posi.pca="topleft")

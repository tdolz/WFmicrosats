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
drabcolors2 <-c("darkgrey","#a6bddb", "#67a9cf", "#1c9099", "#283d3b")
shincolors <-c("#143601","#68b0ab","#538d22","#aad576")
shincolors2 <-c("#29bf12","#abff4f","#3f8efc","#3b28cc")
Mtcolors <-c("#f4d35e","#ee964b","#f95738","#ee4266","#15616d","#0d3b66")
Mtcolors2 <-c("#f72585","#7209b7","#4361ee","#4cc9f0")
Mtcolors3 <-c("#f72585","#7209b7","#4361ee","#4cc9f0","#ff930a","#9ef01a")
Mtcolors4 <-c("#fad2e1","#d8bbff","#ff006e","#8338ec","#4361ee","#4cc9f0") 
lindsaycolors <-c("#ff5d8f","#ff90b3","#ce4257","#8a2846","#84bcda","#ffd166")
admixcols <-c("#8ecae6","#219ebc","#023047","#ffb703","#fb8500")
sunnycolors <-c("#ec4847","#e6df44","#6600ff", "#061283", "#5bd0c8","black")
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

#standard error of the mean
sem <- function(x) sd(x)/sqrt(length(x))

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
# not sure if we're interested in this. 
setPop(wfyoy16) <-~Bay
yoy_dapc <- dapc(wfyoy16) #70, 4
as <- optim.a.score(yoy_dapc)#
yoy_dapc <- dapc(wfyoy16, n.pca = 70, n.da = 4) 
#ggsave("yoyclust_optim.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()
scatter(yoy_dapc, posi.da="bottomright",scree.pca=TRUE,posi.pca="bottomleft", possub="bottomleft", legend=TRUE, col=drabcolors)

#######Shin cohorts#######
setPop(shinco) <-~Bay/Con/Year

# must have column with sample IDs in common to be able to perform joins
setPop(shinco)  <- ~Bay/Con/Year   
df <- genind2df(shinco, usepop=TRUE, oneColPerAll = TRUE)
df$Ind <-rownames(df)
df <-df[,c(38,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37)]
df[df=="NA"] <- 0
#make the strata.scheme
wf.strata <-dplyr::select(df, Ind, pop)
strata.schemes <-wf.strata[,c("Ind","pop")] 
row.names(strata.schemes) <-wf.strata$Ind
GRP <- strata.schemes %>% dplyr::select(Ind, pop)

set.seed(1213)
#xval the assignment. 
X <- scaleGen(shinco, NA.method = "mean")
xval <- xvalDapc(X, GRP$pop, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 50, xval.plot = TRUE)
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`); retain 

sdpc <-filter(xval$`Cross-Validation Results`, n.pca==40)
mean(sdpc$success)
sem(sdpc$success)
sd(sdpc$success)

#redo Shinnecock with fewer PC based on the xval? 
shin_dapc <- dapc(shinco, n.pca = 40) #2
scatter(shin_dapc, bg="white", pch=20, cstar=0,posi.leg = "topleft",scree.pca=TRUE,scree.da=TRUE,posi.da="bottomleft", posi.pca="bottomright",col=shincolors2, clab=0, #turns labels on and off
        legend=TRUE, txt.leg = paste(c("YOY early 2016","YOY early 2017","YOY late 2016","YOY late 2017")))
# They are all on top of each other now... 
## Compare geographic regions of origin and assigned group memberships
grp_membership <- as.data.frame(shin_dapc$posterior) %>% tibble::rownames_to_column("Ind") %>%
  gather(key = GRP, value = MEMBSHIP, 2:5) %>% arrange(MEMBSHIP)%>%
  left_join(GRP)
#NEW compopolot - fewer PC - It looks worse.... 
memprob <-ggplot(grp_membership, aes(x = Ind, y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = shincolors2) +theme_cowplot() ;memprob
#to see how many were assigned to their population of origin, compare how many individuals have a GRP that matches POP
boo <-filter(grp_membership, GRP==pop) %>%n_distinct()
toto <-n_distinct(grp_membership$Ind)
boo/toto
boo<-filter(grp_membership, GRP==pop)
mean(boo$MEMBSHIP) #the percent out of total that are assigned to the correct group
sem(boo$MEMBSHIP) #standard error of the mean percentage
# now we have 44% assignment. 
#### using xval to choose how many PC to retain resulted in worse assignment. 


#you also need to do xval for this and naiive clustering. 
# dataframe with sample information (population designations, hierarchical strata, sex, ...)

#assignment validation
## Compare geographic regions of origin and assigned group memberships
grp_membership <- as.data.frame(shin_dapc$posterior) %>%
  tibble::rownames_to_column("Ind") %>%
  gather(key = GRP, value = MEMBSHIP, 2:5) %>%
  arrange(MEMBSHIP)%>%
  left_join(GRP)

#make your own compoplot
memprob <-ggplot(grp_membership, aes(x = Ind, y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = shincolors2) +
  #ylim(0.5,1)+
  theme_cowplot() ;memprob
#theme(axis.text.x = element_blank());memprob
#ggsave(memprob, file="Shincompoplot.png", width = 10, height = 5)

#to see how many were assigned to their population of origin, compare how many individuals have a GRP that matches POP
boo <-filter(grp_membership, GRP==pop) %>%n_distinct()
toto <-n_distinct(grp_membership$Ind)
boo/toto
boo<-filter(grp_membership, GRP==pop)
mean(boo$MEMBSHIP) #the percent out of total that are assigned to the correct group
sem(boo$MEMBSHIP) #standard error of the mean percentage

############k means clustering for shinnecock? 
gen <- shinco # mtco, wfyoy16
setPop(gen) <- ~Bay/Con/Year
set.seed(1713)
grp_BIC <- find.clusters.genind(gen, stat = "BIC",choose.n.clust = TRUE,max.n.clust = 40)

# plot BIC per K; chosen value for K highlighted in red - try 100, and 3
plot.Kstat(grp_BIC, 4)

bays_kstat <-grp_BIC$Kstat

# determine optimum number of PCs to retain ====
# scale allele frq/missing data replaced with mean
X <- scaleGen(gen, NA.method = "mean")

# perform stratified cross validation
# 90% of individuals used as training individuals, remaining 10% are assigned
xval <- xvalDapc(X, grp_BIC$grp, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$`Mean Successful Assignment by Number of PCs of PCA`
# identify optimum number of PCs to retain  --> 150
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)

# perform DAPC ====
# perform DAPC using k-mean clusters as groups
dapck <-dapc(gen, grp_BIC$grp, n.pca = retain) #not going to specifiy the number of da, retain 2.

scatter(dapck, posi.da="bottomright", bg="white", pch=17:22, cstar=0, scree.pca=TRUE, posi.pca="bottomleft")
#clustering plot looks great. 

#the second thing we want to know is how cluster assignments relate to original group assignments
kgrp <- as.data.frame(dapck$grp) %>% tibble::rownames_to_column(var="Ind") #we will use this later. 
# use group by region ====
GRP <- strata.schemes %>% dplyr::select(Ind, pop)

kgrp <-left_join(kgrp,GRP, by="Ind")
# show cluster distribution across populations
kmeanscols <-c("blue","grey","orange","red")
kgrps <-ddply(kgrp, pop~dapck$grp, summarize, group=n_distinct(Ind))
kgrps %>%
  ggplot(aes(x=pop, y=group, fill=`dapck$grp`)) +
  scale_fill_manual(values=kmeanscols)+ylab("number of individuals")+
  geom_bar(stat="identity",alpha=0.7, position=position_dodge())+
  theme(axis.text.x = element_text(angle = 90),axis.text = element_text(size = 14, face="bold"),axis.title = element_text(size = 14,face="bold"),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),panel.spacing = unit(1, "lines"))
#so they're like evenly mixed across groups. 





###### MATTITUCK CREEK ###########
#Mt cohorts - bay con year
# must have column with sample IDs in common to be able to perform joins
setPop(mtco)  <- ~Bay/Con/Year   
df <- genind2df(mtco, usepop=TRUE, oneColPerAll = TRUE)
df$Ind <-rownames(df)
df <-df[,c(38,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37)]
df[df=="NA"] <- 0
#make the strata.scheme
wf.strata <-dplyr::select(df, Ind, pop)
strata.schemes <-wf.strata[,c("Ind","pop")] 
row.names(strata.schemes) <-wf.strata$Ind
GRP <- strata.schemes %>% dplyr::select(Ind, pop)

set.seed(1213)
#xval the assignment. 
X <- scaleGen(mtco, NA.method = "mean")
xval <- xvalDapc(X, GRP$pop, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 50, xval.plot = TRUE)
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`); retain 

sdpc <-filter(xval$`Cross-Validation Results`, n.pca==40)
mean(sdpc$success)
sem(sdpc$success)
sd(sdpc$success)

#now DAPC
mt_dapc <- dapc(mtco) 
mt_dapc <- dapc(mtco, n.pca = 40) #retaining 3 da
#ggsave("mtclust_optim.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()
scatter(mt_dapc,bg="white", pch=20, cstar=0, scree.pca=TRUE,scree.da=TRUE, posi.da="bottomright",posi.pca="bottomleft",col=Mtcolors4, clab=0, #turns labels on and off
        txt.leg = paste(c("YOYe 2015","YOYe 2016", "YOYl 2015","YOYl 2016","Migrant Adults","Resident Adults")),legend=TRUE)
#ggsave("mtClust_bayconyear.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs");dev.off()


#txt.leg = paste(c("YOYe 2015","YOYe 2016", "YOYl 2015","YOYl 2016","Migrant Adults","Resident Adults"))

#assignment validation
## Compare geographic regions of origin and assigned group memberships
grp_membership <- as.data.frame(mt_dapc$posterior) %>%
  tibble::rownames_to_column("Ind") %>%
  gather(key = GRP, value = MEMBSHIP, 2:7) %>%
  left_join(GRP)

#make your own compoplot
memprob <-ggplot(grp_membership, aes(x = Ind, y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = Mtcolors4) +
  #ylim(0.5,1)+
  theme_cowplot();memprob
#theme(axis.text.x = element_blank());memprob
ggsave(memprob, file="mtcompoplot.png", width = 10, height = 5)

#to see how many were assigned to their population of origin, compare how many individuals have a GRP that matches POP
boo <-filter(grp_membership, GRP==pop) %>%n_distinct()
toto <-n_distinct(grp_membership$Ind)
boo/toto
boo<-filter(grp_membership, GRP==pop)
mean(boo$MEMBSHIP) #61%
sem(boo$MEMBSHIP)
plyr::ddply(boo,~pop,count)
plyr::ddply(grp_membership,~pop, summarize, membership=mean(MEMBSHIP))
plyr::ddply(boo,~pop, summarize, membership=mean(MEMBSHIP))

#xval
X <- scaleGen(mtco, NA.method = "mean")
xval <- xvalDapc(X, GRP$pop, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`); retain #50
### In this case the xval actually reccommended more PC than we previously used. 


#naive clustering
#I can't get naive clustering to work for MT. 
gen <- mtco # mtco, wfyoy16
setPop(gen) <- ~Bay/Con/Year
set.seed(1713)
grp_BIC <- find.clusters.genind(gen, stat = "BIC",choose.n.clust = TRUE,max.n.clust = 40) #40, 2

# plot BIC per K; chosen value for K highlighted in red - try 100, and 3
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
xval$`Mean Successful Assignment by Number of PCs of PCA`
# identify optimum number of PCs to retain  --> 150
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)

# perform DAPC ====
# perform DAPC using k-mean clusters as groups
dapck <-dapc(gen, grp_BIC$grp, n.pca = retain) #not going to specifiy the number of da, retain 2.

scatter(dapck, posi.da="bottomright", bg="white", pch=17:22, cstar=0, scree.pca=TRUE, posi.pca="bottomleft")
#clustering plot looks great. 

#the second thing we want to know is how cluster assignments relate to original group assignments
kgrp <- as.data.frame(dapck$grp) %>% tibble::rownames_to_column(var="Ind") #we will use this later. 
# use group by region ====
GRP <- strata.schemes %>% dplyr::select(Ind, pop)

kgrp <-left_join(kgrp,GRP, by="Ind")
# show cluster distribution across populations
kmeanscols <-c("blue","red")
kgrps <-ddply(kgrp, pop~dapck$grp, summarize, group=n_distinct(Ind))
kgrps %>%
  ggplot(aes(x=pop, y=group, fill=`dapck$grp`)) +
  scale_fill_manual(values=kmeanscols)+ylab("number of individuals")+
  geom_bar(stat="identity",alpha=0.7, position=position_dodge())+
  theme(axis.text.x = element_text(angle = 90),axis.text = element_text(size = 14, face="bold"),axis.title = element_text(size = 14,face="bold"),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),panel.spacing = unit(1, "lines"))
#so they're like evenly mixed across groups. 












###### DAPC Shannon's Code #####

# load supplementary functions (https://gist.github.com/sjoleary/88330fec3bf6b7088bc381441be5c7b3)
source("PCA_DAPCfunctions.R")

#####Edit shannon's gist because some of the code is deprecated ####
plot.Kstat <- function(grp, k_clust){
  Kstat <- as.data.frame(grp$Kstat) %>%
    tibble::rownames_to_column("K") %>%
    separate(K, c("temp", "K"), sep = "=") %>%
    rename(BIC = `grp$Kstat`) %>%
    select(K, BIC)
  
  Kstat$K <- as.numeric(Kstat$K)
  
  ggplot(Kstat, aes(x = K, y = BIC)) +
    geom_line() +
    geom_point(shape = 21, color = "black", fill = "dark grey", size = 4) + 
    geom_point(data = subset(Kstat, K == k_clust), shape = 19, size = 4, color = "red") +
    labs(x = "number of clusters K", y = "BIC") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 16),
      axis.title.y = element_text(vjust = 1.5),
      
      legend.position = "bottom",
      
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 16))
}
## ------------------------------------------------------------------------- ##


## ------------------------------------------------------------------------- ##
# ---------------------- Extract Eigenvalues  ------------------------------- #

# returns dataframe with Eigenvalues and %variance explained by each PC
# argument object class PCA (e.g. generated using adegenet)

eigenvalues <- function(PCA){
  eig <- data.frame(PCA$eig) %>%
    rownames_to_column("PC") %>%
    rename(Eigenvalue = `PCA.eig`) %>%
    mutate(Percent = Eigenvalue/sum(Eigenvalue))
  eig$PC <- as.numeric(eig$PC)
  eig
}
## ------------------------------------------------------------------------- ##


## ------------------------------------------------------------------------- ##
# --------------------------  plot Eigenvalues ------------------------------ #

# plots Eigenvalues of top 50 PCs
# argument data frame with eigenvalues and % variance

plot.eigen <- function(eig){
  ggplot(eig[1:50, ], aes(x = PC, y = Eigenvalue)) +
    geom_bar(stat = "identity", color = "black", fill = "darkorange") +
    labs(x = "Principle Component", y = "Eigenvalue") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 16),
      axis.title.y = element_text(vjust = 1.5),
      
      legend.position = "bottom",
      
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 16))
}
## ------------------------------------------------------------------------- ##


## ------------------------------------------------------------------------- ##
# -----------------  plot %variance explained by PCs ------------------------ #

# plots %variance explained by of top 25 PCs
# argument data frame with eigenvalues and % variance

plot.eigen.variance <- function(eig){
  ggplot(eig[1:25, ], aes(x = PC, y = Percent)) +
    geom_bar(stat = "identity", color = "black", fill = "darkorange") +
    labs(x = "Principle Component", y = "% Variance") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 16),
      axis.title.y = element_text(vjust = 1.5),
      
      legend.position = "bottom",
      
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 16))
}
## ------------------------------------------------------------------------- ##


## ------------------------------------------------------------------------- ##
# ---------------------- PC loadings by individual -------------------------- #


# Individuals' contribution to PCs/calculate Loading per individual and PC
# names samples ID column LIB_ID - can change to Sample ID to be able to better
# join with SampleInfo data frame if needed

PC.ind <- function(PCA){
  PCA$li %>%
    rownames_to_column("LIB_ID") %>%
    mutate(Loading1 = Axis1^2) %>%
    mutate(Loading2 = Axis2^2) %>%
    mutate(Loading3 = Axis3^2)
}

################





######### NOW SHANNON"S WAY ###########

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
df <-df[,c(38,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37)]
df[df=="NA"] <- 0
#make the strata.scheme
wf.strata <-dplyr::select(df, Ind, pop)
strata.schemes <-wf.strata[,c("Ind","pop")] 
row.names(strata.schemes) <-wf.strata$Ind
strata.schemes

# NO ASSUMED POP GEN MODEL/UNSUPERVISED CLUSTERING: K-MEANS CLUSTERING OF INDIVIDUALS + DAPC ----
# k-means-clustering ====
# see documentation to choose K interactively
set.seed(1713)
grp_BIC <- find.clusters.genind(gen, 
                                #n.pca = 100,
                                stat = "BIC", 
                                #choose.n.clust = FALSE, criterion = "min",
                                choose.n.clust = TRUE,
                                max.n.clust = 40)

# plot BIC per K; chosen value for K highlighted in red - try 100, and 3
plot.Kstat(grp_BIC, 4)

bays_kstat <-grp_BIC$Kstat

# determine optimum number of PCs to retain ====
# scale allele frq/missing data replaced with mean
X <- scaleGen(gen, NA.method = "mean")

# perform stratified cross validation
# 90% of individuals used as training individuals, remaining 10% are assigned
xval <- xvalDapc(X, grp_BIC$grp, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$`Mean Successful Assignment by Number of PCs of PCA`
# identify optimum number of PCs to retain  --> 150
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)

# perform DAPC ====
# perform DAPC using k-mean clusters as groups
dapck <-dapc(gen, grp_BIC$grp, n.pca = retain) #not going to specifiy the number of da, retain 2.

scatter(dapck, posi.da="topright", bg="white", pch=17:22, cstar=0, scree.pca=TRUE, posi.pca="topleft")
#clustering plot looks great. 

#the second thing we want to know is how cluster assignments relate to original group assignments
kgrp <- as.data.frame(dapck$grp) %>% tibble::rownames_to_column(var="Ind") #we will use this later. 
# use group by region ====
GRP <- strata.schemes %>% dplyr::select(Ind, pop)

kgrp <-left_join(kgrp,GRP, by="Ind")
# show cluster distribution across populations
kmeanscols <-c("blue","grey","orange","red")
kgrps <-ddply(kgrp, pop~dapck$grp, summarize, group=n_distinct(Ind))
kgrps %>%
  ggplot(aes(x=pop, y=group, fill=`dapck$grp`)) +
  scale_fill_manual(values=kmeanscols)+ylab("number of individuals")+
  geom_bar(stat="identity",alpha=0.7, position=position_dodge())+
  theme(axis.text = element_text(size = 14, face="bold"),axis.title = element_text(size = 14,face="bold"),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),panel.spacing = unit(1, "lines"))
#so they're like evenly mixed across groups. 

#ggsave("mortalitybarplot_drab.png", path="/Users/tdolan/Documents/Proposal/Figures")
#dev.off()
######## A PRIORI GROUPING OF INDIVIDUALS + DAPC ----##############

#https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html

# use group by region ====
GRP <- strata.schemes %>% dplyr::select(Ind, pop)

# run intial DAPC ====
set.seed(4)

# cross validation ====
X <- scaleGen(gen, NA.method = "mean")

xval <- xvalDapc(X, GRP$pop, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)

retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`); retain #says 150

sdpc <-filter(xval$`Cross-Validation Results`, n.pca==retain)
mean(sdpc$success)
sem(sdpc$success)
sd(sdpc$success)

dapcPop <-dapc(gen, GRP$pop, n.pca = retain, n.da =4, var.loadings = TRUE)

#BIPLOTS
scatter(dapcPop, posi.da="bottomright", bg="white", pch=20, cstar=0, scree.pca=TRUE, posi.pca="bottomleft",col=drabcolors2, #clab=0, #turns labels on and off
        legend=TRUE, txt.leg = paste(c("Jamaica","Moriches","Mattituck","Napeague","Shinnecock")))


#loadingplot
set.seed(4)
contrib <- loadingplot(dapcPop$var.contr, axis=2,threshold = 0.02,lab.jitter=2)

loads <-as.data.frame(dapcPop$var.contr) %>% rownames_to_column()%>% mutate(totalLD=LD1+LD2+LD3+LD4) %>% separate(rowname, c("locus","allele"), remove=FALSE)
loads <-arrange(loads,desc(totalLD),locus)

# evaluate results ====
# coordinates of indv/groups used in scatterplot
DAPC_Ind <- as.data.frame(dapcPop$ind.coord) %>%
  tibble::rownames_to_column("Ind") %>%
  left_join(strata.schemes) %>%
  mutate(pop = ordered(pop))

#DF2 by DF1
df12<-ggplot(DAPC_Ind, aes(x = LD1, y = LD2, fill = pop)) +
  geom_point(size = 2, shape = 21, color = "black", alpha = 1) +
  labs(x = "Discriminant Function 1", y = "Discriminant Function2") +
  scale_fill_manual(values = drabcolors) +
  theme_cowplot(); df12
#ggsave(df12, file="df12popprior.png", width = 10, height = 7)

#DF3 by DF2
df23 <-ggplot(DAPC_Ind, aes(x = LD2, y = LD3, fill = pop)) +
  geom_point(size = 2, shape = 21, color = "black", alpha = 1) +
  labs(x = "Discriminant Function 2", y = "Discriminant Function3") +
  scale_fill_manual(values = drabcolors) +
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
  scale_fill_manual(values = drabcolors) +
  theme_cowplot();densplot3
ggsave(densplot3, file="densplotpopprior3.png", width = 10, height = 7)


## Compare geographic regions of origin and assigned group memberships
grp_membership <- as.data.frame(dapcPop$posterior) %>%
  tibble::rownames_to_column("Ind") %>%
  gather(key = GRP, value = MEMBSHIP, 2:6) %>%
  left_join(GRP)

#make your own compoplot
grp_membership <-arrange(grp_membership, MEMBSHIP)
memprob <-ggplot(grp_membership, aes(x = Ind, y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = drabcolors) +
  #ylim(0.5,1)+
  theme_cowplot() +
  theme(axis.text.x = element_blank());memprob
#ggsave(memprob, file="membprobpopprior.png", width = 10, height = 7)

#to see how many were assigned to their population of origin, compare how many individuals have a GRP that matches POP
boo <-filter(grp_membership, GRP==pop) %>%n_distinct()
toto <-n_distinct(grp_membership$Ind)
boo/toto
boo<-filter(grp_membership, GRP==pop)
plyr::ddply(boo,~pop,count)
plyr::ddply(grp_membership,~pop, summarize, membership=mean(MEMBSHIP))
plyr::ddply(boo,~pop, summarize, membership=mean(MEMBSHIP))
mean(boo$MEMBSHIP)
sem(boo$MEMBSHIP)



#compare memembership posterior probabilities
#summary
summary(dapcPop)
round(head(dapcPop$posterior),3)












#Created 3/6/24 to investigate why certain loci are controlling the loadings. 
#notably WF12 and PSY022
#
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
removeloc <- c("WF06", "WF12","PSY022")# create vector containing loci to remove
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

## A priori assignment ###
drabcolors <-c("#d0d1e6","#a6bddb", "#67a9cf", "#1c9099", "#016450")

### A priori assignment ###
bay_dapc <- dapc(wfbays16, n.pca = 150, n.da = 2) 
as <- optim.a.score(bay_dapc)#optimum number of pc is 62
bay_dapc <- dapc(wfbays16, n.pca = 62, n.da = 2) 

#PC plots 62pc
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

#RMSE of PC
xval$`Root Mean Squared Error by Number of PCs of PCA`
# best # of PC
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`); retain 
#100pc
xval$`Mean Successful Assignment by Number of PCs of PCA`

#re-do plots based on xval
bay_dapc <- dapc(wfbays16, n.pca = retain, n.da = 2) 
scatter(bay_dapc, posi.da="topleft",scree.pca=TRUE,posi.pca="topright",col=drabcolors,legend=FALSE,solid=1.0)
#ggsave('yoy16apriori_scatter.png', path="./DAPC/DAPC_figs", width = 4.5, height = 4)
#adjust the solidity value for better plots. 

#LOADING PLOT - diverging alleles between a priori groupings. 
set.seed(4)
contrib <- loadingplot(bay_dapc$var.contr, axis=2,
                       thres=.02, lab.jitter=1)
#save 700x500

# WF22, Pam27, WF517, WF32
#i am not sure I understand why the loci have the subscript. like why Pam27.180 vs. just Pam27?

#PLOT OF allele divergence. 
par(mfrow=c(2,2),  mai = c(0.4, 0.4, 0.4, 0.4))

#PAM 27
freqPam27 <- tab(genind2genpop(wfbays16[loc=c("PAM27")]),freq=TRUE)
#of the Pam27 alleles that are most divergent, are PAM27.180 and PAM27.160 
freqPam27.2 <-freqPam27[,c(2,5)]
matplot(freqPam27.2,pch=c("A","B"), type="b",
        xlab="pop",ylab="allele frequency", xaxt="n",
        cex=1, main="Pam27")
axis(side=1, at=1:5, lab=c("Jam","Mor","Mt","Nap","Shin"))

#WF517
freqWF517 <- tab(genind2genpop(wfbays16[loc=c("WF517")]),freq=TRUE)
#of the Pam27 alleles that are most divergent, are PAM27.180 and PAM27.160 
freqWF517.281 <-freqWF517[,1]
matplot(freqWF517.281,pch=c("A"), type="b",
        xlab="pop",ylab="allele frequency", xaxt="n",
        cex=1, main="WF517")
axis(side=1, at=1:5, lab=c("Jam","Mor","Mt","Nap","Shin"))

#WF22
freqWF22 <- tab(genind2genpop(wfbays16[loc=c("WF22")]),freq=TRUE)
#of the Pam27 alleles that are most divergent, are PAM27.180 and PAM27.160 
freqWF22.191 <-freqWF22[,"WF22.191"]
matplot(freqWF22.191,pch=c("A"), type="b",
        xlab="pop",ylab="allele frequency", xaxt="n",
        cex=1, main="WF22")
axis(side=1, at=1:5, lab=c("Jam","Mor","Mt","Nap","Shin"))

#WF32
freqWF32 <- tab(genind2genpop(wfbays16[loc=c("WF32")]),freq=TRUE)
#of the Pam27 alleles that are most divergent, are PAM27.180 and PAM27.160 
freqWF32.213 <-freqWF32[,"WF32.213"]
matplot(freqWF32.213,pch=c("A"), type="b",
        xlab="pop",ylab="allele frequency", xaxt="n",
        cex=1, main="WF32")
axis(side=1, at=1:5, lab=c("Jam","Mor","Mt","Nap","Shin"))



#Make some kind of cowplot or overlay of all of these. 
# Then, look at the alleles that are driving divergence in the k-means clustering for the n clusters that is best indicated
# by BIC and by the cluster tuning (that is k=2). Are the same alleles driving that divergence? Are the frequency/loadings of 
# of the same magnitude. 

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
ggsave("bays16_apriori_membship_100pc.png", width=10, height=5, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

#SAVE grp_membership under new name for later. 
grp_membership_ap <- mutate(grp_membership, nclust="ap")

#mean membership probability for correct assginment. 
filter(grp_membership, GRP == pop)%>%group_by(pop)%>%summarize(mm=mean(MEMBSHIP))

#mean membership probability for incorrect assignment
filter(grp_membership, GRP != pop)%>%group_by(pop)%>%summarize(mm=mean(MEMBSHIP))

#total percentage of each population colored by each group. 
filter(grp_membership, GRP == pop)%>%group_by(pop)%>%summarize(tA=sum(MEMBSHIP), n_ind=n_distinct(Ind))


########################################### kmeans ###################################################
########################################### 

#changed baycolors to tolcolors
#baycolors <-c("#6932a8","#02edf5","#4b4a5c", "#f78605","#f7ef05", "#2958c4", "#48bf11","#940d0d")
baycolors <-c("#332288","#88CCEE","#44AA99","#999933","#DDCC77","#CC6677","#882255","#AA4499")

## Check that a starting value of 100 pc is reasonable. 
grp_BIC <- find.clusters.genind(wfbays16, stat = "BIC",choose.n.clust = FALSE,max.n.clust = 40)
plot(grp_BIC$Kstat)

#loop to find the best k at different numbers of PC. 
bestk <- c()
k.bic <-c()
pc_num <- seq(5,200,5)
for (i in 1:length(pc_num)){
  grp.bic <-find.clusters.genind(wfbays16, stat = "BIC",choose.n.clust = FALSE,max.n.clust = 20,n.pca=pc_num[i])
  k.bic[i]<-min(grp.bic$Kstat)
  bestk[i]<-match(min(grp.bic$Kstat),grp.bic$Kstat)
}
pltk <- cbind(pc_num,bestk,k.bic)%>% as.data.frame()
pltk %>% ggplot(aes(x=pc_num,y=bestk))+geom_line()+
  xlab("Number of PC retained")+ylab("Number of clusters")+
  ggtitle("Optimal number of clusters as determined by BIC")+
  theme_classic()
ggsave("yoy16_KvsNPC.png", width=5, height=5, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")


#CHOOSE K
#start with 100 pc
grp_BIC <- find.clusters.genind(wfbays16, stat = "BIC",
                                n.clust = 4, #CHOOOSEEEEEE
                                n.pca=100)

# # perform stratified cross validation
X <- scaleGen(wfbays16, NA.method = "mean")
xval <- xvalDapc(X, grp_BIC$grp, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$`Mean Successful Assignment by Number of PCs of PCA`
# identify optimum number of PCs to retain  --> 100, mean successful assignment 
xval$`Root Mean Squared Error by Number of PCs of PCA`
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`)
retain 

#make scatterplot
# perform DAPC using k-mean clusters as groups
dapck <-dapc(wfbays16, grp_BIC$grp, n.pca = retain, n.da=2) #not going to specifiy the number of da, retain 4.
scatter(dapck,grp=grp_BIC$grp, posi.da="topleft", bg="white", pch=17:22, scree.pca=TRUE, posi.pca="bottomleft", col=baycolors)


#LOADING PLOT - diverging alleles between a priori groupings. 
set.seed(4)
contrib <- loadingplot(dapck$var.contr,# axis=2,
                       thres=.05, lab.jitter=1)
# Different loci contributing: WF12 and PSY022. 

tabs <-table(pop(wfbays16), grp_BIC$grp)
#CHANGE THE NUMBER
t8 <- as.data.frame(tabs)%>%dplyr::rename("Bay"="Var1","Clust"="Var2","n_indv"="Freq")

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
grp_membership <-left_join(grp_membership, mjgrp)%>%
  mutate(groupnum=as.factor(groupnum))

#structure plot
#groupnum is the filling, so that you have both the filling associated with the 
memprob <-ggplot(grp_membership, aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot();memprob
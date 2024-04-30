#Created April 30, 2024 from the script DAPC_yoy16.R to experiment with eliminating loci. 


### Here, in this sheet we are trying to do the kmeans tuning. 
### in a second sheet we willcompare bays based on 2016 YOY alone
### in a third sheet, we will compare within Mattituck and Within Shinnecock. 

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
wfpop <- read.genalex("/Users//taradolan/Documents//R-Github//WFmicrosats/data/popcorrect_17_sept20204genalex_doubl0ABC.csv")
wfpop4df <-read.csv("/Users//taradolan/Documents//R-Github//WFmicrosats/data/popcorrect_17_sept2020_doubl0ABC.csv", header = TRUE) #csv version 

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
#bay_dapc <- dapc(wfbays16, n.pca = 150, n.da = 2) 
#as <- optim.a.score(bay_dapc)#optimum number of pc is 62
#bay_dapc <- dapc(wfbays16, n.pca = 62, n.da = 2) 

#PC plots 62pc
#scatter(bay_dapc, posi.da="bottomleft",scree.pca=TRUE,posi.pca="bottomright",col=drabcolors,legend=FALSE,solid=1.0) #adjust the solidity value for better plots. 
#scatter(bay_dapc,1,1, bg="white", scree.da=FALSE, legend=TRUE, solid=.4,col=drabcolors,)

### Apriori with cross validation, by bay ###
setPop(wfbays16)<-~Bay
df <- genind2df(wfbays16, usepop=TRUE, oneColPerAll = TRUE)
df$Ind <-rownames(df)
df <-df[,c(ncol(df),1:ncol(df)-1)]
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

scatter(bay_dapc, posi.da="none", bg="white",
        #pch=17:22, 
        pch=19,
        cstar=0, cellipse = 1.5,
        col=drabcolors, scree.pca=TRUE,
        posi.pca="none",clabel=1)


################LOADINGS ##### see script "loading_plots.R" for visualizations 
#write.csv(bay_dapc$var.contr,file="./DAPC/DAPC_figs/bay16_loadings_apriori.csv")

## Compare geographic regions of origin and assigned group memberships
grp_membership <- as.data.frame(bay_dapc$posterior) %>% tibble::rownames_to_column("Ind") %>%
  gather(key = GRP, value = MEMBSHIP, 2:6) %>% arrange(MEMBSHIP)%>%
  left_join(GRP)
#NEW compopolot - fewer PC - It looks worse.... 
memprob <-ggplot(grp_membership, aes(x = fct_reorder(Ind,MEMBSHIP), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = drabcolors) +theme_cowplot()+
  theme(
    strip.background =element_rect(fill="transparent")
  );memprob
#ggsave("bays16_apriori_membship_100pc.png", width=10, height=5, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

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

#Make a new BIC plot
BIC_100PC <- as.data.frame(grp_BIC$Kstat)
BIC_100PC$k<-seq(1,40,1)
BIC_100PC %>% filter(k < 21)%>%
  ggplot(aes(x=k,y=`grp_BIC$Kstat`))+geom_point()+
  xlab("Number of clusters (k)")+ylab("BIC")+
  ggtitle("Optimal number of clusters as determined by BIC")+
  theme_classic()
#ggsave("yoy16_KvsBIC.png", width=5, height=5, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")


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
#ggsave("yoy16_KvsNPC.png", width=5, height=5, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")


#CHOOSE K
#start with 100 pc
grp_BIC <- find.clusters.genind(wfbays16, stat = "BIC",
                                n.clust = 2, #CHOOOSEEEEEE
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
scatter(dapck,grp=grp_BIC$grp, posi.da="none", bg="white", pch=17:22, scree.pca=TRUE, posi.pca="none", col=baycolors)



################LOADINGS ##### see script "loading_plots.R" for visualizations 
ldvar <-dapck$var.contr
# 
#rownames(ldvar)<-ldvar[,1]

#summarize loadings
lds <- ldvar %>% as.data.frame()%>%rownames_to_column(var="Allele")%>%
  separate(col=Allele, into=c("locus","allele"), sep="[.]",remove=FALSE)%>%mutate(locus=as.factor(locus))

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
ldvar <-as.data.frame(ldvar)#%>%select(-X)
set.seed(4)
#par(mfrow=c(1,2))
#AXIS 1
contrib1 <- loadingplot(ldvar, axis=1,
                        thres=thresh1, lab.jitter=4, main="Axis 1", xlab="Alleles", ylab="DAPC Loadings")

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
  theme_cowplot()+
  theme(
    strip.background =element_rect(fill="transparent")
  );memprob
ggsave("bays16_kmeans_membship_60pc_8clust.png", width=10, height=5, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

#dev.off()
#dim 510x345 gets rid of the white lines.

#mean membership probability. 
mjgrp <-grp_membership%>%group_by(Ind)%>%dplyr::summarize(mxmem =max(MEMBSHIP))
mjgrp<-left_join(grp_membership,mjgrp)%>%
  mutate(major.group=ifelse(mxmem==MEMBSHIP,groupnum,NA))
mjgrp <-mjgrp[!is.na(mjgrp$major.group), ] 
mjgrp%>%group_by(pop)%>%summarize(mm=mean(MEMBSHIP))
mean(mjgrp$MEMBSHIP)

#total percentage of each population colored by each group. 
#grp_membership%>%group_by(pop)%>%summarize(tA=sum(MEMBSHIP), n_ind=n_distinct(Ind))

#SAVE each version of grp membership
#grp_membership_k2 <-mutate(grp_membership,nclust=2)%>%arrange(major.group,Ind)
#grp_membership_k3 <-mutate(grp_membership,nclust=3)%>%arrange(major.group,Ind)
#grp_membership_k4 <-mutate(grp_membership,nclust=4)%>%arrange(major.group,Ind)
#grp_membership_k5 <-mutate(grp_membership,nclust=5)%>%arrange(major.group,Ind)
#grp_membership_k6 <-mutate(grp_membership,nclust=6)%>%arrange(major.group,Ind)
#grp_membership_k7 <-mutate(grp_membership,nclust=7)%>%arrange(major.group,Ind)
#grp_membership_k8 <-mutate(grp_membership,nclust=8)%>%arrange(major.group,Ind)


grp_membership_all <-bind_rows(grp_membership_k2,grp_membership_k3,grp_membership_k4,
                               grp_membership_k5,grp_membership_k6,grp_membership_k7,grp_membership_k8)%>%
  mutate(nclust=as.factor(nclust))%>%bind_rows(grp_membership_ap)
tabs_all <-bind_rows(t2,t3,t4,t5,t6,t7,t8)
write.csv(tabs_all, file="./DAPC/DAPC_figs/yoy16_indbyclust.csv")
write.csv(grp_membership_all, file="./DAPC/DAPC_figs/yoy16dapc.csv")

## Figure out why plots designated by group number are different than plots designated by major group. 
#because groupnum indicates the color and major.group indicates the dominant group. when you're making a bar that 
#is part one color and part another

#I think the table might be the best way to present this. 

#try to reorder the grp membership
grp_membership_all2 <-grp_membership_all %>% 
  mutate(Ind = as.numeric(Ind))%>%
  arrange(major.group,Ind)%>%
  mutate(Ind = as.factor(Ind))

grp_membership_all2%>%
  filter(nclust !="ap")%>%
  arrange(major.group,Ind,MEMBSHIP)%>%
ggplot(aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size=9),
    strip.background =element_rect(fill="transparent"),
    #strip.text = element_text(size=9),
    legend.position = "none",
    panel.spacing = unit(1.5, "mm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = "dark grey"),
    plot.background = element_rect(fill = NA,colour = "black", size=0.5),
    panel.border = element_rect(colour="black",fill=NA, size=0.3)
  )
ggsave("yoy16_allMembship.png", width=8, height=10, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

# you actually cannot do this with facet_grid because with different numbers of k, 
# an individual may be assigned to different major groups! 

#might be possible to hacky workaround with cowplot.Separately graph each row, but this would 
#reorder the individuals in every row.  
# HOWEVER, it is kind of important to show that the same individuals in the same order 
# are assigned to different groups. 



####################### WITHIN MATTITUCK ##########################################
##################################################################################

wfpopLD4 <-wfpopLD
setPop(wfpopLD4) <-~Bay
mt <-popsub(wfpopLD4, sublist =c("Mt"))
setPop(mt) <-~Bay/Con/Year
#exclude unidentified adults but we will include the 2015 YOY for now. 
mt <-popsub(mt,exclude=c("Mt_5_2015", "Mt_5_2016"))
#setPop(mt) <-~Con/Year

### A priori assignment ###
bay_dapc <- dapc(mt, n.pca = 100, n.da = 2) 
as <- optim.a.score(bay_dapc)#optimum number of pc is 1
bay_dapc <- dapc(mt, n.pca = 52, n.da = 2) 

#PC plots 
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
#change the population names here. 
strata.schemes <-mutate(strata.schemes,pop=fct_recode(pop, 
                "Early 2015"="Mt_1_2015","Early 2016"="Mt_1_2016", 
                "Late 2015"="Mt_2_2015","Late 2016"="Mt_2_2016","Migrants"="Mt_3_adults","Residents"="Mt_4_adults"))
row.names(strata.schemes) <-wf.strata$Ind
GRP <- strata.schemes %>% dplyr::select(Ind, pop)

set.seed(1213)
#xval the assignment. 
X <- scaleGen(mt, NA.method = "mean")
xval <- xvalDapc(X, GRP$pop, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`); retain 
xval$`Mean Successful Assignment by Number of PCs of PCA`
xval$`Root Mean Squared Error by Number of PCs of PCA`

#Mtcolors4 <-c("#fb9a99","#66c2a5","#ff006e","#8338ec","#4361ee","#4cc9f0") 
Mtcolors4 <-c("#fad2e1","#d8bbff","#ff006e","#8338ec","#4361ee","#4cc9f0") 

#redo scatter
bay_dapc <- dapc(mt, n.pca = retain, n.da = 2) 
scatter(bay_dapc, posi.da="topleft",scree.pca=TRUE,posi.pca="topright",legend=FALSE,solid=1.0,col=Mtcolors4) #adjust the solidity value for better plots. 
scatter(bay_dapc, label=c("Early 2015","Early 2016","Late 2015","Late 2016","Migrant","Resident"), posi.da="none",scree.pca=TRUE,posi.pca="none",legend=FALSE,solid=1.0,col=Mtcolors4) #adjust the solidity value for better plots. 

scatter(bay_dapc, posi.da="none", bg="white",
        #pch=17:22, 
        pch=19,
        cstar=0, cellipse = 1.5,
        col=Mtcolors4, scree.pca=TRUE,
        posi.pca="none",clabel=0.7, label=c("Early 2015","Early 2016","Late 2015","Late 2016","Migrant","Resident"))


################LOADINGS ##### see script "loading_plots.R" for visualizations 
write.csv(bay_dapc$var.contr,file="./DAPC/DAPC_figs/MT15_loadings_apriori.csv")


## Compare geographic regions of origin and assigned group memberships
grp_membership <- as.data.frame(bay_dapc$posterior) %>% tibble::rownames_to_column("Ind") %>%
  gather(key = GRP, value = MEMBSHIP, 2:7) %>% arrange(MEMBSHIP)%>%
  left_join(GRP)

mjgrp <-grp_membership%>%group_by(Ind)%>%dplyr::summarize(mxmem =max(MEMBSHIP))
mjgrp<-left_join(grp_membership,mjgrp)%>%
  mutate(major.group=ifelse(mxmem==MEMBSHIP,pop,NA))%>%
  dplyr::select(Ind,major.group)
mjgrp <-mjgrp[!is.na(mjgrp$major.group), ] 
grp_membership <-left_join(grp_membership, mjgrp)

#compopolot 
memprob <-ggplot(grp_membership, aes(x = fct_reorder(Ind,MEMBSHIP), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  scale_fill_manual(values = Mtcolors4)+
  facet_grid(. ~ pop, scales = "free", space = "free") +
  theme_cowplot()+
  theme(
  axis.text.x = element_blank(),axis.title.x = element_blank(),
  #strip.background =element_blank(),strip.text.x = element_blank(),
  #legend.position = "none", 
  panel.spacing = unit(1.5, "mm"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = "dark grey"),
  plot.background = element_blank(),
  panel.border = element_rect(colour="black",fill=NA, size=0.3)
  );memprob
ggsave("mt2015_membship_apriori_40pc.png", width=10, height=5, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

#SAVE grp_membership under new name for later. 
grp_membership_MT15ap <- mutate(grp_membership, nclust="ap")

#mean membership probability for correct assginment. 
filter(grp_membership, GRP == pop)%>%group_by(pop)%>%summarize(mm=mean(MEMBSHIP))
#mean membership probability for incorrect assignment
filter(grp_membership, GRP != pop)%>%group_by(pop)%>%summarize(mm=mean(MEMBSHIP))
#total percentage of each population colored by each group. 
filter(grp_membership, GRP == pop)%>%group_by(pop)%>%summarize(tA=sum(MEMBSHIP), n_ind=n_distinct(Ind))


###################### kmeans MT 2015 ###################################
## Let the algo choose k
#loop to find the best k at different numbers of PC. 
bestk <- c()
k.bic <-c()
pc_num <- seq(5,200,5)
for (i in 1:length(pc_num)){
  grp.bic <-find.clusters.genind(mt, stat = "BIC",choose.n.clust = FALSE,max.n.clust = 20,n.pca=pc_num[i])
  k.bic[i]<-min(grp.bic$Kstat)
  bestk[i]<-match(min(grp.bic$Kstat),grp.bic$Kstat)
}
pltk <- cbind(pc_num,bestk,k.bic)%>% as.data.frame()
pltk %>% ggplot(aes(x=pc_num,y=bestk))+geom_line()+
  xlab("Number of PC retained")+ylab("Number of clusters")+
  ggtitle("Optimal number of clusters as determined by BIC")+
  theme_classic()
ggsave("mt15_KvsNPC.png", width=5, height=5, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

#get the BIC from 100 pc
grp_BIC <- find.clusters.genind(mt, stat = "BIC",choose.n.clust = FALSE,max.n.clust = 40)
as.data.frame(grp_BIC$Kstat)

#CHOOSE K
#start with 100 pc
grp_BIC <- find.clusters.genind(mt, stat = "BIC",
                                n.clust = 2, #CHOOOSEEEEEE
                                n.pca=100)

# # perform stratified cross validation
X <- scaleGen(mt, NA.method = "mean")
xval <- xvalDapc(X, grp_BIC$grp, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$`Mean Successful Assignment by Number of PCs of PCA`
xval$`Root Mean Squared Error by Number of PCs of PCA`
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`);retain 

#make scatterplot
# perform DAPC using k-mean clusters as groups
dapck <-dapc(mt, grp_BIC$grp, n.pca = retain, n.da=2) #not going to specifiy the number of da, retain 4.
scatter(dapck,grp=grp_BIC$grp, posi.da="none", bg="white", pch=17:22, scree.pca=TRUE, posi.pca="none", col=baycolors)

################LOADINGS ##### see script "loading_plots.R" for visualizations 
write.csv(dapck$var.contr,file="./DAPC/DAPC_figs/MT15_loadings_k2.csv")

tabs <-table(pop(mt), grp_BIC$grp)
#CHANGE THE NUMBER
t8 <- as.data.frame(tabs)%>%dplyr::rename("Bay"="Var1","Clust"="Var2","n_indv"="Freq")

#cluster assignments relative to original group assignments
kgrp <- strata.schemes %>% dplyr::select(Ind, pop)
kgrp$grp <-grp_BIC$grp
# show cluster distribution across populations
## Compare geographic regions of origin and assigned group memberships
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
ggsave("mt15_kmeans_membship_30pc_8clust.png", width=10, height=5, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

#dev.off()
#dim 510x345 gets rid of the white lines.

#mean membership probability. 
mjgrp <-grp_membership%>%group_by(Ind)%>%dplyr::summarize(mxmem =max(MEMBSHIP))
mjgrp<-left_join(grp_membership,mjgrp)%>%
  mutate(major.group=ifelse(mxmem==MEMBSHIP,groupnum,NA))
mjgrp <-mjgrp[!is.na(mjgrp$major.group), ] 
#mjgrp%>%group_by(pop)%>%summarize(mm=mean(MEMBSHIP))
mean(mjgrp$MEMBSHIP)

#total percentage of each population colored by each group. 
#grp_membership%>%group_by(pop)%>%summarize(tA=sum(MEMBSHIP), n_ind=n_distinct(Ind))

#SAVE each version of grp membership
#grp_membership_k2 <-mutate(grp_membership,nclust=2)%>%arrange(major.group,Ind)
#grp_membership_k3 <-mutate(grp_membership,nclust=3)%>%arrange(major.group,Ind)
#grp_membership_k4 <-mutate(grp_membership,nclust=4)%>%arrange(major.group,Ind)
#grp_membership_k5 <-mutate(grp_membership,nclust=5)%>%arrange(major.group,Ind)
#grp_membership_k6 <-mutate(grp_membership,nclust=6)%>%arrange(major.group,Ind)
#grp_membership_k7 <-mutate(grp_membership,nclust=7)%>%arrange(major.group,Ind)
grp_membership_k8 <-mutate(grp_membership,nclust=8)%>%arrange(major.group,Ind)


grp_membership_all <-bind_rows(grp_membership_k2,grp_membership_k3,grp_membership_k4,
                               grp_membership_k5,grp_membership_k6,grp_membership_k7,grp_membership_k8)%>%
  mutate(nclust=as.factor(nclust))%>%bind_rows(grp_membership_MT15ap)
tabs_all <-bind_rows(t2,t3,t4,t5,t6,t7,t8)
write.csv(tabs_all, file="./DAPC/DAPC_figs/mt2015_indbyclust.csv")
write.csv(grp_membership_all, file="./DAPC/DAPC_figs/mt2015_dapc.csv")

#try to reorder the grp membership
grp_membership_all2 <-grp_membership_all %>% 
  mutate(Ind = as.numeric(Ind))%>%
  arrange(major.group,Ind)%>%
  mutate(Ind = as.factor(Ind))

grp_membership_all2%>%
  filter(nclust !="ap")%>%
  arrange(major.group,Ind,MEMBSHIP)%>%
  ggplot(aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size=9),
    strip.background =element_rect(fill="transparent"),
    #strip.text = element_text(size=9),
    legend.position = "none",
    panel.spacing = unit(1.5, "mm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = "dark grey"),
    plot.background = element_rect(fill = NA,colour = "black", size=0.5),
    panel.border = element_rect(colour="black",fill=NA, size=0.3)
  )
ggsave("mtw2015_allMembship.png", width=8, height=10, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")


############ MATTITUCK WITHOUT 2015 ###########################################
wfpopLD4 <-wfpopLD
setPop(wfpopLD4) <-~Bay
mt <-popsub(wfpopLD4, sublist =c("Mt"))
setPop(mt) <-~Bay/Con/Year
#exclude unidentified adults but we will include the 2015 YOY for now. 
mt <-popsub(mt,exclude=c("Mt_1_2015","Mt_2_2015","Mt_5_2015", "Mt_5_2016"))
unique(mt@pop)


### A priori assignment ###
bay_dapc <- dapc(mt, n.pca = 100, n.da = 2) 
as <- optim.a.score(bay_dapc)#optimum number of pc is 1
bay_dapc <- dapc(mt, n.pca = 100, n.da = 2) 

#PC plots 
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
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`); retain 
xval$`Mean Successful Assignment by Number of PCs of PCA`
xval$`Root Mean Squared Error by Number of PCs of PCA`

Mtcolors4 <-c("#fad2e1","#d8bbff","#ff006e","#8338ec","#4361ee","#4cc9f0") 

#redo scatter
bay_dapc <- dapc(mt, n.pca = retain, n.da = 2) 
scatter(bay_dapc, posi.da="topleft",scree.pca=TRUE,posi.pca="bottomleft",legend=FALSE,solid=1.0,col=Mtcolors4) #adjust the solidity value for better plots. 

#LOADING PLOT - diverging alleles between a priori groupings. 
set.seed(4)
contrib <- loadingplot(bay_dapc$var.contr,# axis=2,
                       thres=.025, lab.jitter=1)
# Different loci contributing: WF16.250, PAM27.178, WF517.299, PAM27.292, WF32.283, WF32.213
#PLOT OF allele divergence. 
par(mfrow=c(3,2),  mai = c(0.4, 0.4, 0.4, 0.4))

#WF16
freqWF16 <- tab(genind2genpop(mt[loc=c("WF16")]),freq=TRUE)
freqWF16 <-freqWF16[,c("WF16.250")]
matplot(freqWF16,pch=c("A","B"), type="b",
        xlab="pop",ylab="allele frequency", xaxt="n",
        cex=1, main="WF16")
axis(side=1, at=1:4, lab=c("Early 2016","Late 2016","Migrants","Residents"))

#WF32
freqWF32 <- tab(genind2genpop(mt[loc=c("WF32")]),freq=TRUE)
freqWF32 <-freqWF32[,c("WF32.283", "WF32.213")]
matplot(freqWF32,pch=c("A","B"), type="b",
        xlab="pop",ylab="allele frequency", xaxt="n",
        cex=1, main="WF32")
axis(side=1, at=1:4, lab=c("Early 2016","Late 2016","Migrants","Residents"))


#PAM27
freqPAM27 <- tab(genind2genpop(mt[loc=c("PAM27")]),freq=TRUE)
freqPAM27 <-freqPAM27[,c("PAM27.178")]
matplot(freqPAM27,pch=c("A","B"), type="b",
        xlab="pop",ylab="allele frequency", xaxt="n",
        cex=1, main="PAM27")
axis(side=1, at=1:4, lab=c("Early 2016","Late 2016","Migrants","Residents"))

#WF517
freqWF517 <- tab(genind2genpop(mt[loc=c("WF517")]),freq=TRUE)
freqWF517 <-freqWF517[,c("WF517.299")]
matplot(freqWF517,pch=c("A","B"), type="b",
        xlab="pop",ylab="allele frequency", xaxt="n",
        cex=1, main="WF517")
axis(side=1, at=1:4, lab=c("Early 2016","Late 2016","Migrants","Residents"))

#WF421
freqWF421 <- tab(genind2genpop(mt[loc=c("WF421")]),freq=TRUE)
freqWF421 <-freqWF421[,c("WF421.292")]
matplot(freqWF421,pch=c("A","B"), type="b",
        xlab="pop",ylab="allele frequency", xaxt="n",
        cex=1, main="WF421")
axis(side=1, at=1:4, lab=c("Early 2016","Late 2016","Migrants","Residents"))

## Compare geographic regions of origin and assigned group memberships
grp_membership <- as.data.frame(bay_dapc$posterior) %>% tibble::rownames_to_column("Ind") %>%
  gather(key = GRP, value = MEMBSHIP, 2:5) %>% arrange(MEMBSHIP)%>%
  left_join(GRP)
#compopolot 
memprob <-ggplot(grp_membership, aes(x = Ind, y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  scale_fill_manual(values = Mtcolors4)+
  facet_grid(. ~ pop, scales = "free", space = "free") +
  theme_cowplot() ;memprob
ggsave("mt16only_membship_apriori_70pc.png", width=10, height=5, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

#SAVE grp_membership under new name for later. 
grp_membership_MT16ap <- mutate(grp_membership, nclust="ap")

#mean membership probability for correct assginment. 
filter(grp_membership, GRP == pop)%>%group_by(pop)%>%summarize(mm=mean(MEMBSHIP))
#mean membership probability for incorrect assignment
filter(grp_membership, GRP != pop)%>%group_by(pop)%>%summarize(mm=mean(MEMBSHIP))
#total percentage of each population colored by each group. 
filter(grp_membership, GRP == pop)%>%group_by(pop)%>%summarize(tA=sum(MEMBSHIP), n_ind=n_distinct(Ind))


###################### kmeans MT NO 2015 ###################################
## Let the algo choose k
#loop to find the best k at different numbers of PC. 
bestk <- c()
k.bic <-c()
pc_num <- seq(5,200,5)
for (i in 1:length(pc_num)){
  grp.bic <-find.clusters.genind(mt, stat = "BIC",choose.n.clust = FALSE,max.n.clust = 20,n.pca=pc_num[i])
  k.bic[i]<-min(grp.bic$Kstat)
  bestk[i]<-match(min(grp.bic$Kstat),grp.bic$Kstat)
}
pltk <- cbind(pc_num,bestk,k.bic)%>% as.data.frame()
pltk %>% ggplot(aes(x=pc_num,y=bestk))+geom_line()+
  xlab("Number of PC retained")+ylab("Number of clusters")+
  ggtitle("Optimal number of clusters as determined by BIC")+
  theme_classic()
ggsave("mt16_KvsNPC.png", width=5, height=5, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

grp_BIC <- find.clusters.genind(mt, stat = "BIC", choose.n.clust = F, n.pca=100)
as.data.frame(grp_BIC$Kstat)

#CHOOSE K
#start with 100 pc
grp_BIC <- find.clusters.genind(mt, stat = "BIC",
                                n.clust = 2, #CHOOOSEEEEEE
                                n.pca=100)

# # perform stratified cross validation
X <- scaleGen(mt, NA.method = "mean")
xval <- xvalDapc(X, grp_BIC$grp, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$`Mean Successful Assignment by Number of PCs of PCA`
xval$`Root Mean Squared Error by Number of PCs of PCA`
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`);retain 

#make scatterplot
# perform DAPC using k-mean clusters as groups
dapck <-dapc(mt, grp_BIC$grp, n.pca = retain, n.da=2) #not going to specifiy the number of da, retain 4.
scatter(dapck,grp=grp_BIC$grp, posi.da="bottomright", bg="white", pch=17:22, scree.pca=TRUE, posi.pca="topright", col=baycolors)

#LOADING PLOT - diverging alleles between a priori groupings. 
set.seed(4)
contrib <- loadingplot(dapck$var.contr,# axis=2,
                       thres=.04, lab.jitter=1)
# Different loci contributing: PAM27.174, PSY022.243, PSY022.235, WF12.264
#PLOT OF allele divergence. 
par(mfrow=c(2,2),  mai = c(0.4, 0.4, 0.4, 0.4))

#PSY022
freqPSY022 <- tab(genind2genpop(mt[loc=c("PSY022")]),freq=TRUE)
freqPSY022 <-freqPSY022[,c("PSY022.235","PSY022.243")]
matplot(freqPSY022,pch=c("A","B"), type="b",
        xlab="pop",ylab="allele frequency", xaxt="n",
        cex=1, main="PSY022")
axis(side=1, at=1:4, lab=c("Early 2016","Late 2016","Migrants","Residents"))

#PAM27
freqPAM27 <- tab(genind2genpop(mt[loc=c("PAM27")]),freq=TRUE)
freqPAM27 <-freqPAM27[,c("PAM27.174")]
matplot(freqPAM27,pch=c("A","B"), type="b",
        xlab="pop",ylab="allele frequency", xaxt="n",
        cex=1, main="PAM27")
axis(side=1, at=1:4, lab=c("Early 2016","Late 2016","Migrants","Residents"))

#WF12
freqWF12 <- tab(genind2genpop(mt[loc=c("WF12")]),freq=TRUE)
freqWF12 <-freqWF12[,c("WF12.264")]
matplot(freqWF12,pch=c("A","B"), type="b",
        xlab="pop",ylab="allele frequency", xaxt="n",
        cex=1, main="WF12")
axis(side=1, at=1:4, lab=c("Early 2016","Late 2016","Migrants","Residents"))

tabs <-table(pop(mt), grp_BIC$grp)
#CHANGE THE NUMBER
t8 <- as.data.frame(tabs)%>%dplyr::rename("Bay"="Var1","Clust"="Var2","n_indv"="Freq")

#cluster assignments relative to original group assignments
kgrp <- strata.schemes %>% dplyr::select(Ind, pop)
kgrp$grp <-grp_BIC$grp
# show cluster distribution across populations
## Compare geographic regions of origin and assigned group memberships
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
ggsave("mt16only_kmeans_membship_30pc_8clust.png", width=10, height=5, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

#dev.off()
#dim 510x345 gets rid of the white lines.

#mean membership probability. 
mjgrp <-grp_membership%>%group_by(Ind)%>%dplyr::summarize(mxmem =max(MEMBSHIP))
mjgrp<-left_join(grp_membership,mjgrp)%>%
  mutate(major.group=ifelse(mxmem==MEMBSHIP,groupnum,NA))
mjgrp <-mjgrp[!is.na(mjgrp$major.group), ] 
#mjgrp%>%group_by(pop)%>%summarize(mm=mean(MEMBSHIP))
mean(mjgrp$MEMBSHIP)

#total percentage of each population colored by each group. 
#grp_membership%>%group_by(pop)%>%summarize(tA=sum(MEMBSHIP), n_ind=n_distinct(Ind))

#SAVE each version of grp membership
#grp_membership_k2 <-mutate(grp_membership,nclust=2)%>%arrange(major.group,Ind)
#grp_membership_k3 <-mutate(grp_membership,nclust=3)%>%arrange(major.group,Ind)
#grp_membership_k4 <-mutate(grp_membership,nclust=4)%>%arrange(major.group,Ind)
#grp_membership_k5 <-mutate(grp_membership,nclust=5)%>%arrange(major.group,Ind)
#grp_membership_k6 <-mutate(grp_membership,nclust=6)%>%arrange(major.group,Ind)
#grp_membership_k7 <-mutate(grp_membership,nclust=7)%>%arrange(major.group,Ind)
#grp_membership_k8 <-mutate(grp_membership,nclust=8)%>%arrange(major.group,Ind)


grp_membership_all <-bind_rows(grp_membership_k2,grp_membership_k3,grp_membership_k4,
                               grp_membership_k5,grp_membership_k6,grp_membership_k7,grp_membership_k8)%>%
  mutate(nclust=as.factor(nclust))%>%bind_rows(grp_membership_MT16ap)
tabs_all <-bind_rows(t2,t3,t4,t5,t6,t7,t8)
write.csv(tabs_all, file="./DAPC/DAPC_figs/mt2016_indbyclust.csv")
write.csv(grp_membership_all, file="./DAPC/DAPC_figs/mt2016_dapc.csv")

#try to reorder the grp membership
grp_membership_all2 <-grp_membership_all %>% 
  mutate(Ind = as.numeric(Ind))%>%
  arrange(major.group,Ind)%>%
  mutate(Ind = as.factor(Ind))

grp_membership_all2%>%
  filter(nclust !="ap")%>%
  arrange(major.group,Ind,MEMBSHIP)%>%
  ggplot(aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size=9),
    strip.background =element_rect(fill="transparent"),
    #strip.text = element_text(size=9),
    legend.position = "none",
    panel.spacing = unit(1.5, "mm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = "dark grey"),
    plot.background = element_rect(fill = NA,colour = "black", size=0.5),
    panel.border = element_rect(colour="black",fill=NA, size=0.3)
  )
ggsave("mtw2016_allMembship.png", width=8, height=10, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")


####################### WITHIN SHINNECOCK ##########################################
##################################################################################

shincolors <-c("#143601","#68b0ab","#538d22","#aad576")
shincolors2 <-c("#29bf12","#abff4f","#3f8efc","#3b28cc")

wfpopLD5 <-wfpopLD
setPop(wfpopLD5) <-~Bay
sh <-popsub(wfpopLD5, sublist =c("Shin"))
setPop(sh) <-~Bay/Con/Year
unique(sh@pop)

### A priori assignment ###
#bay_dapc <- dapc(sh, n.pca = 150, n.da = 2) 
#as <- optim.a.score(bay_dapc)#optimum number of pc is 38
#bay_dapc <- dapc(sh, n.pca = 40, n.da = 2) 

#PC plots 38pc
#scatter(bay_dapc, posi.da="bottomleft",scree.pca=TRUE,posi.pca="bottomright",col=drabcolors,legend=FALSE,solid=1.0) #adjust the solidity value for better plots. 
#scatter(bay_dapc,1,1, bg="white", scree.da=FALSE, legend=TRUE, solid=.4,col=drabcolors,)

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
xval$`Mean Successful Assignment by Number of PCs of PCA`
xval$`Root Mean Squared Error by Number of PCs of PCA`
 
#redo scatter
bay_dapc <- dapc(sh, n.pca = retain, n.da = 2) 
scatter(bay_dapc, posi.da="topleft",scree.pca=TRUE,posi.pca="bottomleft",legend=FALSE,solid=1.0,col=shincolors2) #adjust the solidity value for better plots. 

scatter(bay_dapc, posi.da="none", bg="white",
        #pch=17:22, 
        pch=19,
        cstar=0, cellipse = 1.5,
        col=shincolors2, scree.pca=TRUE,
        posi.pca="none",clabel=0.9, label=c("Early 2016","Early 2017","Late 2016","Late 2017"))


################LOADINGS ##### see script "loading_plots.R" for visualizations 
write.csv(bay_dapc$var.contr,file="./DAPC/DAPC_figs/shin_loadings_apriori.csv")

## Compare geographic regions of origin and assigned group memberships
grp_membership <- as.data.frame(bay_dapc$posterior) %>% tibble::rownames_to_column("Ind") %>%
  gather(key = GRP, value = MEMBSHIP, 2:5) %>% arrange(MEMBSHIP)%>%
  left_join(GRP)
#compopolot 
memprob <-ggplot(grp_membership, aes(x = fct_reorder(Ind,MEMBSHIP), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "INDV", y = "memb. prob") +
  scale_fill_manual(values = shincolors2)+
  facet_grid(. ~ pop, scales = "free", space = "free") +
  theme_cowplot()+
  theme(
    strip.background =element_rect(fill="transparent") 
  );memprob
ggsave("shin_membship_apriori_10pc.png", width=10, height=5, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

#SAVE grp_membership under new name for later. 
grp_membership_shinap <- mutate(grp_membership, nclust="ap")

#mean membership probability for correct assginment. 
filter(grp_membership, GRP == pop)%>%group_by(pop)%>%summarize(mm=mean(MEMBSHIP))
#mean membership probability for incorrect assignment
filter(grp_membership, GRP != pop)%>%group_by(pop)%>%summarize(mm=mean(MEMBSHIP))
#total percentage of each population colored by each group. 
filter(grp_membership, GRP == pop)%>%group_by(pop)%>%summarize(tA=sum(MEMBSHIP), n_ind=n_distinct(Ind))


###################### kmeans SHINNECOCK ###################################
## Let the algo choose k
#loop to find the best k at different numbers of PC. 
bestk <- c()
k.bic <-c()
pc_num <- seq(5,200,5)
for (i in 1:length(pc_num)){
  grp.bic <-find.clusters.genind(sh, stat = "BIC",choose.n.clust = FALSE,max.n.clust = 20,n.pca=pc_num[i])
  k.bic[i]<-min(grp.bic$Kstat)
  bestk[i]<-match(min(grp.bic$Kstat),grp.bic$Kstat)
}
pltk <- cbind(pc_num,bestk,k.bic)%>% as.data.frame()
pltk %>% ggplot(aes(x=pc_num,y=bestk))+geom_line()+
  xlab("Number of PC retained")+ylab("Number of clusters")+
  ggtitle("Optimal number of clusters as determined by BIC")+
  theme_classic()
ggsave("shin_KvsNPC.png", width=5, height=5, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

#find BIC for 100pc
grp_BIC <- find.clusters.genind(sh, stat = "BIC",choose.n.clust = F, n.pca=100)
as.data.frame(grp_BIC$Kstat)

#CHOOSE K
#start with 100 pc
grp_BIC <- find.clusters.genind(sh, stat = "BIC",
                                n.clust = 2, #CHOOOSEEEEEE
                                n.pca=100)

# # perform stratified cross validation
X <- scaleGen(sh, NA.method = "mean")
xval <- xvalDapc(X, grp_BIC$grp, 
                 n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval$`Mean Successful Assignment by Number of PCs of PCA`
xval$`Root Mean Squared Error by Number of PCs of PCA`
retain <- as.numeric(xval$`Number of PCs Achieving Highest Mean Success`);retain 

#make scatterplot
# perform DAPC using k-mean clusters as groups
dapck <-dapc(sh, grp_BIC$grp, n.pca = retain, n.da=2) #not going to specifiy the number of da, retain 4.
scatter(dapck,grp=grp_BIC$grp, posi.da="bottomright", bg="white", pch=17:22, scree.pca=TRUE, posi.pca="none", col=baycolors)

################LOADINGS ##### see script "loading_plots.R" for visualizations 
write.csv(dapck$var.contr,file="./DAPC/DAPC_figs/shin_loadings_k2.csv")


tabs <-table(pop(sh), grp_BIC$grp)
#CHANGE THE NUMBER
t8 <- as.data.frame(tabs)%>%dplyr::rename("Bay"="Var1","Clust"="Var2","n_indv"="Freq")

#cluster assignments relative to original group assignments
kgrp <- strata.schemes %>% dplyr::select(Ind, pop)
kgrp$grp <-grp_BIC$grp
# show cluster distribution across populations
## Compare geographic regions of origin and assigned group memberships
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
grp_membership <-left_join(grp_membership, mjgrp)%>%
  mutate(groupnum=as.factor(groupnum))

#structure plot
#groupnum is the filling, so that you have both the filling associated with the 
memprob <-ggplot(grp_membership, aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(. ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(
    strip.background =element_rect(fill="transparent") 
  );memprob
ggsave("shin_kmeans_membship_40pc_8clust.png", width=10, height=5, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

#dev.off()
#dim 510x345 gets rid of the white lines.

#mean membership probability. 
mjgrp <-grp_membership%>%group_by(Ind)%>%dplyr::summarize(mxmem =max(MEMBSHIP))
mjgrp<-left_join(grp_membership,mjgrp)%>%
  mutate(major.group=ifelse(mxmem==MEMBSHIP,groupnum,NA))
mjgrp <-mjgrp[!is.na(mjgrp$major.group), ] 
#mjgrp%>%group_by(pop)%>%summarize(mm=mean(MEMBSHIP))
mean(mjgrp$MEMBSHIP)

#total percentage of each population colored by each group. 
#grp_membership%>%group_by(pop)%>%summarize(tA=sum(MEMBSHIP), n_ind=n_distinct(Ind))

#SAVE each version of grp membership
#grp_membership_k2 <-mutate(grp_membership,nclust=2)%>%arrange(major.group,Ind)
#grp_membership_k3 <-mutate(grp_membership,nclust=3)%>%arrange(major.group,Ind)
#grp_membership_k4 <-mutate(grp_membership,nclust=4)%>%arrange(major.group,Ind)
#grp_membership_k5 <-mutate(grp_membership,nclust=5)%>%arrange(major.group,Ind)
#grp_membership_k6 <-mutate(grp_membership,nclust=6)%>%arrange(major.group,Ind)
#grp_membership_k7 <-mutate(grp_membership,nclust=7)%>%arrange(major.group,Ind)
grp_membership_k8 <-mutate(grp_membership,nclust=8)%>%arrange(major.group,Ind)


grp_membership_all <-bind_rows(grp_membership_k2,grp_membership_k3,grp_membership_k4,
                               grp_membership_k5,grp_membership_k6,grp_membership_k7,grp_membership_k8)%>%
  mutate(nclust=as.factor(nclust))%>%bind_rows(grp_membership_shinap)
tabs_all <-bind_rows(t2,t3,t4,t5,t6,t7,t8)
write.csv(tabs_all, file="./DAPC/DAPC_figs/shin_indbyclust.csv")
write.csv(grp_membership_all, file="./DAPC/DAPC_figs/shin_dapc.csv")

############## TRY TO REORDER GRP_MEMBSHIP when making the plots. 

#try to reorder the grp membership
grp_membership_all2 <-grp_membership_all %>% 
  arrange(major.group,Ind)

grp_membership_all2%>%
  filter(nclust !="ap")%>%
  arrange(major.group,Ind,MEMBSHIP)%>%
  ggplot(aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size=9),
    strip.background =element_rect(fill="transparent"),
    #strip.text = element_text(size=9),
    legend.position = "none",
    panel.spacing = unit(1.5, "mm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = "dark grey"),
    plot.background = element_rect(fill = NA,colour = "black", size=0.5),
    panel.border = element_rect(colour="black",fill=NA, size=0.3)
  )
ggsave("shin_allMembship.png", width=8, height=10, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")


############################## DAPC FIGS ################################################################




### Relatedness #####

## NOTES ##
## includes 17 loci
## for no half 0 genotypes, use the double0 version of the files 
## population corrected for new cohort assignment
## based on the script "skreportdnaJAN_28_2020MATTITUCK.R"
#install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)


##make a box and whisker plot where we have  mean same bay pair, mean outgroup for each bay multipe boxes in multipe colors
# the color code shows the bay and the comparisons will  be demarcated.  Also Inbreeding on the same plot. 

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
library('reshape2')
library('forcats')

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")

#toolpack
sem <-function(x) sd(x)/sqrt(length(x))


##### Formating the dataset #####
wfpop <- read.genalex("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_17_sept20204genalex_doubl0.csv")
wfpop4df <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_17_sept2020_doubl0.csv", header = TRUE) #csv version 


splitStrata(wfpop) <-~Ocean/Bay/Con/Year
setPop(wfpop) <-~Bay

### Format Data as in Microsats_20.R
wfpopLD <-genclone2genind(wfpop)
all_loci <- locNames(wfpopLD)# create vector of all loci
removeloc <- c("WF27","WF06")# create vector containing loci to remove
keeploc <- setdiff(all_loci, removeloc)# create list of loci to keep
wfpopLD <- wfpopLD[loc = keeploc]# filter loci in genind object
#now re-remove the missing data. 
wfpopLD <-wfpopLD %>% missingno("loci", cutoff=0.30) # removes loci where overall missingness > 30%, so WF01
wfpopLD <- wfpopLD %>% missingno("geno", cutoff=0.25) # remove samples that are > 25% missing
length(locNames(wfpopLD))# check number of loci in genind obj


### Create the clean dataset
#remove loci previously determined to be out of HWE
setPop(wfpop) <-~Ocean
hw.test(wfpop, B=1000) #permutation based

wfpopCLEAN <-genclone2genind(wfpop)
all_loci <- locNames(wfpopCLEAN)# create vector of all loci
removeloc <- c("WF22", "WF33","WF32", "WF517", "WF223","WF01", "WF06", "WF12","PSY022","PSY087","PAM21")# create vector containing loci to remove
keeploc <- setdiff(all_loci, removeloc)# create list of loci to keep
wfpopCLEAN <- wfpopCLEAN[loc = keeploc]# filter loci in genind object
length(locNames(wfpopCLEAN))
wfpopCLEAN <- wfpopCLEAN %>% missingno("geno", cutoff=0.0) 
hw.test(wfpopCLEAN, B=1000)

#Create the YOY 2016 only dataset
popNames(wfpopLD)
setPop(wfpopLD) <-~Bay/Con/Year

wfyoy <-popsub(wfpopLD, blacklist=c("Atl_Mt_3_2015","Atl_Mt_4_2015","Atl_Mt_5_2015","Atl_Mt_3_2016","Atl_Mt_4_2016","Atl_Mt_5_2016"))
setPop(wfyoy) <-~Bay

wfyoy16 <-popsub(wfpopLD, blacklist=c("Atl_Mt_1_2015","Atl_2_2015","Atl_Mt_3_2015","Atl_Mt_4_2015","Atl_Mt_5_2015","Atl_Mt_3_2016","Atl_Mt_4_2016","Atl_Mt_5_2016", "Shin_1_2017", "Shin_2_2017"))
setPop(wfyoy16) <-~Bay

### REMEMEBER TO TURN THIS ON AND OFF. 
wfpopLD <-wfyoy16

######## Pairwise Relatedness ######

#There are many different estimators for relatedness. In previous scripts we have tried more robust ways. 
# for now, Shannon's way. https://gist.github.com/sjoleary/3efd4a7d56b115fad319781298765a31

#convert to the right data format
setPop(wfpopLD) <-~Bay
df <-genind2df(wfpopLD, usepop = FALSE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(35,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34)]
df[df=="NA"] <- 0 # missing data must be 0
write.table(df, "scratch",row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) #write it.
genotypedata <- readgenotypedata("scratch")# import input file as list (gdata, nloci, nalleles, ninds, freqs)

#if you don't want to do them all again, try just doing one... 
simdata <-familysim(genotypedata$freqs,100)
outputTRI <- coancestry(simdata, trioml=1)
simreltri <- cleanuprvals(outputTRI$relatedness , 100)

#custom comparisons;
simdata <-familysim(genotypedata$freqs,100)
output <- coancestry(simdata, lynchrd=1, trioml=1, quellergt=1, wang=1)
simrel <- cleanuprvals(output$relatedness , 100)
riomlpo <- simrel[1:100,5]
triomlfs <-simrel[(100+1):(2*100),5]
triomlhs <-simrel[((2*100) + 1): (3*100),5]
triomlur <-simrel[((3*100)+1):(4*100),5]
wangpo <- simrel[1:100,6]
wangfs <-simrel[(100+1):(2*100),6]
wanghs <-simrel[((2*100) + 1): (3*100),6]
wangur <-simrel[((3*100)+1):(4*100),6]
quellerpo <- simrel[1:100,10]
quellerfs <-simrel[(100+1):(2*100),10]
quellerhs <-simrel[((2*100) + 1): (3*100),10]
quellerur <-simrel[((3*100)+1):(4*100),10]
lynchrdpo <- simrel[1:100,8]
lynchrdfs <-simrel[(100+1):(2*100),8]
lynchrdhs <-simrel[((2*100) + 1): (3*100),8]
lynchrdur <-simrel[((3*100)+1):(4*100),8]

trioml <-rep("tri",100)
wang <-rep("wang",100)
lynchrd <-rep("lynchrd",100)
quellergt <-rep("queller",100)
estimator2 <- c(trioml, wang, quellergt, lynchrd) 
Estimator <- rep(estimator2 , 4)

po <- rep("Parent-Offspring", (4 * 100)) 
fs <- rep("Full-Sibs", (4 * 100))
hs <- rep("Half-Sibs", (4 * 100)) 
ur <- rep("Unrelated", (4 * 100))
relationship <- c(po, fs, hs, ur)

relatednesspo <- c(riomlpo , wangpo , quellerpo , lynchrdpo)
relatednessfs <- c(triomlfs , wangfs , quellerfs , lynchrdfs)
relatednesshs <- c(triomlhs , wanghs , quellerhs , lynchrdhs)
relatednessur <- c(triomlur , wangur , quellerur , lynchrdur)
Relatedness_Value <- c(relatednesspo , relatednessfs , relatednesshs , relatednessur)

combineddata <- as.data.frame(cbind(Estimator , relationship , Relatedness_Value))
combineddata$Relatedness_Value <- as.numeric(as.character( combineddata$Relatedness_Value))

ggplot(combineddata , aes(x = Estimator , y = Relatedness_Value), ylim = c(-0.5, 1.0)) +
  geom_boxplot(fill="light grey") +
  facet_wrap(~ relationship)+
  scale_x_discrete(labels = c("LR","QG","TR","WG"))+
  xlab('Estimator')+ylab("Relatedness")+theme_classic()+guides(fill = FALSE, colour = FALSE)

#calculate correlation coefficient between observed values for each estimator and the expected values. 
urval <- rep(0, 100) 
hsval <- rep (0.25 , 100)
fsval <- rep (0.5 , 100)
poval <- rep (0.5 , 100)
relvals <- c(poval , fsval , hsval , urval)

cor(relvals , simrel[, 5]) #tri
cor(relvals , simrel[, 6]) #wang
cor(relvals , simrel[, 10]) #queller
cor(relvals , simrel[, 8]) #lynch


#extract mean, median & CI values for relatedness from the simulation. 
relsim <- as.data.frame(combineddata)
relply<- ddply(relsim, Estimator~relationship,summarize, mean.rel = mean(Relatedness_Value), LIrel = quantile(Relatedness_Value,0.05), HIrel = quantile(Relatedness_Value, 0.95), medianrel= quantile(Relatedness_Value, 0.5))
relply

#from the abbreviated simulation
relsim <- as.data.frame(simreltri)
relply<- ddply(relsim, ~group,summarize, mean.rel = mean(trioml), LIrel = quantile(trioml,0.05), HIrel = quantile(trioml, 0.95), medianrel= quantile(trioml, 0.5))
relply


#trioml estimate (Wang) - Can't do this because it crashes R studio. 
relatedness_triad <- coancestry(genotypedata$gdata,trioml = 1)
relatedT <- relatedness_triad$relatedness %>%
  dplyr::select(pair.no, ind1.id, ind2.id, trioml)

library("readr")
write_delim(relatedT, "pairwise_relatedness_triomlYOY")

# write inbreeding to file
inbreed <- relatedness_triad$inbreeding %>%
  dplyr::select(ind.id, L3, LH) %>%
  dplyr::rename(INDV = ind.id)
write_delim(inbreed, "inbreedingYOY")

######Inbreeding########
# inbreeding  - you can put the values from the simulation in later. 
ggplot(inbreed, aes(x = L3)) +
  geom_histogram(binwidth = 0.005, color = "light grey", fill = "dark grey") +
  geom_vline(aes(xintercept = mean(LH, na.rm = TRUE)), #dark blue is the mean
             color = "black", linetype = "dashed", size = 0.5) +
  theme_cowplot()+
  labs(x = "inbreeding coefficient (Fis)", y = "individuals") 
ggsave('inbreedldYOY.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")

#inbreeding by bay. 
popinfo <-dplyr::select(wfpop4df, Ind, Pop) %>% separate(Pop, c("Ocean","Bay","Con","Year"), remove=FALSE)
wf.df <-mutate(popinfo,ind1.id=Ind,ind2.id=Ind)
wf.df <-dplyr::rename(wf.df, INDV=Ind)
inbreed <-left_join(inbreed, wf.df, by="INDV")

inbreed %>%
  ggplot(aes(x=fct_rev(Bay),y=L3))+ 
  #geom_boxplot(aes(fill=Bay))+ coord_flip()+
  ggplot2::stat_summary(fun.data = mean_sdl,fun.args = list(mult = 1),geom = "pointrange",position = ggplot2::position_nudge(x = 0.05, y = 0)) +coord_flip()+ 
  geom_flat_violin(aes(fill=Bay),position = position_nudge(x = .1, y = 0),adjust=2, trim = FALSE)+
  geom_hline(aes(yintercept = mean(L3, na.rm = TRUE)),color = "black", linetype = "dashed", size = 0.5) +
  scale_fill_manual(name = "Bay",values = drabcolors)+
  xlab(' ')+ylab("Inbreeding Coefficient (Fis)")+ 
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('indvinbreeding.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 8, height = 10)
#it may be more appropriate to compare YOY only. 

## anova to compare inbreeding between bays. 
library("agricolae")
anoin <-lm(L3~Bay, na.action=na.omit, data=inbreed)
ano2 <-car::Anova(anoin)
ano2
ano3<-df.residual(anoin)
MSerror<-deviance(anoin)/ano3
comparison <- HSD.test(anoin,c("Bay"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05, group=TRUE)
comparison



# relatedness (TRIOML)- remember to change these values depending on what dataset you are using. 
relatedT %>%
  filter(trioml > 0) %>%
ggplot(aes(x = trioml)) +
  geom_histogram(binwidth = 0.001, color = "light grey", fill = "light grey") +
  geom_vline(aes(xintercept = mean(trioml, na.rm = TRUE)),
             color = "black", linetype = "dashed", size = 0.5) +
  #geom_vline(aes(xintercept = quantile(lynchrd, 0.95)),
             #color = "darkred", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = 0.04), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for unrelated in the simulation
  geom_vline(aes(xintercept = 0.25), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for half sibs in the simulation
  geom_vline(aes(xintercept = 0.52), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for full sibs in the simulation
  #geom_vline(aes(xintercept = 0.54), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for parent offspring in the simulation
  #geom_vline(aes(xintercept = 0.08),color = "darkblue", linetype = "dotted", size = 0.5) + #the Lynchrd 0.05 quantile on the simulated half sibs relationship
  #geom_vline(aes(xintercept = 0.49), color = "darkblue", linetype = "dotted", size = 0.5) + #red is the 0.95 quantile
  labs(x = "relatedness", y = "number of pairs")+
  theme_cowplot()
ggsave("pairwiserelatenessTRIOMLYOY_17.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
dev.off()


#attach strata information to the relatedness estimation. 
h <-melt(relatedT, id=c("pair.no","ind1.id","ind2.id"))
h <-dplyr::rename(h, Estimator=variable,Relatedness_Value=value)

#how many of each? (based on estimator mean)
filter(relatedT, trioml > 0.099575) %>% n_distinct()
filter(relatedT, trioml > 0.25) %>% n_distinct()
filter(relatedT,trioml > 0.25 & trioml <= 0.52) %>% n_distinct() #how many half sibs
filter(relatedT, trioml > 0.52) %>% n_distinct() #full sibs or parent offspring. 
filter(relatedT, trioml > 0.54) %>% n_distinct() #parent offspring
n_distinct(relatedT)


#create a dataset that combines fish pairs with information about their population
popinfo <-dplyr::select(wfpop4df, Ind, Pop) %>% separate(Pop, c("Ocean","Bay","Con","Year"), remove=FALSE)
wf.df <-mutate(popinfo,ind1.id=Ind,ind2.id=Ind)
hs <-left_join(h,wf.df, by=c("ind1.id"))
hs <-dplyr::rename(hs, Bay1=Bay,Con1=Con,Year1=Year)
hs <-dplyr::select(hs,-ind2.id.y) %>% dplyr::rename(ind2.id =ind2.id.x)
hs <-left_join(hs,wf.df,by=c("ind2.id")) %>% dplyr::rename(Bay2=Bay,Con2=Con,Year2=Year)
hs <-dplyr::rename(hs,ind1.id=ind1.id.x) %>% dplyr::select(-ind1.id.y,-Ind.y, -Ind.x)
#get rid of "sketchy" individuals
hs <-filter(hs,!is.na(Bay2))
hs <-filter(hs,!is.na(Bay1))

#find related fish from the same bay!
same.bay <-filter(hs, Bay1==Bay2)
#mean relatedness value of same bay pairs
ddply(same.bay, ~Bay1, summarize, avr = mean(Relatedness_Value), sdrel =sd(Relatedness_Value),se= sem(Relatedness_Value),LCI = quantile(Relatedness_Value, 0.025), UCI= quantile(Relatedness_Value, 0.975))
filter(same.bay,Relatedness_Value >= 0.25 & Relatedness_Value <= 0.52) %>% n_distinct() #how many half sibs
filter(same.bay,Relatedness_Value > 0.52) %>% n_distinct() #how many full sibs

#relatedness diff bays
diff.bay <-filter(hs, Bay1!=Bay2)
#mean relatedness value of diff bay pairs
ddply(diff.bay, Bay2~Bay1, summarize, avr = mean(Relatedness_Value), sdrel =sd(Relatedness_Value),se= sem(Relatedness_Value),LCI = quantile(Relatedness_Value, 0.025), UCI= quantile(Relatedness_Value, 0.975))
filter(diff.bay,Relatedness_Value >= 0.25 & Relatedness_Value <= 0.52) %>% n_distinct() #how many half sibs
filter(diff.bay,Relatedness_Value > 0.52) %>% n_distinct() #how many full sibs

##barplot of mean relatedness from same bay pairs. 
drabcolors <-c("#d0d1e6","#a6bddb", "#67a9cf", "#1c9099", "#016450")

same.bay %>%
  ggplot(aes(x=fct_rev(Bay1),y=Relatedness_Value),fill=Bay1)+
  #geom_boxplot(aes(fill=Bay1))+ 
  ggplot2::stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "pointrange",
    position = ggplot2::position_nudge(x = 0.05, y = 0)
  ) +
  geom_flat_violin(aes(fill=Bay1),position = position_nudge(x = .1, y = 0),adjust=2, trim = FALSE)+
  scale_fill_manual(name = "Bay",values = drabcolors)+coord_flip()+ 
  geom_hline(aes(yintercept = 0.25), color = "black", linetype = "dotted", size = 0.7) + #the estimated mean value for half sibs in the simulation
  geom_hline(aes(yintercept = 0.52), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for full sibs in the simulation
  #geom_hline(aes(yintercept = 0.50), color = "darkgreen", linetype = "dashed", size = 0.7) + #the estimated mean value for parent offspring in the simulation
  xlab(' ')+ylab("Mean Relatedness")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('samebay_YOY.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

##violin plot of diff bay pairs
diff.bay <-mutate(diff.bay, pair1 = pmin(Bay1,Bay2), pair2 =pmax(Bay1,Bay2)) %>% 
  unite(BayPair, pair1,pair2, sep="-", remove=TRUE) %>% unique() 
  
diff.bay %>%
  #ggplot(aes(x=fct_rev(BayPair),y=Relatedness_Value))+ # bay pairs
  ggplot(aes(x=fct_rev(Bay1),y=Relatedness_Value))+ # bay pairs
  #geom_boxplot(aes(fill=Bay1))+ 
  ggplot2::stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "pointrange",
    position = ggplot2::position_nudge(x = 0.05, y = 0)
  ) +coord_flip()+ 
  geom_flat_violin(aes(fill=Bay1),position = position_nudge(x = .1, y = 0),adjust=2, trim = FALSE)+
  scale_fill_manual(name = "Bay",values = drabcolors)+
  geom_hline(aes(yintercept = 0.25), color = "black", linetype = "dotted", size = 0.7) + #the estimated mean value for half sibs in the simulation
  geom_hline(aes(yintercept = 0.52), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for full sibs in the simulation
  #geom_hline(aes(yintercept = 0.50), color = "darkgreen", linetype = "dashed", size = 0.7) + #the estimated mean value for parent offspring in the simulation
  xlab(' ')+ylab("Mean Relatedness")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('diffbay_YOY.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

#in groups and out groups
same.bay <-mutate(same.bay, io = "same")
diff.bay <-mutate(diff.bay, io = "diff")
inout <- bind_rows(same.bay, diff.bay)
inout <-unite(inout, ioo, Bay1, io, sep="_", remove=FALSE)

inout %>%
  ggplot(aes(x=fct_rev(Bay1),y=Relatedness_Value))+ # bay pairs
  #geom_boxplot(aes(fill=Bay1))+ coord_flip()+
  ggplot2::stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "pointrange",
    position = ggplot2::position_nudge(x = 0.05, y = 0)
  ) +coord_flip()+ 
  geom_flat_violin(aes(fill=Bay1),position = position_nudge(x = .1, y = 0),adjust=2, trim = FALSE)+
  scale_fill_manual(name = "Bay",values = drabcolors)+
  geom_hline(aes(yintercept = 0.04), color = "grey", linetype = "dotted", size = 0.7) + #the estimated mean value for unrelated indv in the simulation
  geom_hline(aes(yintercept = 0.25), color = "grey", linetype = "dotted", size = 0.7) + #the estimated mean value for half sibs in the simulation
  geom_hline(aes(yintercept = 0.52), color = "grey", linetype = "dotted", size = 0.7) + #the estimated mean value for full sibs in the simulation
  #geom_hline(aes(yintercept = 0.54), color = "grey", linetype = "dotted", size = 0.7) + #the estimated mean value for full sibs in the simulation
   xlab(' ')+ylab("Mean Relatedness")+ 
  facet_grid(~io, scales="free")+
  theme(axis.text = element_text(size = 16),axis.title = element_text(size = 16),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('inandoutgroupsYOY.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 6)
#it may be more appropriate to compare YOY only. 

## anova to compare relatedness and ingroup outgroup
library("agricolae")
anorel <-lm(Relatedness_Value~io*Bay1, na.action=na.omit, data=inout)
ano2 <-car::Anova(anorel)
ano2
ano3<-df.residual(anorel)
MSerror<-deviance(anorel)/ano3
comparison <- HSD.test(anorel,c("io","Bay1"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05, group=TRUE)
comparison

## anova to compare relatedness and ingroup outgroup
anorel <-lm(Relatedness_Value~Bay1, na.action=na.omit, data=same.bay)
ano2 <-car::Anova(anorel)
ano2
ano3<-df.residual(anorel)
MSerror<-deviance(anorel)/ano3
comparison <- HSD.test(anorel,c("Bay1"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05, group=TRUE)
comparison




############ group relatedness ##############
#format data so that first two letters of IND denote group name, then into gdata format
wf.df <-dplyr::select(wf.df, -ind1.id, -ind2.id) 
wf.df <-dplyr::rename(wf.df, Ind=INDV)
group.df <-left_join(df, wf.df, by="Ind") %>% mutate(pop.prefix = substr(Bay, 1,2)) %>% 
  unite(name, pop.prefix, Ind, sep="_", remove =TRUE) %>% dplyr::select(-Pop, -Ocean, -Bay, -Con, -Year)
write.table(group.df, "scratch",row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) #write it.
genotypedata <- readgenotypedata("scratch")# import input file as list (gdata, nloci, nalleles, ninds, freqs)


#re-run overnight
relbay <- grouprel(genotypes = genotypedata$gdata, estimatorname = "trioml", usedgroups = "all", iterations= 100)

#re run 
shinrel <-grouprel(genotypes = genotypedata$gdata, estimatorname = "trioml", usedgroups = "Sh", iterations= 100)

#for now...
exprel <-read.csv("expectedrel.csv",header=TRUE)
names(exprel) <-c("sim","Napeague","Moriches","Jamaica","Mattituck","Shinnecock", "ALL")

filter(exprel, Napeague >= 0.03695569) %>% n_distinct #61/100
filter(exprel, Moriches >= 0.03703134) %>% n_distinct #64/100
filter(exprel, Jamaica >= 0.03799505) %>% n_distinct #59
filter(exprel, Shinnecock >= 0.04365348) %>% n_distinct #1  <_this seems wrong. 
filter(exprel, Mattituck >= 0.04013806) %>% n_distinct #42

exprel %>%
  ggplot(aes(x=ALL))+
  geom_histogram(bins=20)
 

## make a version of this analysis done with the clean dataset and run overnight. 

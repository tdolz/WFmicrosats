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
#geom flat_violin:
  #https://gist.github.com/dgrtwo/eb7750e74997891d7c20
drabcolors <-c("#d0d1e6","#a6bddb", "#67a9cf", "#1c9099", "#016450")
shincolors2 <-c("#29bf12","#abff4f","#3f8efc","#3b28cc")
shincolors <-c("#143601","#68b0ab","#538d22","#aad576")
Mtcolors <-c("#f4d35e","#ee964b","#f95738","#ee4266","#15616d","#0d3b66")
Mtcolors2 <-c("#f72585","#7209b7","#4361ee","#4cc9f0")


##### Formating the dataset #####
wfpop <- read.genalex("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_17_sept20204genalex_doubl0ABC.csv")
wfpop4df <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_17_sept2020_doubl0ABC.csv", header = TRUE) #csv version 


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
#hw.test(wfpop, B=1000) #permutation based

wfpopCLEAN <-genclone2genind(wfpop)
all_loci <- locNames(wfpopCLEAN)# create vector of all loci
removeloc <- c("WF22", "WF33","WF32", "WF517", "WF223","WF01", "WF06", "WF12","PSY022","PSY087","PAM21")# create vector containing loci to remove
keeploc <- setdiff(all_loci, removeloc)# create list of loci to keep
wfpopCLEAN <- wfpopCLEAN[loc = keeploc]# filter loci in genind object
length(locNames(wfpopCLEAN))
wfpopCLEAN <- wfpopCLEAN %>% missingno("geno", cutoff=0.0) 
#hw.test(wfpopCLEAN, B=1000)

#Create the YOY 2016 only dataset
popNames(wfpopLD)
setPop(wfpopLD) <-~Bay
setPop(wfpopLD) <-~Ocean/Bay/Con/Year

wfyoy <-popsub(wfpopLD, blacklist=c("Atl_Mt_3_adults","Atl_Mt_4_adults","Atl_Mt_5_2015","Atl_Mt_5_2016"))
setPop(wfyoy) <-~Bay

wfyoy16 <-popsub(wfpopLD, blacklist=c("Atl_Mt_1_2015","Atl_Mt_2_2015","Atl_Mt_3_adults","Atl_Mt_4_adults","Atl_Mt_5_2015","Atl_Mt_5_2016", "Atl_Shin_1_2017", "Atl_Shin_2_2017"))
setPop(wfyoy16) <-~Bay

shinyoy <- popsub(wfpopLD, sublist=c("Atl_Shin_1_2016", "Atl_Shin_2_2016", "Atl_Shin_1_2017","Atl_Shin_2_2017"))
setPop(shinyoy) <-~Bay/Con/Year

setPop(wfpopLD) <-~Bay/Con
mtpop <-popsub(wfpopLD, sublist=c("Mt_1","Mt_2","Mt_3","Mt_4"))
setPop(mtpop) <-~Bay/Con/Year

### REMEMEBER TO TURN THIS ON AND OFF. 
wfpopLD <-wfpopLD
setPop(wfpopLD) <-~Bay/Year
######## Pairwise Relatedness ######

#There are many different estimators for relatedness. In previous scripts we have tried more robust ways. 
# for now, Shannon's way. https://gist.github.com/sjoleary/3efd4a7d56b115fad319781298765a31

#convert to the right data format
setPop(wfpopLD) <-~Bay/Year
df <-genind2df(wfpopLD, usepop = FALSE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(35,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34)]
df[df=="NA"] <- 0 # missing data must be 0
write.table(df, "scratch",row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) #write it.
genotypedata <- readgenotypedata("scratch")# import input file as list (gdata, nloci, nalleles, ninds, freqs)

### for the clean dataset analysis only. 
#setPop(wfpopLD) <-~Bay
#df <-genind2df(wfpopLD, usepop = FALSE, oneColPerAll = TRUE) 
#df$Ind <- rownames(df)
#df <-df[,c(19,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)]
#df[df=="NA"] <- 0 # missing data must be 0
#write.table(df, "scratch",row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) #write it.
#genotypedata <- readgenotypedata("scratch")# import input file as list (gdata, nloci, nalleles, ninds, freqs)



#if you don't want to do them all again, try just doing one... 
simdata <-familysim(genotypedata$freqs,100)
outputTRI <- coancestry(simdata, trioml=2)
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
  relsim2 <- as.data.frame(simreltri)
  relply2<- ddply(relsim2, ~group,summarize, mean.rel = mean(trioml),semrel =sem(trioml), LIrel = quantile(trioml,0.05), HIrel = quantile(trioml, 0.95), medianrel= quantile(trioml, 0.5))
  relply2

#thresholds - you want to know the mean but also the thresholds.
halfsibs <- relply2[1,2]
halfsibs.LCI <-relply2[1,4]
halfsibs.HCI <-relply2[1,5]
fullsibs <-relply2[3,2]
fullsibs.LCI <-relply2[3,4]
fullsibs.HCI <-relply2[3,5]
parentoffspring <-relply2[2,2] #not interesting when only doing yoy
unrelated <-relply2[4,2]
halfsibs.sem <-relply2[1,3]
fullsibs.sem <-relply2[3,3]
po.sem <-relply2[2,3]
unre.sem <-relply2[4,3]

#trioml estimate (Wang) 
relatedness_triad <- coancestry(genotypedata$gdata,trioml = 2) #remember to turn on 1 or 2, for no CI or w/ CI. 
relatedT <- relatedness_triad$relatedness %>%
  dplyr::select(pair.no, ind1.id, ind2.id, trioml)

relatedCI <- relatedness_triad$relatedness.ci95 %>%
  dplyr::select(pair.no, ind1.id, ind2.id, trioml.low, trioml.high)

relatedT <-left_join(relatedT, relatedCI)

library("readr")
write_delim(relatedT, "pairwise_relatedness_trioml")


# relatedness (TRIOML)- remember to change these values depending on what dataset you are using. 
relatedT %>%
  #filter(trioml > 0) %>%
ggplot(aes(x = trioml)) +
  geom_histogram(binwidth = 0.001, color = "light grey", fill = "light grey") +
  geom_vline(aes(xintercept = mean(trioml, na.rm = TRUE)),
             color = "black", linetype = "dashed", size = 0.5) +
  ylim(0, 500)+
  geom_vline(aes(xintercept = unrelated), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for unrelated in the simulation
  geom_vline(aes(xintercept = halfsibs), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for half sibs in the simulation
  geom_vline(aes(xintercept = fullsibs), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for full sibs in the simulation
  labs(x = "relatedness", y = "number of pairs")+
  theme_cowplot()
ggsave("pairwiserelatenessTRIall.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
dev.off()


#attach strata information to the relatedness estimation. 
hT <-pivot_longer(relatedT, cols=c("trioml","trioml.low","trioml.high"))
h <-dplyr::rename(relatedT, Relatedness_Value=trioml)

#how many of each? (based on estimator mean)
n_distinct(relatedT)
#percentage SIGNFICANTLY different from ZERO
(filter(relatedT, trioml.low > 0) %>% n_distinct())/n_distinct(relatedT$pair.no)
#relatedness values significantly LESS than the mean halfsibs value
(filter(relatedT, trioml.high < (halfsibs-halfsibs.sem)) %>% n_distinct())/n_distinct(relatedT$pair.no)

#create a dataset that combines fish pairs with information about their population
popinfo <-dplyr::select(wfpop4df, Ind, Pop) %>% 
  separate(Pop, c("Ocean","Bay","Con","Year"), remove=FALSE)%>%
  unite(ConYear, Con, Year, remove=FALSE)
wf.df <-mutate(popinfo,ind1.id=Ind,ind2.id=Ind)
hs <-left_join(h,wf.df, by=c("ind1.id"))
hs <-dplyr::rename(hs, Bay1=Bay,Con1=Con,Year1=Year, ConYear1 = ConYear)
hs <-dplyr::select(hs,-ind2.id.y) %>% dplyr::rename(ind2.id =ind2.id.x)
hs <-left_join(hs,wf.df,by=c("ind2.id")) %>% dplyr::rename(Bay2=Bay,Con2=Con,Year2=Year,ConYear2 = ConYear)
hs <-dplyr::rename(hs,ind1.id=ind1.id.x) %>% dplyr::select(-ind1.id.y,-Ind.y, -Ind.x)
#get rid of "sketchy" individuals
hs <-filter(hs,!is.na(Bay2))
hs <-filter(hs,!is.na(Bay1))
hs <-mutate(hs, pair1 = pmin(Bay1,Bay2), pair2 =pmax(Bay1,Bay2)) %>% 
  unite(BayPair, pair1,pair2, sep="-", remove=TRUE) %>% unique() 
hs <-unite(hs, BayYear1, Bay1, Year1, sep="_", remove=FALSE) %>% unite(BayYear2, Bay2, Year2, sep="_", remove=FALSE)
hs <-mutate_at(hs, vars(Pop.x, Bay1,ConYear1,Year1,BayYear1, BayYear2, Pop.y, Bay2, ConYear2, Year2), funs(as.factor(.)))

same.bay <-filter(hs, Bay1==Bay2)
diff.bay <-filter(hs, Bay1!=Bay2)

#different bays full siblings
(filter(diff.bay, trioml.low > (fullsibs - fullsibs.sem)) %>% n_distinct())/n_distinct(relatedT$pair.no)
#different bays half siblings: significant
(filter(diff.bay, trioml.low > (halfsibs - halfsibs.sem)) %>% n_distinct())/n_distinct(relatedT$pair.no)
#different bays half siblings: all, but DIFFERENT BAYS
(filter(diff.bay, Relatedness_Value > halfsibs) %>% n_distinct())/n_distinct(relatedT$pair.no)

#Within group pairwise relatedness: Different between bays and years? 
same.bay <-mutate(same.bay, GRP = Bay1) %>% mutate(level="bay")
same.year <-filter(same.bay, Year1==Year2)%>% mutate(GRP=BayYear1, level="year")
same.CY <-filter(same.bay, Pop.x==Pop.y) %>% mutate(GRP=Pop.x,level=ifelse(Bay1=="Mt","Mt","Sh"))
same.grp <-bind_rows(same.bay, same.year, same.CY)
same.grp <-mutate(same.grp, GRP=as.factor(GRP))

allcolors <-c("#d0d1e6","#a6bddb","#67a9cf","#1c9099","#016450","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#016450","#016450","#016450","#016450","#016450")
same.grp <-filter(same.grp,GRP %in% c("Nap","Mor","Jam","Shin","Shin_2016","Shin_2017","Atl_Shin_1_2016","Atl_Shin_2_2016","Atl_Shin_1_2017","Atl_Shin_2_2017","Mt","Mt_2015","Mt_adults","Mt_2016", "Atl_Mt_3_adults","Atl_Mt_4_adults","Atl_Mt_1_2015","Atl_Mt_2_2015","Atl_Mt_1_2016","Atl_Mt_2_2016"))
same.grp$GRP <-fct_relevel(same.grp$GRP, c("Jam","Mor","Mt","Mt_2015","Mt_2016","Mt_adults","Nap","Shin","Shin_2016","Shin_2017","Atl_Mt_1_2015","Atl_Mt_2_2015","Atl_Mt_1_2016","Atl_Mt_2_2016", "Atl_Mt_3_adults","Atl_Mt_4_adults","Atl_Shin_1_2016","Atl_Shin_2_2016","Atl_Shin_1_2017","Atl_Shin_2_2017"))
yearcols <- c("#d0d1e6","#a6bddb","#67a9cf","#67a9cf","#67a9cf","#1c9099","#016450","#016450")

##barplot of mean relatedness from same bay pairs. 
dodge <- position_dodge(width = 0.9)
same.year %>%
  ggplot(aes(x=fct_rev(GRP),y=Relatedness_Value,fill=GRP))+
  geom_point(aes(colour = GRP), position = position_jitterdodge(dodge.width = 0.9), alpha = 0.8) +
  #geom_boxplot(outlier.colour = NA, position = dodge)+ 
  ggplot2::stat_summary(fun.data = mean_sdl,fun.args = list(mult = 1),geom = "pointrange",position = ggplot2::position_nudge(x = 0.05, y = 0)) +
  geom_flat_violin(aes(fill=GRP),position = position_nudge(x = .1, y = 0),adjust=2, trim = FALSE)+
  scale_fill_manual(name = "GRP",values = yearcols)+
  scale_color_manual(name = "GRP",values = yearcols)+
  geom_hline(aes(yintercept = halfsibs), color = "black", linetype = "dotted", size = 0.7) + #the estimated mean value for half sibs in the simulation
  geom_hline(aes(yintercept = fullsibs), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for full sibs in the simulation
  xlab(' ')+ylab("Pairwise Relatedness")+ coord_flip()+
  #facet_wrap(~level, scales="free")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('sameyearboxJitter.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)

##violin plot of diff bay pairs
diff.bay %>%
  ggplot(aes(x=fct_rev(BayPair),y=Relatedness_Value))+ # bay pairs
  geom_point(aes(colour = "BayPair"), position = position_jitterdodge(dodge.width = 0.9), alpha = 0.8) +
  #geom_boxplot(fill="light grey")+ 
  ggplot2::stat_summary(fun.data = mean_sdl,fun.args = list(mult = 1),geom = "pointrange",position = ggplot2::position_nudge(x = 0.05, y = 0))+
  geom_flat_violin(aes(fill="light grey"),fill="light grey",position = position_nudge(x = .1, y = 0),adjust=2, trim = FALSE)+
  coord_flip()+ 
  scale_color_manual(name = "BayPair", values = c(rep("light grey",10)))+
  geom_hline(aes(yintercept = halfsibs), color = "black", linetype = "dotted", size = 0.7) + #the estimated mean value for half sibs in the simulation
  geom_hline(aes(yintercept = fullsibs), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for full sibs in the simulation
  xlab(' ')+ylab("Pairwise Relatedness")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('diffbayjitter.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)

### diff conyears
diff.CY <-filter(same.bay, Pop.x!=Pop.y) %>% unite(CYpair, ConYear1, ConYear2, sep="-", remove=FALSE)
diff.CY %>%
  ggplot(aes(x=fct_rev(CYpair),y=Relatedness_Value))+ # CY pairs
  geom_boxplot(fill="light grey")+ coord_flip()+ 
  geom_hline(aes(yintercept = halfsibs), color = "black", linetype = "dotted", size = 0.7) + #the estimated mean value for half sibs in the simulation
  geom_hline(aes(yintercept = fullsibs), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for full sibs in the simulation
  xlab(' ')+ylab("Pairwise Relatedness")+
  facet_wrap(~Bay1, scales="free")+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"))+guides(fill = FALSE, colour = FALSE) 
ggsave('diffCY.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 12)

### diff bayears
diff.BY <-filter(hs, BayYear1!=BayYear2) %>% unite(BYpair, BayYear1, BayYear2, sep="-", remove=FALSE)
diff.BY %>%
  ggplot(aes(x=fct_rev(BYpair),y=Relatedness_Value))+ # CY pairs
  geom_boxplot(fill="light grey")+ coord_flip()+ 
  geom_hline(aes(yintercept = halfsibs), color = "black", linetype = "dotted", size = 0.7) + #the estimated mean value for half sibs in the simulation
  geom_hline(aes(yintercept = fullsibs), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for full sibs in the simulation
  xlab(' ')+ylab("Pairwise Relatedness")+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"))+guides(fill = FALSE, colour = FALSE) 
ggsave('diffBAYYEAR.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 12)



######in groups and out groups#####
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
  #scale_fill_manual(name = "Bay",values = shincolors2)+
  scale_fill_grey()+
  geom_hline(aes(yintercept = unrelated), color = "grey", linetype = "dotted", size = 0.7) + #the estimated mean value for unrelated indv in the simulation
  geom_hline(aes(yintercept = halfsibs), color = "grey", linetype = "dotted", size = 0.7) + #the estimated mean value for half sibs in the simulation
  geom_hline(aes(yintercept = fullsibs), color = "grey", linetype = "dotted", size = 0.7) + #the estimated mean value for full sibs in the simulation
  xlab(' ')+ylab("Pairwise Relatedness")+ 
  facet_grid(~io, scales="free")+
  theme(axis.text = element_text(size = 16),axis.title = element_text(size = 16),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('inandoutgroupsMT.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 6)
#### not sure we will include ###

# same bays
anorel <-lm(Relatedness_Value~Bay1*Year1, na.action=na.omit, data=same.bay)
ano2 <-car::Anova(anorel)
ano2
ano3<-df.residual(anorel)
MSerror<-deviance(anorel)/ano3
comparison <- HSD.test(anorel,c("Bay1","Year1"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/8, group=TRUE)
comparison

#same bays just bay
anorel <-lm(Relatedness_Value~Bay1, na.action=na.omit, data=same.bay)
ano2 <-car::Anova(anorel)
ano2
ano3<-df.residual(anorel)
MSerror<-deviance(anorel)/ano3
comparison <- HSD.test(anorel,c("Bay1"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/5, group=TRUE)
comparison

ddply(same.bay, Bay1~Year1, summarize, averel=mean(Relatedness_Value),serel=sem(Relatedness_Value))

#diff bays
anorel <-lm(Relatedness_Value~BayPair, na.action=na.omit, data=diff.bay)
ano2 <-car::Anova(anorel)
ano2
ano3<-df.residual(anorel)
MSerror<-deviance(anorel)/ano3
comparison <- HSD.test(anorel,c("BayPair"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/10, group=TRUE)
comparison

#hs2 <-filter(hs, Bay1=="Shin" & Bay2=="Shin")
hs2 <-filter(hs, Bay1=="Mt" & Bay2=="Mt")

half_sibs <-filter(hs2, Relatedness_Value > halfsibs & Relatedness_Value < 0.40)
full_sibs <-filter(hs2, Relatedness_Value > 0.40)
all_sibs <-filter(hs2, Relatedness_Value > halfsibs)
all_count <-ddply(hs2, ConYear1~ConYear2, summarize, ALL = n_distinct(pair.no))
hs_count <-ddply(half_sibs, ConYear1~ConYear2, summarize, HS = n_distinct(pair.no)) 
fs_count <-ddply(full_sibs, ConYear1~ConYear2, summarize, FS = n_distinct(pair.no))
as_count <-ddply(all_sibs, ConYear1~ConYear2, summarize, AS = n_distinct(pair.no))
sum_rel <-ddply(hs2, ConYear1~ConYear2, summarize, avr = mean(Relatedness_Value), sdrel =sd(Relatedness_Value),se= sem(Relatedness_Value),LCI = quantile(Relatedness_Value, 0.025), UCI= quantile(Relatedness_Value, 0.975))
sum_rel <-full_join(sum_rel, hs_count, by = c("ConYear1","ConYear2"))
sum_rel <-full_join(sum_rel, fs_count, by = c("ConYear1","ConYear2"))
sum_rel <-full_join(sum_rel, as_count, by = c("ConYear1","ConYear2"))
sum_rel <-full_join(sum_rel, all_count, by = c("ConYear1","ConYear2"))
sum_rel <-mutate(sum_rel, percentsib = round((HS/ALL)*1000,0),
                 ConYear1 = as.character(ConYear1), ConYear2=as.character(ConYear2)) #in a population of 1000 individuals, how many half siblings. 

CONYEAR <- ddply(hs, ~ConYear1, summarize, avr = mean(Relatedness_Value), sdrel =sd(Relatedness_Value),se= sem(Relatedness_Value),LCI = quantile(Relatedness_Value, 0.025), UCI= quantile(Relatedness_Value, 0.975))

#how many individuals in Mt SIGNIFICANTLY exceed the half sibling threshold?
filter(hs2, trioml.low > halfsibs) %>%n_distinct()
#how many individuals significantly exceed the full sibling mean?
filter(hs2, trioml.low > fullsibs) %>%n_distinct()

#Bay pairs sig diffs. 
hs2 <-unite(hs2, ConPair, ConYear1, ConYear2, sep="-", remove=FALSE)
anorel <-lm(Relatedness_Value~ConPair, na.action=na.omit, data=hs2)
ano2 <-car::Anova(anorel)
ano2
ano3<-df.residual(anorel)
MSerror<-deviance(anorel)/ano3
comparison <- HSD.test(anorel,c("ConPair"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/10, group=TRUE)
comparison

#heatmap
sum_rel <-filter(sum_rel, ConYear1 %in% c("1_2015","2_2015","3_adults","4_adults", "1_2016","2_2016"))%>%
  filter(ConYear2 %in% c("1_2015","2_2015","3_adults","4_adults", "1_2016","2_2016"))
  
sum_rel <-mutate(sum_rel, cy1 = pmin(ConYear1,ConYear2), cy2 =pmax(ConYear1,ConYear2)) %>% arrange(cy1)%>%
 mutate(ConYear1 = as.factor(ConYear1), ConYear2=as.factor(ConYear2))
cy1 <- fct_recode(sum_rel$ConYear1, c("1_2015","2_2015", "1_2016","2_2016","3_adults","4_adults"))
cy2 <- fct_recode(sum_rel$ConYear2, c("1_2015","2_2015", "1_2016","2_2016","3_adults","4_adults"))

ggplot(sum_rel, aes(cy1, cy2, label= AS))+ #switch between AS and percentsib
  geom_tile(aes(fill=avr), show.legend = TRUE)+
  scale_fill_viridis(option="A", direction = -1)+  #A is magma, D is regular viridis
  #scale_fill_gradient(low ="#ffeda0" , high = "#800026")+ #warm colors
  #scale_fill_gradient(low ="#023858" , high = "#ece7f2")+ #cool colors
   geom_text(color="white", size=6)+
  coord_fixed(ratio = 1) +
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 14),axis.title = element_text(size =14),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
ggsave('Mtrelheatmag.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)

#now for shinnecock. 
hs2 <-filter(hs, Bay1=="Shin" & Bay2=="Shin")
half_sibs <-filter(hs2, Relatedness_Value > halfsibs & Relatedness_Value < 0.40)
full_sibs <-filter(hs2, Relatedness_Value > 0.40)
all_sibs <-filter(hs2, Relatedness_Value > halfsibs)
all_count <-ddply(hs2, ConYear1~ConYear2, summarize, ALL = n_distinct(pair.no))
hs_count <-ddply(half_sibs, ConYear1~ConYear2, summarize, HS = n_distinct(pair.no)) 
fs_count <-ddply(full_sibs, ConYear1~ConYear2, summarize, FS = n_distinct(pair.no))
as_count <-ddply(all_sibs, ConYear1~ConYear2, summarize, AS = n_distinct(pair.no))
sum_rel <-ddply(hs2, ConYear1~ConYear2, summarize, avr = mean(Relatedness_Value), sdrel =sd(Relatedness_Value),se= sem(Relatedness_Value),LCI = quantile(Relatedness_Value, 0.025), UCI= quantile(Relatedness_Value, 0.975))
sum_rel <-full_join(sum_rel, hs_count, by = c("ConYear1","ConYear2"))
sum_rel <-full_join(sum_rel, fs_count, by = c("ConYear1","ConYear2"))
sum_rel <-full_join(sum_rel, as_count, by = c("ConYear1","ConYear2"))
sum_rel <-full_join(sum_rel, all_count, by = c("ConYear1","ConYear2"))
sum_rel <-mutate(sum_rel, percentsib = round((HS/ALL)*1000,0),
                 ConYear1 = as.character(ConYear1), ConYear2=as.character(ConYear2)) #in a population of 1000 individuals, how many half siblings. 
CONYEAR <- ddply(hs, ~ConYear1, summarize, avr = mean(Relatedness_Value), sdrel =sd(Relatedness_Value),se= sem(Relatedness_Value),LCI = quantile(Relatedness_Value, 0.025), UCI= quantile(Relatedness_Value, 0.975))
#heatmap
sum_rel <-filter(sum_rel, ConYear1 %in% c("1_2016","2_2016", "1_2017","2_2017"))%>%
  filter(ConYear2 %in% c("1_2016","2_2016", "1_2017","2_2017"))
sum_rel <-mutate(sum_rel, cy1 = pmin(ConYear1,ConYear2), cy2 =pmax(ConYear1,ConYear2)) %>% arrange(cy1)%>%
  mutate(ConYear1 = as.factor(ConYear1), ConYear2=as.factor(ConYear2))
cy1 <- fct_recode(sum_rel$ConYear1, c("1_2016","2_2016", "1_2017","2_2017"))
cy2 <- fct_recode(sum_rel$ConYear2, c("1_2016","2_2016", "1_2017","2_2017"))
ggplot(sum_rel, aes(cy1, cy2, label= AS))+ #switch between AS and percentsib
  geom_tile(aes(fill=avr), show.legend = TRUE)+
  #scale_fill_viridis(option="D", direction = -1)+  #A is magma, D is regular viridis
  scale_fill_gradient(low ="#ffeda0" , high = "#800026")+ #warm colors
  #scale_fill_gradient(low ="#023858" , high = "#ece7f2")+ #cool colors
  geom_text(color="white", size=6)+
  coord_fixed(ratio = 1) +
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 14),axis.title = element_text(size =14),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
ggsave('Shinrelheat.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)


#Bay pairs sig diffs. 
hs2 <-unite(hs2, ConPair, ConYear1, ConYear2, sep="-", remove=FALSE)
anorel <-lm(Relatedness_Value~ConPair, na.action=na.omit, data=hs2)
ano2 <-car::Anova(anorel)
ano2
ano3<-df.residual(anorel)
MSerror<-deviance(anorel)/ano3
comparison <- HSD.test(anorel,c("ConPair"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/10, group=TRUE)
comparison





#who is producing the half sibs???
# network plot??? 
#https://github.com/sctyner/geomnet
#for now, 
stack_ind <-unique(c(all_sibs$ind1.id, all_sibs$ind2.id)) 

count_pairs1 <-ddply(all_sibs, ~ind1.id, summarize, pair.num =n_distinct(pair.no)) %>% dplyr::rename(INDV = ind1.id)
count_pairs2 <-ddply(all_sibs, ~ind2.id, summarize, pair.num =n_distinct(pair.no)) %>% dplyr::rename(INDV = ind2.id)
stack_pairs <-bind_rows(count_pairs1,count_pairs2) %>%arrange(INDV)
mean(stack_pairs$pair.num)
sd(stack_pairs$pair.num)
max(stack_pairs$pair.num)
sum(stack_pairs$pair.num)

#correlation with inbreeding. 
inbreed <-unite(inbreed, ConYear, Con, Year)
in_rel <-ddply(inbreed, ~ConYear, summarize, mean.L3 = mean(L3), mean.LH=mean(LH))
#idk which one to use, so we'll use L3
in_rel <-mutate(in_rel, ConYear1=ConYear, ConYear2=ConYear)
sum_rel2 <-left_join(sum_rel, in_rel, by=c("ConYear1")) %>% dplyr::select(-ConYear2.y)%>%
  dplyr::rename(L3.cy1 = mean.L3, LH.cy1 =mean.LH, ConYear2 = ConYear2.x)
sum_rel2 <-left_join(sum_rel2, in_rel, by=c("ConYear2")) %>% dplyr::select(-ConYear.x, -ConYear.y, -ConYear1.y)%>%
  dplyr::rename(L3.cy2 = mean.L3, LH.cy2 =mean.LH, ConYear1 = ConYear1.x)

cor.test(sum_rel2$percentsib, sum_rel2$avr,  method = "pearson", exact=TRUE)
cor.test(sum_rel2$percentsib, sum_rel2$L3.cy1,method = "pearson", exact=TRUE) #inbreeding of first conyear to percent siblings
cor.test(sum_rel2$percentsib, sum_rel2$L3.cy2,method = "pearson", exact=TRUE) #inbreeding of second conyear to percent siblings




#mean relatedness value of same bay pairs
ddply(same.bay, ~ConYear1, summarize, avr = mean(Relatedness_Value), sdrel =sd(Relatedness_Value),se= sem(Relatedness_Value),LCI = quantile(Relatedness_Value, 0.025), UCI= quantile(Relatedness_Value, 0.975))
n_distinct(same.bay)
filter(same.bay,Relatedness_Value >= halfsibs & Relatedness_Value <= fullsibs) %>% n_distinct() #how many half sibs
filter(same.bay,Relatedness_Value > fullsibs) %>% n_distinct() #how many full sibs
filter(same.bay,Relatedness_Value > parentoffspring) %>% n_distinct()


ddply(diff.bay, ConYear2~ConYear1, summarize, avr = mean(Relatedness_Value), sdrel =sd(Relatedness_Value),se= sem(Relatedness_Value),LCI = quantile(Relatedness_Value, 0.025), UCI= quantile(Relatedness_Value, 0.975))
n_distinct(diff.bay)
filter(diff.bay,Relatedness_Value >= halfsibs & Relatedness_Value <= fullsibs) %>% n_distinct() #how many half sibs
filter(diff.bay,Relatedness_Value > fullsibs) %>% n_distinct() #how many full sibs
filter(diff.bay,Relatedness_Value > parentoffspring) %>% n_distinct() #how many full sibs








## anova to compare relatedness and ingroup outgroup
library("agricolae")
anorel <-lm(Relatedness_Value~Bay1*io, na.action=na.omit, data=inout)
ano2 <-car::Anova(anorel)
ano2
ano3<-df.residual(anorel)
MSerror<-deviance(anorel)/ano3
comparison <- HSD.test(anorel,c("io","Bay1"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05, group=TRUE)
comparison

## anova to compare relatedness and ingroup outgroup
library("agricolae")
anorel <-lm(Relatedness_Value~ioo, na.action=na.omit, data=inout)
ano2 <-car::Anova(anorel)
ano2
ano3<-df.residual(anorel)
MSerror<-deviance(anorel)/ano3
comparison <- HSD.test(anorel,"ioo",MSerror=MSerror, unbalanced=TRUE,alpha=0.05, group=TRUE)
comparison

# same bays
anorel <-lm(Relatedness_Value~Con1*Year1, na.action=na.omit, data=same.bay)
ano2 <-car::Anova(anorel)
ano2
ano3<-df.residual(anorel)
MSerror<-deviance(anorel)/ano3
comparison <- HSD.test(anorel,c("Con1","Year1"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05, group=TRUE)
comparison

#diff bay
anorel <-lm(Relatedness_Value~ConYear1, na.action=na.omit, data=diff.bay)
ano2 <-car::Anova(anorel)
ano2
ano3<-df.residual(anorel)
MSerror<-deviance(anorel)/ano3
comparison <- HSD.test(anorel,c("ConYear1"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05, group=TRUE)
comparison

#pairwise tests within each bay between same bay and diff bay
pairwise.t.test(same.bay$Relatedness_Value,same.bay$Bay1, pool.sd=FALSE)


########################################################################################
# write inbreeding to file
inbreed <- relatedness_triad$inbreeding %>%
  dplyr::select(ind.id, L3, LH) %>%
  dplyr::rename(INDV = ind.id) %>% mutate(INDV=as.factor(INDV))
write_delim(inbreed, "inbreeding")

inbreed.ci <-relatedness_triad$inbreeding.ci95%>%
  dplyr::select(ind.id, L3.low, L3.high, LH.low, LH.high, LR.low, LR.high) %>%
  dplyr::rename(INDV = ind.id)%>%mutate(INDV=as.factor(INDV))
#this is problematic because the individual ID doesn't match up. But we won't worry because I don't think we will use the CI> 

######Inbreeding########
# inbreeding  - you can put the values from the simulation in later. 
ggplot(inbreed, aes(x = L3)) +
  geom_histogram(binwidth = 0.005, color = "light grey", fill = "dark grey") +
  geom_vline(aes(xintercept = mean(LH, na.rm = TRUE)), #dark blue is the mean
             color = "black", linetype = "dashed", size = 0.5) +
  theme_cowplot()+
  labs(x = "inbreeding coefficient (Fis)", y = "individuals") 
ggsave('inbreedldyoyonly.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")

#inbreeding by bay/year 
popinfo <-dplyr::select(wfpop4df, Ind, Pop) %>% separate(Pop, c("Ocean","Bay","Con","Year"), remove=FALSE)
wf.df <-mutate(popinfo,ind1.id=Ind,ind2.id=Ind)
wf.df <-dplyr::rename(wf.df, INDV=Ind)
inbreed <-left_join(inbreed, wf.df, by="INDV") %>% unite(BayYear,Bay,Year,sep=" ", remove=FALSE)

inbreed %>%
  ggplot(aes(x=fct_rev(BayYear),y=L3))+ 
  #geom_boxplot(aes(fill=Bay))+ coord_flip()+
  ggplot2::stat_summary(fun.data = mean_sdl,fun.args = list(mult = 1),geom = "pointrange",position = ggplot2::position_nudge(x = 0.05, y = 0)) +coord_flip()+ 
  geom_flat_violin(aes(fill=Bay),position = position_nudge(x = .1, y = 0),adjust=2, trim = FALSE)+
  geom_hline(aes(yintercept = mean(L3, na.rm = TRUE)),color = "black", linetype = "dashed", size = 0.5) +
  scale_fill_manual(name = "Bay",values = drabcolors)+
  xlab(' ')+ylab("Inbreeding Coefficient (Fis)")+ 
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('indvinbreedingbayear.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 8, height = 10)


## anova to compare inbreeding between bays. 
library("agricolae")
anoin <-lm(L3~BayYear, na.action=na.omit, data=inbreed)
ano2 <-car::Anova(anoin)
ano2
ano3<-df.residual(anoin)
MSerror<-deviance(anoin)/ano3
comparison <- HSD.test(anoin,c("Con","Year"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05, group=TRUE)
comparison





############ group relatedness ##############
#format data so that first two letters of IND denote group name, then into gdata format
wf.df <-dplyr::select(wf.df, -ind1.id, -ind2.id) 
#wf.df <-dplyr::rename(wf.df, Ind=INDV)
group.df <-left_join(df, wf.df, by="Ind") %>% 
  #mutate(pop.prefix = substr(Bay, 1,2)) %>% 
  #mutate(pop.prefix = ifelse(ConYear == "1_2016", "FE", ifelse(ConYear=="1_2017","SE",ifelse(ConYear=="2_2016","FL","SL"))))%>%
  mutate(pop.prefix = ifelse(ConYear == "1_2015", "FE", ifelse(ConYear=="1_2016","SE",ifelse(ConYear=="2_2015","FL",
                                  ifelse(ConYear=="2_2016","SL",ifelse(ConYear=="3_2016","AM",ifelse(ConYear=="3_2015","AM",ifelse(ConYear=="4_2015","AR","AR"))))))))%>%
  unite(name, pop.prefix, Ind, sep="_", remove =TRUE) %>% dplyr::select(-Pop, -Ocean, -Bay, -Con, -Year, -ConYear)
write.table(group.df, "scratch",row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) #write it.
genotypedata <- readgenotypedata("scratch")# import input file as list (gdata, nloci, nalleles, ninds, freqs)


#re-run overnight
relbay <- grouprel(genotypes = genotypedata$gdata, estimatorname = "trioml", usedgroups = "all", iterations= 100)



## make a version of this analysis done with the clean dataset and run overnight. 

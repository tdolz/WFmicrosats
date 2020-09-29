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

##### Formating the dataset #####
wfpop <- read.genalex("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_17_sept20204genalex_doubl0.csv")

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
  geom_boxplot() +
  facet_wrap(~ relationship)

#calculate correlation coefficient between observed values for each estimator and the expected values. 
urval <- rep(0, 100) 
hsval <- rep (0.25 , 100)
fsval <- rep (0.5 , 100)
poval <- rep (0.5 , 100)
relvals <- c(poval , fsval , hsval , urval)

cor(relvals , simrel[, 5])
cor(relvals , simrel[, 6])
cor(relvals , simrel[, 10])
cor(relvals , simrel[, 11])

#the simulation --REMEMBER TO TURN OFF AND ON
#change bootstrap back to 100 because 1000 just takes too long. 
cosim <-compareestimators(relatedness_triad, ninds=100)
cosim
#ggsave("relatendess_simulation.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()
#extract mean, median & CI values for relatedness from the simulation. 
relsim <- as.data.frame(cosim$data)
relply<- ddply(relsim, Estimator~relationship,summarize, mean.rel = mean(Relatedness_Value), LIrel = quantile(Relatedness_Value,0.05), HIrel = quantile(Relatedness_Value, 0.95), medianrel= quantile(Relatedness_Value, 0.5))
relply




#trioml estimate (Wang) - Can't do this because it crashes R studio. 
relatedness_triad <- coancestry(genotypedata$gdata,trioml = 1)
relatedT <- relatedness_triad$relatedness %>%
  dplyr::select(pair.no, ind1.id, ind2.id, trioml)

# write relatedness to file
relatedn <- relatedness_lynchrd$relatedness %>%
  dplyr::select(pair.no, ind1.id, ind2.id, lynchrd)

#not sure what's going on here. 
library("readr")
write_delim(relatedn, "pairwise_relatedness")
write_delim(relatedT, "pairwise_relatedness_trioml")

# write inbreeding to file
inbreed <- relatedness_triad$inbreeding %>%
  dplyr::select(ind.id, L3, LH) %>%
  dplyr::rename(INDV = ind.id)
write_delim(inbreed, "inbreeding")



#make your own boxplot of just lynch and ritland or whatever. 
relsim %>%
  filter(Estimator == "W") %>%
  ggplot(aes(x=fct_rev(relationship),y=Relatedness_Value))+ 
  geom_boxplot()+ 
  xlab('relationship ')+ylab("Pairwise Relatedness")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('simulationLR.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)










# relatedness (LYNCHRD)- remember to change these values depending on what dataset you are using. 
ggplot(relatedn, aes(x = lynchrd)) +
  geom_histogram(binwidth = 0.001, color = "light grey", fill = "light grey") +
  geom_vline(aes(xintercept = mean(lynchrd, na.rm = TRUE)),
             color = "black", linetype = "dashed", size = 0.5) +
  #geom_vline(aes(xintercept = quantile(lynchrd, 0.95)),
             #color = "darkred", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = 0.23), color = "darkblue", linetype = "dashed", size = 0.5) + #the estimated mean value for half sibs in the simulation
  geom_vline(aes(xintercept = 0.48), color = "darkred", linetype = "dashed", size = 0.5) + #the estimated mean value for full sibs in the simulation
  geom_vline(aes(xintercept = 0.47), color = "darkgreen", linetype = "dashed", size = 0.5) + #the estimated mean value for parent offspring in the simulation
  #geom_vline(aes(xintercept = 0.08),color = "darkblue", linetype = "dotted", size = 0.5) + #the Lynchrd 0.05 quantile on the simulated half sibs relationship
  #geom_vline(aes(xintercept = 0.49), color = "darkblue", linetype = "dotted", size = 0.5) + #red is the 0.95 quantile
  labs(x = "relatedness", y = "number of pairs")+
  theme_cowplot()
#ggsave("pairwiserelatenessLYNCHRDALL_20.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()

#attach strata information to the relatedness estimation. 
#h <-melt(relatedn, id=c("pair.no","ind1.id","ind2.id"))
#h <-dplyr::rename(h, Estimator=variable,Relatedness_Value=value)

#how many of each? (based on estimator mean)
filter(relatedn,lynchrd > 0.23 & lynchrd <= 0.49) %>% n_distinct() #how many half sibs
filter(relatedn, lynchrd > 0.48) %>% n_distinct() #full sibs or parent offspring
n_distinct(relatedn)
filter(relatedn, lynchrd <= 0.23) %>% n_distinct()
filter(relatedn, lynchrd < 0.08) %>% n_distinct()

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
ddply(same.bay, ~Bay1, summarize, avr = mean(Relatedness_Value), sdrel =sd(Relatedness_Value),LCI = quantile(Relatedness_Value, 0.025), UCI= quantile(Relatedness_Value, 0.975))
filter(same.bay,Relatedness_Value >= 0.23 & Relatedness_Value <= 0.48) %>% n_distinct() #how many half sibs
filter(same.bay,Relatedness_Value > 0.48) %>% n_distinct() #how many full sibs

#relatedness diff bays
diff.bay <-filter(hs, Bay1!=Bay2)
#mean relatedness value of same bay pairs
ddply(diff.bay, Bay2~Bay1, summarize, avr = mean(Relatedness_Value), sdrel =sd(Relatedness_Value),LCI = quantile(Relatedness_Value, 0.025), UCI= quantile(Relatedness_Value, 0.975))
filter(diff.bay,Relatedness_Value >= 0.23 & Relatedness_Value <= 0.48) %>% n_distinct() #how many half sibs
filter(diff.bay,Relatedness_Value > 0.48) %>% n_distinct() #how many full sibs

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
  geom_hline(aes(yintercept = 0.26), color = "black", linetype = "dotted", size = 0.7) + #the estimated mean value for half sibs in the simulation
  geom_hline(aes(yintercept = 0.50), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for full sibs in the simulation
  #geom_hline(aes(yintercept = 0.50), color = "darkgreen", linetype = "dashed", size = 0.7) + #the estimated mean value for parent offspring in the simulation
  xlab(' ')+ylab("Mean Relatedness")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('samebay_wAdults.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

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
  geom_hline(aes(yintercept = 0.26), color = "black", linetype = "dotted", size = 0.7) + #the estimated mean value for half sibs in the simulation
  geom_hline(aes(yintercept = 0.50), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for full sibs in the simulation
  #geom_hline(aes(yintercept = 0.50), color = "darkgreen", linetype = "dashed", size = 0.7) + #the estimated mean value for parent offspring in the simulation
  xlab(' ')+ylab("Mean Relatedness")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('diffbay_wAdults.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

#in groups and out groups
same.bay <-mutate(same.bay, io = "same")
diff.bay <-mutate(diff.bay, io = "diff")
inout <- bind_rows(same.bay, diff.bay)
inout <-unite(inout, ioo, Bay1, io, sep="_", remove=FALSE)

inout %>%
  ggplot(aes(x=fct_rev(ioo),y=Relatedness_Value))+ # bay pairs
  #geom_boxplot(aes(fill=Bay1))+ coord_flip()+
  ggplot2::stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "pointrange",
    position = ggplot2::position_nudge(x = 0.05, y = 0)
  ) +coord_flip()+ 
  geom_flat_violin(aes(fill=Bay1),position = position_nudge(x = .1, y = 0),adjust=2, trim = FALSE)+
  scale_fill_manual(name = "Bay",values = drabcolors)+
  geom_hline(aes(yintercept = 0.21), color = "black", linetype = "dotted", size = 0.7) + #the estimated mean value for half sibs in the simulation
  geom_hline(aes(yintercept = 0.07), color = "grey", linetype = "dotted", size = 0.7) + #the estimated mean value for half sibs in the simulation
  geom_hline(aes(yintercept = 0.49), color = "grey", linetype = "dotted", size = 0.7) + #the estimated mean value for half sibs in the simulation
  geom_hline(aes(yintercept = 0.47), color = "black", linetype = "dashed", size = 0.5) + #the estimated mean value for full sibs in the simulation
  xlab(' ')+ylab("Mean Relatedness")+ 
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('inandoutgroups.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 8, height = 10)
#it may be more appropriate to compare YOY only. 

##overlaying histograms of same and different bays, you can see they are not different. 
plot_multi_histogram <- function(df, feature, label_column) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
    geom_density(alpha=0.7) +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density")
  plt + guides(fill=guide_legend(title=label_column))
}
plot_multi_histogram(inout, 'Relatedness_Value', 'io')


######Inbreeding########
# inbreeding  - you can put the values from the simulation in later. 
ggplot(inbreed, aes(x = L3)) +
  geom_histogram(binwidth = 0.005, color = "light grey", fill = "light grey") +
  geom_vline(aes(xintercept = mean(LH, na.rm = TRUE)), #dark blue is the mean
             color = "black", linetype = "dashed", size = 0.5) +
  theme_cowplot()+
  labs(x = "inbreeding coefficient (Fis)", y = "individuals") 
ggsave('inbreedld.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")

#inbreeding by bay. 
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

## anova to compare relatedness
library("agricolae")
anorel <-lm(Relatedness_Value~io*Bay1, na.action=na.omit, data=inout)
ano2 <-car::Anova(anorel)
ano2
df<-df.residual(anorel)
MSerror<-deviance(anorel)/df
comparison <- HSD.test(anorel,c("io","Bay1"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05, group=TRUE)
comparison

## anova to compare inbreeding
anoin <-lm(L3~Bay*Con, na.action=na.omit, data=inbreed)
ano2 <-car::Anova(anoin)
ano2
df<-df.residual(anoin)
MSerror<-deviance(anoin)/df
comparison <- HSD.test(anoin,c("Bay","Con"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05, group=TRUE)
comparison


### Remove individuals that had super high relatedness, are they the same individual? think about removing them. 
## I am not going to remove them because they're from the same bay and the'yre not out of the range of full sibs from the simulation (which is also from the data so..idk)

############### Now Do without the mattituck adults ###################
#First remove the mattituck adults, because you don't want potential parents messing up the siblings relationship. 
popNames(wfpopLD)
setPop(wfpopLD) <-~Bay/Con/Year
wfyoy <-popsub(wfpopLD, blacklist=c("Atl_Mt_3_2015","Atl_Mt_4_2015","Atl_Mt_5_2015","Atl_Mt_3_2016","Atl_Mt_4_2016","Atl_Mt_5_2016"))


#convert to the right data format
setPop(wfyoy) <-~Bay
df <-genind2df(wfyoy, usepop = FALSE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(33,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32)]
df[df=="NA"] <- 0 # missing data must be 0
write.table(df, "scratchyoy",row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) #write it.
genotypedata <- readgenotypedata("scratchyoy")# import input file as list (gdata, nloci, nalleles, ninds, freqs)

# pairwise relatedness (Lynch & Ritland 1999)
relatedness_lynchrd <- coancestry(genotypedata$gdata,
                                  lynchrd = 1)

# write relatedness to file
relatedn <- relatedness_lynchrd$relatedness %>%
  dplyr::select(pair.no, ind1.id, ind2.id, lynchrd)

#not sure what's going on here. 
library("readr")
write_delim(relatedn, "pairwise_relatedness")
 #write inbreeding to file
inbreed <- relatedness_lynchrd$inbreeding %>%
  dplyr::select(ind.id, L3, LH) %>%
  dplyr::rename(INDV = ind.id)
write_delim(inbreed, "inbreeding")

#the simulation
#change bootstrap back to 100 because 1000 just takes too long. 
cosim <-compareestimators(relatedness_lynchrd, ninds=100)
cosim
ggsave("relatendess_simulation_yoyonly.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
dev.off()


#extract mean, median & CI values for relatedness from the simulation. 
relsim <- as.data.frame(cosim$data)
relply<- ddply(relsim, Estimator~relationship,summarize, mean.rel = mean(Relatedness_Value), LIrel = quantile(Relatedness_Value,0.05), HIrel = quantile(Relatedness_Value, 0.95), medianrel= quantile(Relatedness_Value, 0.5))
relply

# relatedness - remember to change these values depending on what dataset you are using. 
ggplot(relatedn, aes(x = lynchrd)) +
  geom_histogram(binwidth = 0.001, color = "light grey", fill = "light grey") +
  geom_vline(aes(xintercept = mean(lynchrd, na.rm = TRUE)),
             color = "black", linetype = "dashed", size = 0.5) +
  #geom_vline(aes(xintercept = quantile(lynchrd, 0.95)),
  #color = "darkred", linetype = "dashed", size = 0.5) +
  geom_vline(aes(xintercept = 0.23), color = "darkblue", linetype = "dashed", size = 0.7) + #the estimated mean value for half sibs in the simulation
  geom_vline(aes(xintercept = 0.48), color = "darkred", linetype = "dashed", size = 0.7) + #the estimated mean value for full sibs in the simulation
  #geom_vline(aes(xintercept = 0.49), color = "darkgreen", linetype = "dashed", size = 0.7) + #the estimated mean value for parent offspring in the simulation
  #geom_vline(aes(xintercept = 0.08),color = "darkblue", linetype = "dotted", size = 0.5) + #the Lynchrd 0.05 quantile on the simulated half sibs relationship
  #geom_vline(aes(xintercept = 0.49), color = "darkblue", linetype = "dotted", size = 0.5) + #red is the 0.95 quantile
  labs(x = "relatedness", y = "number of pairs")+
  theme_cowplot()
ggsave("pairwiserelatenessALL_20_yoy.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
dev.off()

#attach strata information to the relatedness estimation. 
h <-melt(relatedn, id=c("pair.no","ind1.id","ind2.id"))
h <-dplyr::rename(h, Estimator=variable,Relatedness_Value=value)

#how many of each? (based on estimator mean)
filter(h,Relatedness_Value > 0.23 & Relatedness_Value <= 0.48) %>% n_distinct()
filter(h,Relatedness_Value > 0.48) %>% n_distinct() 

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
ddply(same.bay, ~Bay1, summarize, avr = mean(Relatedness_Value),sdrel =sd(Relatedness_Value),LCI = quantile(Relatedness_Value, 0.025), UCI= quantile(Relatedness_Value, 0.975))
filter(same.bay,Relatedness_Value >= 0.23 & Relatedness_Value <= 0.48) %>% n_distinct() 
filter(same.bay,Relatedness_Value > 0.48) %>% n_distinct()  
  #very high relatedness value between 14 & 18 in Napeague and between S55 and S58 in Shinnecock. You should probably remove those individuals. 


#mean relatedness value of different bay pairs
diff.bay <-filter(hs, Bay1!=Bay2)
ddply(diff.bay, Bay1~Bay2, summarize, avr = mean(Relatedness_Value),sdrel =sd(Relatedness_Value),LCI = quantile(Relatedness_Value, 0.025), UCI= quantile(Relatedness_Value, 0.975))
filter(diff.bay,Relatedness_Value >= 0.23 & Relatedness_Value < 0.48) %>% n_distinct() 
filter(diff.bay,Relatedness_Value > 0.48) %>% n_distinct() 

#Shinnecock only
shin.only <-filter(hs, Bay1=="Shin" & Bay2=="Shin")
shin.halfsibs <-filter(shin.only, Relatedness_Value >= 0.23 & Relatedness_Value < 0.48) #half sibs. 
shin.fullsibs <-filter(shin.only,  Relatedness_Value >= 0.48) #full sibs. 
shin.sibs <-filter(shin.only, Relatedness_Value >= 0.23) #half or full
all_ssibs <- as.vector(rbind(shin.sibs$ind1.id, shin.sibs$ind2.id))
n_distinct(all_ssibs) #how many individuals involved in at least one pair  #37, so some are involved in multiple pairs. 
shinsibs <-unique(all_ssibs) #which individuals. 

#could make a tile plot visualizing all the individuals in a bay and their pairwise relatedness? 
## unique color ramp if they are siblings, half siblings.  # we could always change out the estimator later. 
## could also make plots for diff bay and same bay pairs, with a column called IND BAY YEAR as the identifier? 
# this would at the very least allow you to visualize problematic individuals. 

cols <- c("(-0.13,0]"="#fff7ec", "(0,0.27]" = "#fdd49e", "(0.27,0.49]" = "#fc8d59", "(0.49,1]" = "#990000")
#shin.only <- unite(shin.only, ind_bay, c("bayyear1", "bayyear2"),sep="-")
#shin.only<-mutate(shin.only, p_if_sig=ifelse(sigp=="sig",p.value,NA), starsig=ifelse(sigp=="sig","*",NA)) 
shin.only <- shin.only %>% arrange(Pop.x, Pop.y)
shin.only<-mutate(shin.only, relbreaks = cut(Relatedness_Value, breaks=c(-0.13,0, 0.27,0.49,1)))
ggplot(shin.only, aes(ind1.id,ind2.id, label=round(Relatedness_Value),3))+
  geom_tile(aes(fill=relbreaks), show.legend = TRUE)+
  scale_fill_manual(values=cols)+
  geom_text(color="white", size=1)+
  #facet_grid(Pop.x~Pop.y, scales="free")+
  xlab("")+ylab("")+
  theme(text=element_text(size=8), axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
#ggsave("shinonly11_tiles_facet.png", path="/Users/tdolan/documents/WF SK PROJ/Survey data/Field Survey Paper/final figures")
#ggsave("shinonly11_tiles.png", path="/Users/tdolan/documents/WF SK PROJ/Survey data/Field Survey Paper/final figures")
#dev.off()

#all of them, by bay. 
hs<-mutate(hs, relbreaks = cut(Relatedness_Value, breaks=c(-0.13,0, 0.27,0.49,1)))
ggplot(hs, aes(ind1.id,ind2.id, label=round(Relatedness_Value),3))+
  geom_tile(aes(fill=relbreaks), show.legend = TRUE)+
  scale_fill_manual(values=cols)+
  geom_text(color="white", size=1)+
  facet_grid(Bay1~Bay2, scales="free")+
  xlab("")+ylab("")+
  theme(text=element_text(size=8), axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
#ggsave("baypairs11_tiles_facet.png", path="/Users/tdolan/documents/WF SK PROJ/Survey data/Field Survey Paper/final figures")
#dev.off()

#let's go back to the diff bay pair barplot scenario.
diff.bay <-tidyr::unite(diff.bay,col="bay.pair", c("Bay1","Bay2"),sep="_",remove=FALSE )
diff.bay <-mutate(diff.bay, bay.pair =ifelse(bay.pair=="Jam_Nap" | bay.pair=="Nap_Jam", "Jam:Nap",bay.pair))
diff.bay <-mutate(diff.bay, bay.pair =ifelse(bay.pair=="Shin_Nap"| bay.pair=="Nap_Shin","Shin:Nap",bay.pair))
diff.bay <-mutate(diff.bay, bay.pair =ifelse(bay.pair=="Mor_Nap"| bay.pair=="Nap_Mor","Mor:Nap",bay.pair))
diff.bay <-mutate(diff.bay, bay.pair =ifelse(bay.pair=="Jam_Mor"|bay.pair=="Mor_Jam","Jam:Mor",bay.pair))
diff.bay <-mutate(diff.bay, bay.pair =ifelse(bay.pair=="Shin_Mor"|bay.pair=="Mor_Shin","Shin:Mor",bay.pair))
diff.bay <-mutate(diff.bay, bay.pair =ifelse(bay.pair=="Shin_Jam"| bay.pair=="Jam_Shin","Shin:Jam",bay.pair))
diff.bay <-mutate(diff.bay, as.factor(bay.pair))

#what's the mean out group relatedness
ddply(diff.bay, ~bay.pair, summarize, avr = mean(Relatedness_Value))

#half siblings
hs.bays <-filter(diff.bay, Relatedness_Value >= 0.27 & Relatedness_Value <= 0.49)
#how many pairs
ddply(hs.bays, Bay1~Bay2, summarize, num_halfsibs = dplyr::n_distinct(pair.no))
all_sibs <- as.vector(rbind(hs.bays$ind1.id, hs.bays$ind2.id))
n_distinct(all_sibs) #103 individuals involved in at least one cross-bay half sib pair
diffbayhalfsibs <-unique(all_sibs) #which individuals. 

#what about a heatmap of these idiots? - the spread is very interesting. 
#hs<-mutate(hs, relbreaks = cut(Relatedness_Value, breaks=c(-0.13,0, 0.27,0.49,1)))
ggplot(hs.bays, aes(ind1.id,ind2.id, label=round(Relatedness_Value),3))+
  #geom_tile(aes(fill=relbreaks), show.legend = TRUE)+
  geom_tile(aes(fill=Relatedness_Value), show.legend = TRUE)+
 # scale_fill_manual(values=cols)+
  geom_text(color="white", size=1)+
  facet_grid(Bay1~Bay2, scales="free")+
  xlab("")+ylab("")+
  theme(text=element_text(size=8), axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
#ggsave("diffbaypairs11_tiles_facet.png", path="/Users/tdolan/documents/WF SK PROJ/Survey data/Field Survey Paper/final figures")
#dev.off()

#look at potential full sibs in different bays: There are none. Which is good!
fs.bays <-filter(diff.bay, Relatedness_Value >= 0.49)
#how many pairs
ddply(fs.bays, Bay1~Bay2, summarize, num_halfsibs = dplyr::n_distinct(pair.no))
fall_sibs <- as.vector(rbind(fs.bays$ind1.id, fs.bays$ind2.id))
n_distinct(fall_sibs) #how many individuals involved in at least one pair
diffbayfullsibs <-unique(fall_sibs) #which individuals. 

# inbreeding  - you can put the values from the simulation in later. 
ggplot(inbreed, aes(x = L3)) +
  geom_histogram(binwidth = 0.005, color = "black", fill = "darkorange") +
  geom_vline(aes(xintercept = mean(LH, na.rm = TRUE)), #dark blue is the mean
             color = "darkblue", linetype = "dashed", size = 1) +
  theme_cowplot()+
  #geom_vline(aes(xintercept = 0.22), color = "green", linetype = "dashed", size = 1) + #the estimated mean value for half sibs in the simulation
  #geom_vline(aes(xintercept = 0.48), color = "green", linetype = "dashed", size = 1) + #the estimated mean value for full sibs in the simulation
  #geom_vline(aes(xintercept = 0.5), color = "green", linetype = "dashed", size = 1) + #the estimated mean value for parent offspring in the simulation
  #geom_vline(aes(xintercept = 0.04),color = "purple", linetype = "dashed", size = 0.5) + #the Lynchrd 0.05 quantile on the simulated half sibs relationship
  #geom_vline(aes(xintercept = quantile(LH, 0.95)), color = "darkred", linetype = "dashed", size = 1) + #red is the 0.95 quantile
  labs(x = "inbreeding coefficient (Fis)", y = "individuals") 
#ggsave('inbreedld.png', width = 7, height = 7)




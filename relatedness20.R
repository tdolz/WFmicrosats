### Relatedness #####

## NOTES ##
## includes 16 loci
## for no half 0 genotypes, use the double0 version of the files 
## population corrected for new cohort assignment
## based on the script "skreportdnaJAN_28_2020MATTITUCK.R"
#install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)
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

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")

##### Formating the dataset #####
wfpop <- read.genalex("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_20_sept20204genalex_doubl0.csv")

splitStrata(wfpop) <-~Ocean/Bay/Con/Year
setPop(wfpop) <-~Bay

### Format Data as in Microsats_20.R
wfpopLD <-genclone2genind(wfpop)
all_loci <- locNames(wfpopLD)# create vector of all loci
removeloc <- c("WF27","WF06","WF32")# create vector containing loci to remove
keeploc <- setdiff(all_loci, removeloc)# create list of loci to keep
wfpopLD <- wfpopLD[loc = keeploc]# filter loci in genind object
#now re-remove the missing data. 
wfpopLD <-wfpopLD %>% missingno("loci", cutoff=0.30) # removes loci where overall missingness > 30%, so WF01
wfpopLD <- wfpopLD %>% missingno("geno", cutoff=0.25) # remove samples that are > 25% missing
length(locNames(wfpopLD))# check number of loci in genind obj

######## Pairwise Relatedness (SIBLINGS ONLY) ######

#There are many different estimators for relatedness. In previous scripts we have tried more robust ways. 
# for now, Shannon's way. https://gist.github.com/sjoleary/3efd4a7d56b115fad319781298765a31

#convert to the right data format
setPop(wfpopLD) <-~Bay
df <-genind2df(wfpopLD, usepop = FALSE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(33,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32)]
df[df=="NA"] <- 0 # missing data must be 0
write.table(df, "scratch",row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) #write it.
genotypedata <- readgenotypedata("scratch")# import input file as list (gdata, nloci, nalleles, ninds, freqs)

# pairwise relatedness (Lynch & Ritland 1999)
relatedness_lynchrd <- coancestry(genotypedata$gdata,lynchrd = 2)

#trioml estimate (Wang) - Can't do this because it crashes R studio. 
#relatedness_triad <- coancestry(genotypedata$gdata,trioml = 2)

# write relatedness to file
relatedn <- relatedness_lynchrd$relatedness %>%
  dplyr::select(pair.no, ind1.id, ind2.id, lynchrd)

#not sure what's going on here. 
library("readr")
write_delim(relatedn, "pairwise_relatedness")
# write inbreeding to file
inbreed <- relatedness_lynchrd$inbreeding %>%
  dplyr::select(ind.id, L3, LH) %>%
  dplyr::rename(INDV = ind.id)
write_delim(inbreed, "inbreeding")

#the simulation
#change bootstrap back to 100 because 1000 just takes too long. 
cosim <-compareestimators(relatedness_lynchrd, ninds=100)
cosim
ggsave("relatendess_simulation.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
dev.off()


#extract mean, median & CI values for relatedness from the simulation. 
relsim <- as.data.frame(cosim$data)
relply<- ddply(relsim, Estimator~relationship,summarize, mean.rel = mean(Relatedness_Value), LIrel = quantile(Relatedness_Value,0.05), HIrel = quantile(Relatedness_Value, 0.95), medianrel= quantile(Relatedness_Value, 0.5))
relply

# relatedness - remember to change these values depending on what dataset you are using. 
ggplot(relatedn, aes(x = lynchrd)) +
  geom_histogram(binwidth = 0.001, color = "black", fill = "darkorange") +
  geom_vline(aes(xintercept = mean(lynchrd, na.rm = TRUE)),
             color = "darkblue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = quantile(lynchrd, 0.95)),
             color = "darkred", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = 0.23), color = "green", linetype = "dashed", size = 1) + #the estimated mean value for half sibs in the simulation
  geom_vline(aes(xintercept = 0.48), color = "green", linetype = "dashed", size = 1) + #the estimated mean value for full sibs in the simulation
  geom_vline(aes(xintercept = 0.47), color = "green", linetype = "dashed", size = 1) + #the estimated mean value for parent offspring in the simulation
  geom_vline(aes(xintercept = 0.03),color = "purple", linetype = "dashed", size = 0.5) + #the Lynchrd 0.05 quantile on the simulated half sibs relationship
  geom_vline(aes(xintercept = 0.53), color = "darkred", linetype = "dashed", size = 1) + #red is the 0.95 quantile
  labs(x = "relatedness", y = "number of pairs")+
  theme_cowplot()
#ggsave("pairwiserelatenessALL_11.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()


#attach strata information to the relatedness estimation. 
h <-melt(relatedn, id=c("pair.no","ind1.id","ind2.id"))
h <-dplyr::rename(h, Estimator=variable,Relatedness_Value=value)

#how many of each? (based on estimator mean)
filter(h,Relatedness_Value > 0.23 & Relatedness_Value <= 0.49) %>% n_distinct() #how many half sibs
filter(h,Relatedness_Value > 0.47) %>% n_distinct() #full sibs or parent offspring

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
ddply(same.bay, ~Bay1, summarize, avr = mean(Relatedness_Value))
filter(same.bay,Relatedness_Value >= 0.27 & Relatedness_Value <= 0.49) %>% n_distinct() #how many half sibs, 27
filter(same.bay,Relatedness_Value > 0.49) %>% n_distinct() #how many full sibs, 2, 
#very high relatedness value between 14 & 18 in Napeague and between S55 and S58 in Shinnecock. You should probably remove those individuals. 

### Remove individuals that had super high relatedness, are they the same individual? think about removing them. 
## I am not going to remove them because they're from the same bay and the'yre not out of the range of full sibs from the simulation (which is also from the data so..idk)

######Inbreeding########
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




############### Now Do without the mattituck adults ###################
#First remove the mattituck adults, because you don't want potential parents messing up the siblings relationship. 
popNames(wfpopLD)
setPop(wfpopLD) <-~Bay/Con/Year
wfyoy <-popsub(wfpopLD, blacklist=c("Atl_Mt_3_2015","Atl_Mt_4_2015","Atl_Mt_5_2015","Atl_Mt_3_2016","Atl_Mt_4_2016","Atl_Mt_5_2016"))



#convert to the right data format
setPop(wfyoy) <-~Bay
df <-genind2df(wfyoy, usepop = FALSE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(23,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)]
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
# write inbreeding to file
inbreed <- relatedness_lynchrd$inbreeding %>%
  dplyr::select(ind.id, L3, LH) %>%
  dplyr::rename(INDV = ind.id)
write_delim(inbreed, "inbreeding")

#the simulation
#change bootstrap back to 100 because 1000 just takes too long. 
cosim <-compareestimators(relatedness_lynchrd, ninds=100)
cosim
#ggsave("relatendess_simulation_yoyonly.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()


#extract mean, median & CI values for relatedness from the simulation. 
relsim <- as.data.frame(cosim$data)
relply<- ddply(relsim, Estimator~relationship,summarize, mean.rel = mean(Relatedness_Value), LIrel = quantile(Relatedness_Value,0.05), HIrel = quantile(Relatedness_Value, 0.95), medianrel= quantile(Relatedness_Value, 0.5))
relply

# relatedness - remember to change these values depending on what dataset you are using. 
ggplot(relatedn, aes(x = lynchrd)) +
  geom_histogram(binwidth = 0.001, color = "black", fill = "darkorange") +
  geom_vline(aes(xintercept = mean(lynchrd, na.rm = TRUE)),
             color = "darkblue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = quantile(lynchrd, 0.95)),
             color = "darkred", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = 0.27), color = "green", linetype = "dashed", size = 1) + #the estimated mean value for half sibs in the simulation
  geom_vline(aes(xintercept = 0.49), color = "green", linetype = "dashed", size = 1) + #the estimated mean value for full sibs in the simulation
  geom_vline(aes(xintercept = 0.04),color = "purple", linetype = "dashed", size = 0.5) + #the Lynchrd 0.05 quantile on the simulated half sibs relationship
  geom_vline(aes(xintercept = 0.57), color = "darkred", linetype = "dashed", size = 1) + #red is the 0.95 quantile
  labs(x = "relatedness", y = "number of pairs")+
  theme_cowplot()
#ggsave("pairwiserelatenessALL_11_yoy.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()

#attach strata information to the relatedness estimation. 
h <-melt(relatedn, id=c("pair.no","ind1.id","ind2.id"))
h <-dplyr::rename(h, Estimator=variable,Relatedness_Value=value)

#how many of each? (based on estimator mean)
filter(h,Relatedness_Value > 0.27 & Relatedness_Value <= 0.49) %>% n_distinct() #how many half sibs, 95
filter(h,Relatedness_Value > 0.49) %>% n_distinct() #how many full sibs, 2

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
ddply(same.bay, ~Bay1, summarize, avr = mean(Relatedness_Value))
filter(same.bay,Relatedness_Value >= 0.27 & Relatedness_Value <= 0.49) %>% n_distinct() #how many half sibs, 27 pairs
filter(same.bay,Relatedness_Value > 0.49) %>% n_distinct() #how many full sibs, 2, 
  #very high relatedness value between 14 & 18 in Napeague and between S55 and S58 in Shinnecock. You should probably remove those individuals. 


#mean relatedness value of different bay pairs
diff.bay <-filter(hs, Bay1!=Bay2)
ddply(diff.bay, ~Bay1, summarize, avr = mean(Relatedness_Value))
filter(diff.bay,Relatedness_Value >= 0.27 & Relatedness_Value < 0.49) %>% n_distinct() #how many half sibs, 69
filter(diff.bay,Relatedness_Value > 0.49) %>% n_distinct() #how many full sibs 0

#Shinnecock only
shin.only <-filter(hs, Bay1=="Shin" & Bay2=="Shin")
shin.halfsibs <-filter(shin.only, Relatedness_Value >= 0.27 & Relatedness_Value < 0.49) #half sibs. 
shin.fullsibs <-filter(shin.only,  Relatedness_Value >= 0.49) #full sibs. 
shin.sibs <-filter(shin.only, Relatedness_Value >= 0.27) #half or full
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




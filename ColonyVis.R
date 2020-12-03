### visualizing colony data
#Oct 10, 2020

#notes
#for next time 
#https://rdrr.io/github/thierrygosselin/radiator/man/write_colony.html

library('plyr')
library("dplyr")
library("tidyr")
library('purrr')
library("ggplot2")
library("lattice")
library("related")
library("cowplot")
library('reshape2')
library('forcats')

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")

#toolpack
sem <-function(x) sd(x)/sqrt(length(x))

#data
clust <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/colonyclust.csv", header = TRUE)
halfs <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/colonyhalfsibs.csv", header = TRUE)
wfpop4df <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_17_sept2020_doubl0ABC.csv", header = TRUE) #csv version 
rel_mt <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/relmt.csv",header = TRUE)
rel_shin <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/relshin.csv",header = TRUE)

#join pop info to clusts 
popinfo <-dplyr::select(wfpop4df, Ind, Pop) %>% separate(Pop, c("Ocean","Bay","Con","Year"), remove=FALSE)
clust <-left_join(clust, popinfo, by=c("Ind","Bay"))

#separate the runs
clust17 <-filter(clust, run=="17clust") %>% dplyr::select(-run)
clustHWE <-filter(clust, run =="HWEclust") %>% dplyr::select(-run)
halfs17 <-filter(halfs, run =="17clust") %>% dplyr::select(-run)
halfsHWE <-filter(halfs, run =="HWE") %>% dplyr::select(-run)
#do this for 17 first

#join pop info to clust
wf.df <-mutate(popinfo,ind1=Ind,ind2=Ind) %>% unite(ConYear, Con, Year, remove=FALSE)
fam <-dplyr::select(clust17, Ind, clust)
wf.df <-left_join(wf.df, fam, by=c("Ind"))
wf.df <-dplyr::rename(wf.df, clust1=clust)
fam <-dplyr::rename(fam, ind2=Ind)
clust17new <-dplyr::select(wf.df, -ind1, -ind2)
#join clust & pop to halfsibs
wf.df <-left_join(wf.df, fam, by=c("ind2")) %>% dplyr::rename(clust2=clust)
hs <-left_join(halfs17,wf.df, by=c("ind1","Bay"))
hs <-dplyr::rename(hs, Bay1=Bay,Con1=Con,Year1=Year, ConYear1 = ConYear)
hs <-dplyr::select(hs,-ind2.y) %>% dplyr::rename(ind2 =ind2.x)
hs <-left_join(hs,wf.df,by=c("ind2")) %>% dplyr::rename(Bay2=Bay,Con2=Con,Year2=Year,ConYear2 = ConYear)
hs <-dplyr::rename(hs,ind1=ind1.x, clust2=clust2.y, clust1=clust1.x) %>% dplyr::select(-ind1.y,-Ind.y, -Ind.x, -clust2.x, -clust1.y)
hs <-dplyr::select(hs, -Pop.x, -Ocean.x, -Pop.y, -Ocean.y)

#split by bay
bay17 <- hs %>% base::split(.$Bay1)
clust17 <-unite(clust17, ConYear, Con, Year, remove=FALSE) %>%mutate(clust=as.factor(clust))
parentcols <-c("1"="#d9d9d9","2"="#969696","3"="#525252","4"="#000000")

#MATTITUCK
Mt17 <-bay17$Mt
Mt17 <-mutate(Mt17, clust1=as.factor(clust1), clust2=as.factor(clust2))
Mt17$ConYear1 <-revalue(Mt17$ConYear1,c("1_2015"="2015 early","2_2015"="2015 late","1_2016"="2016 early","2_2016"="2016 late") )
Mt17$ConYear2 <-revalue(Mt17$ConYear2,c("1_2015"="2015 early","2_2015"="2015 late","1_2016"="2016 early","2_2016"="2016 late") )

Mtclust <-filter(clust17, Bay=="Mt")
Mtclust$ConYear <-revalue(Mtclust$ConYear,c("1_2015"="2015 early","2_2015"="2015 late","1_2016"="2016 early","2_2016"="2016 late") )


#scatter ind1 vs. ind2 - This one is good. 
ggplot(Mt17, aes(ind1, ind2))+
  geom_point(aes(color=clust1, size=prob))+ #options: shape = ConYear1 within aes, size=3 outside of aes
  #scale_color_brewer(palette ="Paired")+
  facet_grid(ConYear1~ConYear2, scales="free", drop=FALSE)+ #margins =TRUE is interesting. 
  xlab("offspring 1")+ylab("offspring 2")+
  theme(axis.text = element_text(size = 5),axis.title = element_text(size =10),axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = "white", colour = "black"),strip.text = element_text(colour = 'black', face="bold"))
ggsave('Mtbubbles.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 8, height = 8)


#barplot count of family within conyear. 
#clust 1 conyear 1  --- This one is good. 
#barplot count of family within conyear. 
Mtbbs <- ddply(Mtclust,~clust, summarize, countgroups= n_distinct(ConYear))
Mtbb17 <-left_join(Mtclust, Mtbbs, by="clust") %>%mutate(countgroups=as.factor(countgroups), ConYear=as.factor(ConYear))

Mtbb17 %>%
  ggplot(aes(clust, fill=countgroups))+# fill by number of appearances across groups
  scale_fill_manual(values=parentcols, guide=FALSE)+ #turn this off if you're doing family colors
  geom_histogram(stat="count")+
  facet_wrap(~ConYear, nrow=4)+
  xlab("family")+ylab("number of offspring")+
  theme(strip.background =element_rect(fill="white"),strip.text = element_text(colour = 'black', face="bold"),
        axis.text.y = element_text(size = 12),axis.text.x=element_text(size = 12),axis.title = element_text(size =12),
        panel.background = element_rect(fill = "white", colour = "black"))
ggsave('Mtfam17.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 4, height = 8)

#fathers
Mtdads <- ddply(Mtclust, ~father, summarize, countgroups= n_distinct(ConYear))
Mtdad17 <-Mtclust %>% left_join(Mtdads, by="father") %>%mutate(countgroups=as.factor(countgroups))
Mtdad17 %>%
ggplot(aes(father, fill=countgroups))+
  geom_histogram(stat="count")+
  scale_y_discrete(limits=c("1","2","3"))+
  scale_fill_manual(values=parentcols,guide=FALSE)+
  facet_wrap(~ConYear, nrow=4)+
  xlab("father")+ylab("number of offspring")+
  theme(strip.background =element_rect(fill="white"),strip.text = element_text(colour = 'black', face="bold"),
        axis.text.y = element_text(size = 12),axis.text.x=element_text(size = 12),axis.title = element_text(size =12),
        panel.background = element_rect(fill = "white", colour = "black"))
ggsave('Mtdads17.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 4, height = 8)


#mothers
Mtmoms <-Mtclust %>%  ddply(~mother, summarize, countgroups= n_distinct(ConYear))
Mtmom17 <-Mtclust %>% left_join(Mtmoms, by="mother") %>%mutate(countgroups=as.factor(countgroups))
Mtmom17 %>%
  ggplot(aes(mother, fill=countgroups))+
  geom_histogram(stat="count")+
  scale_y_discrete(limits=c("1","2", "3"))+
  scale_fill_manual(values=parentcols,guide=FALSE)+
  facet_wrap(~ConYear, nrow=4)+
  xlab("mother")+ylab("number of offspring")+
  theme(strip.background =element_rect(fill="white"),strip.text = element_text(colour = 'black', face="bold"),
        axis.text.y = element_text(size = 12),axis.text.x=element_text(size = 12),axis.title = element_text(size =12),
        panel.background = element_rect(fill = "white", colour = "black"))
ggsave('Mtmoms17.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 4, height = 8)

#compare families to relatedness
#join relatedness info to colony info
rel_mt <-mutate(rel_mt, ind1.id=as.factor(ind1.id), ind2.id=as.factor(ind2.id))
sibjoin <-dplyr::select(Mtclust, clust, Ind, father, mother, prob) %>% mutate(Ind=as.factor(Ind))
sj1 <-dplyr::rename(sibjoin, ind1.id=Ind, clust1=clust, father1=father, mother1=mother, prob1=prob) 
sj <-left_join(sj1, rel_mt, by=c("ind1.id"))
sj2 <-dplyr::rename(sibjoin, ind2.id=Ind, clust2=clust, father2=father, mother2=mother, prob2=prob)
sj <-left_join(sj2,sj, by=c("ind2.id"))

#find full siblings
siblings <-filter(sj, father1==father2 & mother1==mother2)








#SHINNECOCK BAY
Shin17 <-bay17$Shin
Shin17 <-mutate(Shin17, clust1=as.factor(clust1), clust2=as.factor(clust2))
Shin17$ConYear1 <-revalue(Shin17$ConYear1,c("1_2016"="2016 early","2_2016"="2016 late","1_2017"="2017 early","2_2017"="2017 late") )
Shin17$ConYear2 <-revalue(Shin17$ConYear2,c("1_2016"="2016 early","2_2016"="2016 late","1_2017"="2017 early","2_2017"="2017 late") )

Shiclust <-filter(clust17, Bay=="Shin")
Shiclust$ConYear <-revalue(Shiclust$ConYear,c("1_2016"="2016 early","2_2016"="2016 late","1_2017"="2017 early","2_2017"="2017 late") )


#scatter ind1 vs. ind2 - This one is good. 
ggplot(Shin17, aes(ind1, ind2))+
  geom_point(aes(color=clust1, size=prob))+ #options: shape = ConYear1 within aes, size=3 outside of aes
  #scale_color_brewer(palette ="Paired")+
  facet_grid(ConYear1~ConYear2, scales="free", drop=FALSE)+ #margins =TRUE is interesting. 
  xlab("offspring 1")+ylab("offspring 2")+
  theme(axis.text = element_text(size = 5),axis.title = element_text(size =10),axis.text.x = element_text(angle = 90),
        panel.background = element_rect(fill = "white", colour = "black"),strip.text = element_text(colour = 'black', face="bold"))
ggsave('Shibubbles.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 8, height = 8)

#barplot count of family within conyear. 
Shibbs <-Shiclust %>% ddply(~clust, summarize, countgroups= n_distinct(ConYear))
Shibb17 <-Shiclust %>% left_join(Shibbs, by="clust") %>%mutate(countgroups=as.factor(countgroups))
Shibb17 %>%
  ggplot(aes(clust,fill=countgroups))+# fill by number of appearances across groups
  scale_fill_manual(values=parentcols,guide=FALSE)+ #turn this off if you're doing family colors
  geom_histogram(stat="count")+
  facet_wrap(~ConYear, nrow=4)+
  xlab("family")+ylab("number of offspring")+
  theme(strip.background =element_rect(fill="white"),strip.text = element_text(colour = 'black', face="bold"),
        axis.text.y = element_text(size = 12),axis.text.x=element_text(size = 12),axis.title = element_text(size =12),panel.background = element_rect(fill = "white", colour = "black"))
ggsave('Shifam17.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 4, height = 8)

#fathers
Shidads <-Shiclust %>%ddply(~father, summarize, countgroups= n_distinct(ConYear))
Shidad17 <-Shiclust %>% left_join(Shidads, by="father") %>%mutate(countgroups=as.factor(countgroups))
Shidad17 %>%
  ggplot(aes(father, fill=countgroups))+
  geom_histogram(stat="count")+
  scale_y_discrete(limits=c("1","2","3"))+
  scale_fill_manual(values=parentcols,guide=FALSE)+
  facet_wrap(~ConYear, nrow=4)+
  xlab("father")+ylab("number of offspring")+
  theme(strip.background =element_rect(fill="white"),strip.text = element_text(colour = 'black', face="bold"),
        axis.text.y = element_text(size = 10),axis.text.x=element_text(size = 5),axis.title = element_text(size =10),
        panel.background = element_rect(fill = "white", colour = "black"))
ggsave('Shidads17.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 4, height = 8)

#mothers
Shimoms <-Shiclust %>% ddply(~mother, summarize, countgroups= n_distinct(ConYear))
Shimom17 <-Shiclust %>% left_join(Shimoms, by="mother") %>%mutate(countgroups=as.factor(countgroups))
Shimom17 %>%
  ggplot(aes(mother, fill=countgroups))+
  geom_histogram(stat="count")+
  scale_y_discrete(limits=c("1","2","3"))+
  scale_fill_manual(values=parentcols, guide=FALSE)+
  facet_wrap(~ConYear, nrow=4)+
  xlab("mother")+ylab("number of offspring")+
  theme(strip.background =element_rect(fill="white"),strip.text = element_text(colour = 'black', face="bold"),
        axis.text.y = element_text(size = 12),axis.text.x=element_text(size = 12),axis.title = element_text(size =10),
        panel.background = element_rect(fill = "white", colour = "black"))
ggsave('Shimoms17.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 4, height = 8)










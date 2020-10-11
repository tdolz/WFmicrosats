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
Mt17 <-bay17$Mt

#try some kind of heatmap
ggplot(Mt17, aes(ind1, ind2))+ #switch between AS and percentsib
  geom_tile(aes(fill=prob), show.legend = TRUE)+ #fill is continuos probability
  scale_fill_gradient(low ="#ffeda0" , high = "#800026")+ #warm colors
  #geom_text(color="white", size=5)+
  coord_fixed(ratio = 1) +
  #facet_grid(ConYear1~ConYear2)+
  #facet_grid(clust1~clust2)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size =12),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))

#average tile
sum_prob <-ddply(Mt17, ConYear1~ConYear2, summarize, avr = mean(prob))
sum_prob <-mutate(sum_prob, cy1 = pmin(ConYear1,ConYear2), cy2 =pmax(ConYear1,ConYear2)) %>% arrange(cy1)
ggplot(sum_prob, aes(cy1, cy2))+ 
  geom_tile(aes(fill=avr), show.legend = TRUE)+
  scale_fill_gradient(low ="#ffeda0" , high = "#800026")+ #warm colors
  #geom_text(color="white", size=5)+
  coord_fixed(ratio = 1) +
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size =12),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))

Mt17 <-mutate(Mt17, clust1=as.factor(clust1), clust2=as.factor(clust2))

#scatter ind1 vs. ind2 - This one is good. 
ggplot(Mt17, aes(ind1, ind2))+
  geom_point(aes(color=clust1, size=prob))+ #options: shape = ConYear1 within aes, size=3 outside of aes
  facet_grid(ConYear1~ConYear2, scales="free", drop=FALSE)+ #margins =TRUE is interesting. 
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 5),axis.title = element_text(size =5),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))

#scatter ind1 vs. ind1
ggplot(Mt17, aes(ind1, ind1))+
  geom_point(aes(color=clust1), size=3)+ #you can add shape = ConYear1 within aes
  facet_wrap(~ConYear1, scales="free", drop=FALSE)+ #margins =TRUE is interesting. 
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 5),axis.title = element_text(size =5),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
#this tells you that families are spread across clusters for the most part. 

#barplot count of family within conyear. 
#clust 1 conyear 1  --- This one is good. 
ggplot(Mt17, aes(clust1))+
geom_histogram(stat="count")+
  facet_grid(~ConYear1)
#switch
ggplot(Mt17, aes(clust2))+
  geom_histogram(stat="count")+
  facet_grid(~ConYear2)

#try this
exMt <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/experimentalMTmatrix.csv", header = FALSE) #csv version 
exMT <-data.matrix(exMt)
exMT[is.na(exMT)] <-0

heatmap(exMT, Rowv=Mt17$ind1)

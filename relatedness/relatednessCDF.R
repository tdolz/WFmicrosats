#### RELATEDNESS CDF investigation #########
### November 2, 2020 #####

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

#convert to the right data format
setPop(wfpopLD) <-~Bay/Year
df <-genind2df(wfpopLD, usepop = FALSE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(35,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34)]
df[df=="NA"] <- 0 # missing data must be 0
write.table(df, "scratch",row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) #write it.
genotypedata <- readgenotypedata("scratch")# import input file as list (gdata, nloci, nalleles, ninds, freqs)

#if you don't want to do them all again, try just doing one... 
simdata <-familysim(genotypedata$freqs,100)
outputTRI <- coancestry(simdata, trioml=2)
simreltri <- cleanuprvals(outputTRI$relatedness , 100)

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
parentoffspring <-relply2[2,2] #not interesting when only doing yoy
unrelated <-relply2[4,2]
halfsibs.sem <-relply2[1,3]

#trioml estimate (Wang) 
relatedness_triad <- coancestry(genotypedata$gdata,trioml = 2) #remember to turn on 1 or 2, for no CI or w/ CI. 
relatedT <- relatedness_triad$relatedness %>%
  dplyr::select(pair.no, ind1.id, ind2.id, trioml)
relatedCI <- relatedness_triad$relatedness.ci95 %>%
  dplyr::select(pair.no, ind1.id, ind2.id, trioml.low, trioml.high)
relatedT <-left_join(relatedT, relatedCI)

#attach strata information to the relatedness estimation. 
h <-melt(relatedT, id=c("pair.no","ind1.id","ind2.id"))
h <-dplyr::rename(h, Estimator=variable,Relatedness_Value=value)
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

#find related fish from the same bay!
same.bay <-filter(hs, Bay1==Bay2)
diff.bay <-filter(hs, Bay1!=Bay2)

#cohorts/contingents in MT
Mtpairs <-filter(same.bay, Bay1=="Mt")
Shipairs <-filter(same.bay, Bay1=="Shin")

######## CDFs##########

#Simulation CDF
HScd <-filter(simreltri, group=="HSHS")%>% arrange(trioml) #half sibs
HSCDF = cumsum(HScd$trioml)
SUMhs<-sum(HScd$trioml)
HScd <- mutate(HScd, SurveyCDFFinal=HSCDF/SUMhs) 
#full sibs CDF
SBcd <-filter(simreltri, group=="SBSB")%>% arrange(trioml) #full sibs
SBCDF = cumsum(SBcd$trioml)
SUMsb<-sum(SBcd$trioml)
SBcd <- mutate(SBcd, SurveyCDFFinal=SBCDF/SUMsb) 
##Parent-offspring CDF
POcd <-filter(simreltri, group=="POPO")%>% arrange(trioml) #parent-offspring
POCDF = cumsum(POcd$trioml)
SUMpo<-sum(POcd$trioml)
POcd <- mutate(POcd, SurveyCDFFinal=POCDF/SUMpo) 
##Unrelated CDF
URcd <-filter(simreltri, group=="URUR")%>% arrange(trioml) #unrelated
URCDF = cumsum(URcd$trioml)
SUMur<-sum(URcd$trioml)
URcd <- mutate(URcd, SurveyCDFFinal=URCDF/SUMur) 

#### Real Data CDFs ####
#bay CDFS
#Jamaica
JAMcd <-filter(same.bay, Bay1=="Jam")%>% arrange(Relatedness_Value) #Jamaica
JAMCDF = cumsum(JAMcd$Relatedness_Value)
SUMJAM<-sum(JAMcd$Relatedness_Value)
JAMcd <- mutate(JAMcd, SurveyCDFFinal=JAMCDF/SUMJAM) 
#Moriches
MORcd <-filter(same.bay, Bay1=="Mor")%>% arrange(Relatedness_Value) #Moriches
MORCDF = cumsum(MORcd$Relatedness_Value)
SUMMOR<-sum(MORcd$Relatedness_Value)
MORcd <- mutate(MORcd, SurveyCDFFinal=MORCDF/SUMMOR)
#Mattituck
MTcd <-filter(same.bay, Bay1=="Mt")%>% arrange(Relatedness_Value) #Mattituck
MTCDF = cumsum(MTcd$Relatedness_Value)
SUMMT<-sum(MTcd$Relatedness_Value)
MTcd <- mutate(MTcd, SurveyCDFFinal=MTCDF/SUMMT) 
#Napeague
NAPcd <-filter(same.bay, Bay1=="Nap")%>% arrange(Relatedness_Value) #NAPaica
NAPCDF = cumsum(NAPcd$Relatedness_Value)
SUMNAP<-sum(NAPcd$Relatedness_Value)
NAPcd <- mutate(NAPcd, SurveyCDFFinal=NAPCDF/SUMNAP) 
#Shinnecock
SHIcd <-filter(same.bay, Bay1=="Shin")%>% arrange(Relatedness_Value) #SHIaica
SHICDF = cumsum(SHIcd$Relatedness_Value)
SUMSHI<-sum(SHIcd$Relatedness_Value)
SHIcd <- mutate(SHIcd, SurveyCDFFinal=SHICDF/SUMSHI) 

drabcolors <-c("#d0d1e6","#a6bddb", "#67a9cf", "#1c9099", "#016450")
#plot
a <-  ggplot(HScd, aes(x=trioml, y=SurveyCDFFinal))
a+geom_line(color="#666666",linetype="dotted")+ #halfsibs
  geom_line(aes(x=trioml, y=SurveyCDFFinal), color="#000000", linetype="dotted", data=SBcd)+ #full sibs
  geom_line(aes(x=trioml, y=SurveyCDFFinal), color="#333333", linetype="dotted",data=POcd)+ #parent offspring
  geom_line(aes(x=trioml, y=SurveyCDFFinal), color="#999999", linetype="dotted",data=URcd)+ #unrelated
  geom_line(aes(x=Relatedness_Value, y=SurveyCDFFinal), color="#d0d1e6", data=JAMcd)+ #jamaica
  geom_line(aes(x=Relatedness_Value, y=SurveyCDFFinal), color="#a6bddb", data=MORcd)+ #Moriches
  geom_line(aes(x=Relatedness_Value, y=SurveyCDFFinal), color="#67a9cf", data=MTcd)+ #Mattituck
  geom_line(aes(x=Relatedness_Value, y=SurveyCDFFinal), color="#1c9099", data=NAPcd)+ #Napeague
  geom_line(aes(x=Relatedness_Value, y=SurveyCDFFinal), color="#016450", data=SHIcd)+ #Shinnecock
  ylab("CDF")+xlab("Relatedness")+
theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
      panel.grid.major = element_line(colour = "white"))+guides(fill = FALSE, colour = FALSE) 

#We could repeat this graph for cohort pairs.
#another way we could do this is to look at different bay pairs? 
# but instead of doing that lets look up the quantile of every pair in each CDF. 

hs.quant <-(HScd[HScd$SurveyCDFFinal=='0.5',])

hs.quant <-c()
for (i in 1:length(relatedT$Relatedness_Value)){
  hs.quant[i] <-quantile(HScd, )
}

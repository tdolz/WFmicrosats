#format data for microdrop and NE estimator. 

###September 19, 2020

#install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)
library('plyr')
library("dplyr")
library('readr')
library('ggplot2')

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")
wfpop4df <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_20_sept2020_doubl0.csv", header = TRUE) #csv version 

allele1 <- dplyr::select(wfpop4df, Ind,Pop, "J42.1", "WF22.1","PAM27.1","PAM79.1","WF27.1","WF16.1", "WF33.1","WF3.1","WF517.1","WF196.1","WF223.1","WF421.1",
                      "A441.1","PAM21.1", "Psy087.1", "Psy022.1", "WF06.1","WF12.1","WF01.1","WF32.1")

names(allele1) <-c("Ind","Pop","J42", "WF22","PAM27","PAM79","WF27","WF16", "WF33","WF3","WF517","WF196","WF223","WF421",
                   "A441","PAM21", "Psy087", "Psy022", "WF06","WF12","WF01","WF32")


allele2 <-dplyr::select(wfpop4df, Ind,Pop, "J42.2", "WF22.2","PAM27.2","PAM79.2","WF27.2","WF16.2", "WF33.2","WF3.2","WF517.2","WF196.2","WF223.2","WF421.2",
                         "A441.2","PAM21.2", "Psy087.2", "Psy022.2", "WF06.2","WF12.2","WF01.2","WF32.2")

names(allele2) <-c("Ind","Pop","J42", "WF22","PAM27","PAM79","WF27","WF16", "WF33","WF3","WF517","WF196","WF223","WF421",
                   "A441","PAM21", "Psy087", "Psy022", "WF06","WF12","WF01","WF32")

microdrop <-full_join(allele1,allele2)
microdrop <-arrange(microdrop, Ind)
microdrop[microdrop=="0"] <- -9 # missing data must be -9
microdrop[microdrop=="NA"]<- -9

#check for places where one allele is filled in an the other isn't
check <-filter(microdrop, WF27=="-9")%>%count(Ind); filter(check, n < 2)
#WF27 was the only problematic one before, and we eliminated that loci in our previous analysis. 

write_delim(microdrop, "microdrop_datafile")

### format for NE estimator
nepop <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_17_sept2020_doubl0ABC.csv", header = TRUE) #csv version 

nepop <-tidyr::separate(nepop, Pop, c("Oc", "Bay", "Con", "Year"))
nepop16 <-filter(nepop, Year=="2016") %>% dplyr::select(-Oc, -Con, -Year) %>% dplyr::rename(Pop=Bay)
nepop16 <- dplyr::select(nepop16, -WF01.1, -WF01.2, -WF06.1, -WF06.2, -WF27.1, -WF27.2) #exclude loci with missing data and LD. 
nepop16 <- tidyr::unite(nepop16, J42, J42.1, J42.2, sep="") %>% 
  tidyr::unite(WF22, WF22.1, WF22.2, sep="") %>% tidyr::unite(PAM27, PAM27.1, PAM27.2, sep="")%>%
  tidyr::unite(PAM79, PAM79.1, PAM79.2, sep="") %>% tidyr::unite(WF16, WF16.1, WF16.2, sep="")%>%
  tidyr::unite(WF33, WF33.1, WF33.2, sep="") %>% tidyr::unite(WF3, WF3.1, WF3.2, sep="")%>%
  tidyr::unite(WF517, WF517.1, WF517.2, sep="") %>% tidyr::unite(WF196, WF196.1, WF196.2, sep="")%>%
  tidyr::unite(WF223, WF223.1, WF223.2, sep="") %>% tidyr::unite(WF421, WF421.1, WF421.2, sep="")%>%
  tidyr::unite(A441, A441.1, A441.2, sep="") %>% tidyr::unite(PAM21, PAM21.1, PAM21.2, sep="")%>%
  tidyr::unite(Psy087, Psy087.1, Psy087.2, sep="") %>% tidyr::unite(Psy022, Psy022.1, Psy022.2, sep="")%>%
  tidyr::unite(WF12, WF12.1, WF12.2, sep="") %>% tidyr::unite(WF32, WF32.1, WF32.2, sep="")
  
nepop16[nepop16=="00"] <- "000000" # missing data must be -9
nepop16 <-dplyr::select(nepop16,-Pop)

write_delim(nepop16, "nepop_datafile")

#nepop16 <-filter(nepop, Pop %in% c("Atl_Jam_6_2016","Atl_Mor_6_2016","Atl_Mt_1_2016","Atl_Mt_2_2016","Atl_Nap_6_2016","Atl_Shin_1_2016","Atl_Shin_2_2016"))


### Format data for COLONY
### format for NE estimator
nepop <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_17_sept2020_doubl0ABC.csv", header = TRUE) #csv version 
mtsex <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/fishsex.csv", header = TRUE)
mtadults <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/mtinfo.csv", header = TRUE)

mtadults <-dplyr::rename(mtadults, sex2=sex, Ind=samp) %>% dplyr::select(-sex2)
mtadults <-left_join(mtadults, mtsex, by="Ind") %>%mutate(Ind=as.factor(Ind))

nepop<-separate(nepop, Pop, c("Oc","Bay","Con","Year"), remove=FALSE)%>% dplyr::select(-Oc)
mtpop <- filter(nepop, Bay=="Mt") %>% dplyr::select(-Bay,-Year) 
mtpop <-left_join(mtpop, mtadults, by=c("Ind"))
#write.csv(mtpop, file="/Users//tdolan/Documents//WIP research//microsats//colony2//mtpop.csv")

cfs <- filter(mtpop, sex=="F")
cms <- filter(mtpop, sex=="M")
ofs <-filter(mtpop, Con %in% c("1","2"))

#all the hatch and spawning dates
mtage <-read.csv(file="mtcohortreassignAug182020.csv", header=TRUE)
mtage <-mutate(mtage, bday=as.Date(hatch.date), cohort=as.factor(new.cohort)) %>%slice(1:49) %>% dplyr::select(-old.pop)
spawndate <-read.csv(file="mtspawndate.csv", header=TRUE)
#merge spawn date and fish information. 
spawndate <-dplyr::rename(spawndate, fishID=Tag.ID)
mtage2 <-left_join(mtage,spawndate, by="fishID")

#join the spawning info to mtadults
mtage2 <- dplyr::select(mtage2, -cohort)%>% dplyr::rename(Ind=samp, FishID=fishID, location=Location) %>% mutate(Ind=as.factor(Ind))
mtadults <-dplyr::select(mtadults, -Fin.Clip, -stage, -put.age, -year.class)
mtinfo <-full_join(mtage2, mtadults, by=c("Ind","Year","FishID","location","TL","catch.date"))


### plot NE comparison
NEests <-read.csv(file="compareNE.csv", header=TRUE)

NEests %>%
  ggplot(aes(x=bay,y=est, group=run, color=run))+
  scale_color_viridis(discrete=TRUE)+
  geom_pointrange(aes(ymin=low,ymax=high))+
  xlab("")+ylab("estimated Nb")+
  #ylim(0,500)+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))



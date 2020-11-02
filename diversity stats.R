#### Diversity Stats ######
# made on October 28, 2020 from Microsats 2020
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
library('readr')
library("reshape2")
library("strataG")
library("viridis")
library("purrr")

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")

#toolpack
sem <-function(x) sd(x)/sqrt(length(x))



##### Formating the dataset #####
wfpop <- read.genalex("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_17_sept20204genalex_doubl0ABC.csv")
wfpop4df <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_17_sept2020_doubl0ABC.csv", header = TRUE) #csv

splitStrata(wfpop) <-~Ocean/Bay/Con/Year
setPop(wfpop) <-~Bay

#clean dataset. 
wfpopLD <-genclone2genind(wfpop)
all_loci <- locNames(wfpopLD)# create vector of all loci
removeloc <- c("WF06","WF27")# create vector containing loci to remove
keeploc <- setdiff(all_loci, removeloc)# create list of loci to keep
wfpopLD <- wfpopLD[loc = keeploc]
wfpopLD <-wfpopLD %>% missingno("loci", cutoff=0.30) # removes loci where overall missingness > 30%, so WF01
wfpopLD <- wfpopLD %>% missingno("geno", cutoff=0.25) # remove samples that are > 25% missing

## Subset databases
#shinnecock only database
setPop(wfpopLD) <-~Bay
shinco <-popsub(wfpopLD, sublist=c("Shin"))
setPop(shinco) <-~Bay/Con/Year

#mattituck only database, exclude MT_5
setPop(wfpopLD) <-~Bay/Con
mtco <-popsub(wfpopLD, sublist=c("Mt_1","Mt_2","Mt_3","Mt_4"))
setPop(mtco) <-~Bay/Con/Year

####HWE - New tests with Benjamini-Hochberg.####
#bh adjust function
pradj <-function(df){
  df <-as.data.frame(df) %>%tibble::rownames_to_column() %>%
    mutate(Pr.adj=p.adjust(Pr.exact, method="BH"))%>%mutate(sig=ifelse(Pr.adj < 0.05, "sig","not sig"))}

setPop(wfpopLD) <-~Ocean
wfhwe.pop <- seppop(wfpopLD) %>% lapply(hw.test)
Atlh <-pradj(wfhwe.pop$Atl)

setPop(wfpopLD) <-~Bay
wfhwe.pop <- seppop(wfpopLD) %>% lapply(hw.test) 
#%>%pmap(pradj)
# i can't figure out how to do this at all. it would be easier to just bonferroni correct everything... 

#using the bejamnini
jamh <-pradj(wfhwe.pop$Jam)
morh <-pradj(wfhwe.pop$Mor)
Mth <-pradj(wfhwe.pop$Mt)
Naph <-pradj(wfhwe.pop$Nap)
Shih <-pradj(wfhwe.pop$Shin)
#we could make a new levelplot figure for this, but hold off for now.  Code for that is in Microsats_20.R

#Mattituck HWE
setPop(mtco) <-~Bay/Con/Year
mpop <-seppop(mtco) %>% lapply(hw.test)
(pradj(mpop$Mt_1_2015))
(pradj(mpop$Mt_2_2015))
(pradj(mpop$Mt_1_2016))
(pradj(mpop$Mt_2_2016))
(pradj(mpop$Mt_3_adults))
(pradj(mpop$Mt_4_adults))

#Shinnecock HWE
setPop(shinco) <-~Bay/Con/Year
spop <-seppop(shinco) %>% lapply(hw.test)
(pradj(spop$Shin_1_2016))
(pradj(spop$Shin_2_2016))
(pradj(spop$Shin_1_2017))
(pradj(spop$Shin_2_2017))


## Diversity stats ### 
# Ho is observed heterozygosity
# Hs is expected heterozygosity
# Fis inbreeding coeffcient
# Ht? Nei's gene diversity

library(tidyverse)
library(ggplot2)
library(adegenet)
library(hierfstat)

# genind object of genotypes
wfpopLD

#dataframe with sample information
setPop(wfpopLD) <- ~Bay/Con/Year
df <-genind2df(wfpopLD, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(36,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)] #17 loci
df[df=="NA"] <- 0 # missing data must be 0
SampleInfo <- dplyr::select(df, Ind, pop)
SampleInfo <- separate(SampleInfo, pop, c("Bay","Con","Year"))
SampleInfo <-mutate(SampleInfo, Ocean="Atl")
SampleInfo <-SampleInfo[,c(1,5,2,3,4)]

setPop(wfpopLD)<-~Ocean
gen_oce <-seppop(wfpopLD)

setPop(wfpopLD)<-~Bay
gen_bay <-seppop(wfpopLD)

setPop(shinco) <-~Bay/Con/Year
gen_shin <-seppop(shinco)

setPop(mtco) <-~Bay/Con/Year
gen_mt <-seppop(mtco)


gen_grp <-c(gen_oce,gen_bay,gen_shin,gen_mt)
gen_grp[["ALL"]] <-wfpopLD

# calculate allelic richness/diversity stats ====
loc_stats <- list()
for (p in names(gen_grp)) {
  locA <- locus_table(gen_grp[[p]], index = "shannon") %>%
    as.data.frame() %>%
    rownames_to_column("LOCUS")
  locB <- locus_table(gen_grp[[p]], index = "simpson") %>%
    as.data.frame() %>%
    rownames_to_column("LOCUS")
  temp <- left_join(locA, locB)
  locC <- locus_table(gen_grp[[p]], index = "invsimpson") %>%
    as.data.frame() %>%
    rownames_to_column("LOCUS")
  loc_stats[[p]] <- left_join(temp, locC)
}
#Now loc.stats is a nice little list of all these cool metrics for all the groups! 

#here it is in one dataframe
loc_stats <- plyr::ldply(loc_stats, data.frame) %>%
  dplyr::select(-Hexp) %>%
  dplyr::rename(GRP = `.id`,
                SIMPSON_IDX = `X1.D`,
                N_ALLELES = allele,
                SHANNON_IDX = H,
                STODD_TAYLOR_IDX = G,
                EVENNESS = Evenness)

# calculate genetic diversity stats (heterozygosity-based) ====
loc_stats_2 <- list()

for (p in names(gen_grp)) {
  dat <- hierfstat:::.genind2hierfstat(gen_grp[[p]])
  stats <- basic.stats(dat)
  loc_stats_2[[p]] <- stats$perloc %>%
    rownames_to_column("LOCUS")
}
#Ht and Hs are already the same at this point... 


# combine into single data frame ====
loc_stats_2 <- plyr::ldply(loc_stats_2, data.frame) %>%
  dplyr::rename(GRP = `.id`)

loc_stats <- left_join(loc_stats, loc_stats_2) %>%
  dplyr::select(GRP, LOCUS, N_ALLELES, EVENNESS, Ho, Hs, Ht, Fis, SHANNON_IDX, SIMPSON_IDX, STODD_TAYLOR_IDX) %>%
  filter(LOCUS != "mean")

loc_stats <-mutate(loc_stats, Hdef=Ho-Hs)

loc_stats[is.na(loc_stats)] <- NA
#write.csv(loc_stats,file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/loc_stats17.csv")
pivlocstats <-pivot_longer(loc_stats,cols = c("N_ALLELES", "SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis","SHANNON_IDX","STODD_TAYLOR_IDX", "Hdef"),
                           names_to="variable", values_to="value")
sumlocstats <-ddply(pivlocstats, GRP~variable, summarize, meanvar=mean(value), sdvar=sd(value))
#write.csv(sumlocstats,file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/sumlocstats17.csv")

###rareified allelic richness###
# overall ====
setPop(wfpopLD) <- ~Ocean
dat <- hierfstat:::.genind2hierfstat(wfpopLD)
ar <- allelic.richness(dat,diploid = TRUE)

ar <- as.data.frame(ar$Ar) %>%
  rownames_to_column("LOCUS") %>%
  dplyr::rename(ALL = V1)

# By Bay ===
# by region ====
setPop(wfpopLD) <- ~Bay
dat <- hierfstat:::.genind2hierfstat(wfpopLD, pop = wfpopLD@pop)
df <- allelic.richness(dat,diploid = TRUE)
df <- as.data.frame(df$Ar) %>%
  rownames_to_column("LOCUS") %>%
  dplyr::rename(Nap = V1,Mor = V2,Jam = V3,Shin = V4,Mt = V5)
ar <- left_join(ar, df)
#for rarified alleles, you have to make a new Allelic Richness analysis for each group you are comparing. 
meltar2 <- pivot_longer(ar, cols=c("Mt","Shin","Nap","Mor","Jam","ALL"),names_to="variable", values_to="value")
meltar3 <-filter(meltar2, variable != "ALL") #it doenst like it when you exclude the "all" category... 
meltar_all <-filter(meltar2, variable == "ALL")
meltar4 <-ddply(meltar2, ~variable, summarize,meanar=mean(value), sdar=sd(value)) #excel version

## AR mattituck
setPop(mtco) <- ~Bay/Con/Year
dat <- hierfstat:::.genind2hierfstat(mtco, pop = mtco@pop)
df <- allelic.richness(dat,diploid = TRUE)
df <- as.data.frame(df$Ar) %>%
  rownames_to_column("LOCUS") %>%
  dplyr::rename(Mt_1_2015=V1, Mt_1_2016=V2, Mt_2_2015=V3, Mt_2_2016=V4, Mt_3_adults=V5, Mt_4_adults=V6)
meltar2 <- pivot_longer(df, cols=c("Mt_1_2015", "Mt_1_2016", "Mt_2_2015", "Mt_2_2016", "Mt_3_adults", "Mt_4_adults"),names_to="variable", values_to="value")
mtar <-filter(meltar2, variable != "ALL") #it doenst like it when you exclude the "all" category... 
mtdd <-ddply(mtar, ~variable, summarize,meanar=mean(value), sdar=sd(value)) #excel version

## AR Shinnecock
setPop(shinco) <- ~Bay/Con/Year
dat <- hierfstat:::.genind2hierfstat(shinco, pop = shinco@pop)
df <- allelic.richness(dat,diploid = TRUE)
df <- as.data.frame(df$Ar) %>%
  rownames_to_column("LOCUS") %>%
  dplyr::rename(Shin_1_2016=V1, Shin_1_2017=V2, Shin_2_2016=V3, Shin_2_2017=V4)
meltar2 <- pivot_longer(df, cols=c("Shin_1_2016", "Shin_1_2017", "Shin_2_2016", "Shin_2_2017"),names_to="variable", values_to="value")
shar <-filter(meltar2, variable != "ALL") #it doenst like it when you exclude the "all" category... 
shdd <-ddply(shar, ~variable, summarize,meanar=mean(value), sdar=sd(value)) #excel version

#csv of allellic richness
AR <- bind_rows(meltar4,mtdd,shdd)
#write_excel_csv(AR,path="/Users/tdolan/Documents/WIP research/microsats/microsats_results/allelicrichness.csv")

####Private alleles ####
#bay
df <-genind2df(wfpopLD, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(36,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)]
df[df=="NA"] <- 0 # missing data must be 0
wf.g <-df2gtypes(df,ploidy=2)
library("data.table")
pA <-as.data.frame(privateAlleles(wf.g)) 
setDT(pA,keep.rownames=TRUE)
colnames(pA)[1] <- "LOCUS"
meltpA <- pivot_longer(pA, cols=c("Mt","Shin","Nap","Mor","Jam"),names_to="variable", values_to="value")

#mattituck
df <-genind2df(mtco, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(36,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)]
df[df=="NA"] <- 0 # missing data must be 0
wf.g <-df2gtypes(df,ploidy=2)
pA <-as.data.frame(privateAlleles(wf.g)) 
setDT(pA,keep.rownames=TRUE)
colnames(pA)[1] <- "LOCUS"
meltMTpa <- pivot_longer(pA, cols=c("Mt_1_2015", "Mt_1_2016", "Mt_2_2015", "Mt_2_2016", "Mt_3_adults", "Mt_4_adults"),names_to="variable", values_to="value")

#Shinnecock
df <-genind2df(shinco, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(36,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)]
df[df=="NA"] <- 0 # missing data must be 0
wf.g <-df2gtypes(df,ploidy=2)
pA <-as.data.frame(privateAlleles(wf.g)) 
setDT(pA,keep.rownames=TRUE)
colnames(pA)[1] <- "LOCUS"
meltSHpa <- pivot_longer(pA, cols=c("Shin_1_2016", "Shin_1_2017", "Shin_2_2016", "Shin_2_2017"),names_to="variable", values_to="value")


#csv of private alleles
PA <- bind_rows(meltpA,meltMTpa,meltSHpa)
PA <-as.data.frame(PA)%>% mutate(variable=as.factor(variable))
sumPA <-ddply(PA,~variable, summarize, meanPA=mean(value), sdPA=sd(value), sumPA=sum(value))
#write_excel_csv(sumPA,path="/Users/tdolan/Documents/WIP research/microsats/microsats_results/privatealleles.csv")

#heatmap of private alleles
#Giant barplot
PA2 <-filter(PA,variable %in% c("Nap","Mor","Jam","Shin","Shin_1_2016","Shin_2_2016","Shin_1_2017","Shin_2_2017","Mt","Mt_3_adults","Mt_4_adults","Mt_1_2015","Mt_2_2015","Mt_1_2016","Mt_2_2016"))
PA2$variable <-fct_relevel(PA2$variable, c("Jam","Mor","Mt","Nap","Shin","Mt_1_2015","Mt_2_2015","Mt_1_2016","Mt_2_2016", "Mt_3_adults","Mt_4_adults","Shin_1_2016","Shin_2_2016","Shin_1_2017","Shin_2_2017"))

PA2 %>%
  #filter(variable %in% c("Mt","Shin","Nap","Mor","Jam"))%>%
ggplot(aes(x = LOCUS, y = fct_rev(variable), fill=value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = ifelse(value > 0, value, "")), color="grey") +
  scale_fill_viridis(option="plasma") +
  coord_fixed(ratio = 1) +
  ylab("")+
  #theme_standard() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave('pabayheat17.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

#also calculate within bays. 
#bay
setPop(wfpopLD)<-~Bay/Con/Year
df <-genind2df(wfpopLD, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(36,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)]
df[df=="NA"] <- 0 # missing data must be 0
wf.g <-df2gtypes(df,ploidy=2)
library("data.table")
pA <-as.data.frame(privateAlleles(wf.g)) 
setDT(pA,keep.rownames=TRUE)
colnames(pA)[1] <- "LOCUS"
meltpA22 <- pivot_longer(pA, cols=c("Jam_6_2016", "Mor_6_2016", "Mt_1_2015", "Mt_1_2016", "Mt_2_2015", "Mt_2_2016", "Mt_3_adults", "Mt_4_adults", "Mt_5_2015", "Mt_5_2016", "Nap_6_2016", "Shin_1_2016", "Shin_1_2017", "Shin_2_2016", "Shin_2_2017"),names_to="variable", values_to="value")

meltpA22$variable <-fct_relevel(meltpA22$variable, c("Jam_6_2016", "Mor_6_2016","Nap_6_2016", "Mt_1_2015", "Mt_1_2016", "Mt_2_2015", "Mt_2_2016", "Mt_3_adults", "Mt_4_adults", "Mt_5_2015", "Mt_5_2016",  "Shin_1_2016", "Shin_1_2017", "Shin_2_2016", "Shin_2_2017"))

meltpA22 %>%
  filter(variable %in% c("Jam_6_2016", "Mor_6_2016", "Mt_1_2015", "Mt_1_2016", "Mt_2_2015", "Mt_2_2016", "Mt_3_adults", "Mt_4_adults", "Nap_6_2016", "Shin_1_2016", "Shin_1_2017", "Shin_2_2016", "Shin_2_2017"))%>%
  ggplot(aes(x = LOCUS, y = fct_rev(variable), fill=value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = ifelse(value > 0, value, "")), color="grey") +
  scale_fill_viridis(option="plasma") +
  coord_fixed(ratio = 1) +
  ylab("")+
  #theme_standard() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave('paheat17byc.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

#Combined heatmap. 
meltpA23 <- filter(meltpA22, variable %in% c("Mt_1_2016", "Mt_2_2015", "Mt_2_2016", "Mt_3_adults", "Mt_4_adults",  "Shin_1_2016", "Shin_1_2017", "Shin_2_2016", "Shin_2_2017")) %>% dplyr::rename(valueBYC=value)
PA3 <-filter(PA2,variable %in% c("Nap","Mor","Jam","Shin","Shin_1_2016","Shin_2_2016","Shin_1_2017","Shin_2_2017","Mt","Mt_3_adults","Mt_4_adults","Mt_1_2015","Mt_2_2015","Mt_1_2016","Mt_2_2016")) 
PA4=full_join(meltpA23,PA3)
PA4 <-mutate(PA4, valueBYC=ifelse(is.na(valueBYC),value,valueBYC))
PA4$variable <-fct_relevel(PA4$variable, c("Jam","Mor","Mt","Nap","Shin","Mt_1_2015","Mt_2_2015","Mt_1_2016","Mt_2_2016", "Mt_3_adults","Mt_4_adults","Shin_1_2016","Shin_2_2016","Shin_1_2017","Shin_2_2017"))


PA4 %>%
  ggplot(aes(x = LOCUS, y = fct_rev(variable), fill=valueBYC)) +
  geom_tile(color = "black") +
  geom_text(aes(label = ifelse(value > 0, value, "")), color="grey") +
  scale_fill_viridis(option="plasma") +
  coord_fixed(ratio = 1) +
  ylab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text = element_text(size = 12),axis.title = element_text(size = 12),panel.background = element_rect(fill = 'white', colour = 'black'),
      panel.grid.major = element_line(colour = "white"))
#



#### Recalculate HS vecause Hs and Ht are the same.###  
#I think an HDef column would be good where we calculate Ho-Hs for all of them.... but I am not sure Hs is correct because it's the same as Ht on 
#so we have this other method to calculate hs, but it's complicated so i gave up. 
g2h_all <-genind2hierfstat(wfpopLD)
g2hall <-basic.stats(g2h_all)
global <-as.data.frame(g2hall$perloc) %>% tibble::rownames_to_column("locus") %>% mutate(Bay="ALL")
Hsglobal <-dplyr::select(global, locus, Hs, Bay)

setPop(wfpopLD) <-~Bay
pops <-seppop(wfpopLD) %>% lapply(genind2hierfstat)
Hsj <-basic.stats(pops$Jam)
Hsj <-as.data.frame(Hsj$perloc)%>% tibble::rownames_to_column("locus") %>%dplyr::select(locus, Hs)%>%mutate(GRP="Jam")


#### Inbreeding - Internal Relatedness ####
library("Rhh")
setPop(wfpopLD) <-~Bay
wfp <-genind2df(wfpopLD, usepop = TRUE, oneColPerAll = TRUE) %>% tibble::rownames_to_column("Ind")
wfp <-mutate_at(wfp,vars(-pop, -Ind),as.numeric)
wfp[is.na(wfp)] <- 0
rel <- dplyr::select(wfp, -1, -2)
rel <-ir(rel)
rel2 <-as.data.frame(rel)
names(rel2) <-c("IR")
rel2 <-cbind(wfp,rel2) 
rel2 <-dplyr::select(rel2,Ind,pop,IR)
ddply(rel2, ~pop, summarize, meanir=mean(IR),sdir=sd(IR))
#all
setPop(wfpopLD) <-~Ocean
wfp <-genind2df(wfpopLD, usepop = TRUE, oneColPerAll = TRUE) %>% tibble::rownames_to_column("Ind")
wfp <-mutate_at(wfp,vars(-pop, -Ind),as.numeric)
wfp[is.na(wfp)] <- 0
rel <- dplyr::select(wfp, -1, -2)
rel <-ir(rel)
orel2 <-as.data.frame(rel)
names(orel2) <-c("IR")
orel2 <-cbind(wfp,orel2) 
orel2 <-dplyr::select(orel2,Ind,pop,IR)
ddply(orel2, ~pop, summarize, meanir=mean(IR),sdir=sd(IR))

#mattituck IR 
setPop(mtco) <-~Bay/Con/Year
wfp <-genind2df(mtco, usepop = TRUE, oneColPerAll = TRUE) %>% tibble::rownames_to_column("Ind")
wfp <-mutate_at(wfp,vars(-pop, -Ind),as.numeric)
wfp[is.na(wfp)] <- 0
mrel <- dplyr::select(wfp, -1, -2)
mrel <-ir(mrel)
mrel2 <-as.data.frame(mrel)
names(mrel2) <-c("IR")
mrel2 <-cbind(wfp,mrel2) 
mrel2 <-dplyr::select(mrel2,Ind,pop,IR)
ddply(mrel2, ~pop, summarize, meanir=mean(IR),sdir=sd(IR))
#Shinnecock IR 
setPop(shinco) <-~Bay/Con/Year
wfp <-genind2df(shinco, usepop = TRUE, oneColPerAll = TRUE) %>% tibble::rownames_to_column("Ind")
wfp <-mutate_at(wfp,vars(-pop, -Ind),as.numeric)
wfp[is.na(wfp)] <- 0
srel <- dplyr::select(wfp, -1, -2)
srel <-ir(srel)
srel2 <-as.data.frame(srel)
names(srel2) <-c("IR")
srel2 <-cbind(wfp,srel2) 
srel2 <-dplyr::select(srel2,Ind,pop,IR)
ddply(srel2, ~pop, summarize, meanir=mean(IR),sdir=sd(IR))

rel3 <-bind_rows(rel2,orel2,mrel2,srel2) 


#anova on relatedness: BAY
Irov <- lm(IR~pop, data=rel2,  na.action=na.omit)
car::Anova(Irov)
df<-df.residual(Irov)
MSerror<-deviance(Irov)/df
comparison <- HSD.test(Irov,c("pop"),MSerror=MSerror,alpha=0.05/5,  group=TRUE)
comparison

#anova on relatedness: mattituck
Irov <- lm(IR~pop, data=mrel2,  na.action=na.omit)
car::Anova(Irov)
df<-df.residual(Irov)
MSerror<-deviance(Irov)/df
comparison <- HSD.test(Irov,c("pop"),MSerror=MSerror,alpha=0.05/6,  group=TRUE)
comparison

#anova on relatedness: shinnecock
Irov <- lm(IR~pop, data=srel2,  na.action=na.omit)
car::Anova(Irov)
df<-df.residual(Irov)
MSerror<-deviance(Irov)/df
comparison <- HSD.test(Irov,c("pop"),MSerror=MSerror,alpha=0.05/4,  group=TRUE)
comparison

#####Giant bar plot####
#get the dataframe together. 

#HWE 
#atlhwe #all, wfhwe.pop #by bay, mthwe #mt cohorts, shihwe #shincohorts
#main stats
allstats <-as.data.frame(pivlocstats) # GRP, LOCUS, variable, value
#allelic richness - not including all because it screws up the graph proportions. 
rar <-bind_rows(meltar3,mtar,shar) %>%as.data.frame()%>%dplyr::rename(GRP=variable)%>%mutate(variable="Ar")
#private alleles
PA <- bind_rows(meltpA,meltMTpa,meltSHpa) %>%as.data.frame()%>%dplyr::rename(GRP=variable)%>%mutate(variable="PA")
#IR
#is not by locus, so you can't really do it. That's ok. 
allstats <-full_join(allstats,rar) %>% full_join(PA) %>%mutate(GRP=as.factor(GRP),LOCUS=as.factor(LOCUS), variable=as.factor(variable))

#you can actually tricky force IR to be in there. just rename the columns. 
rel3 <-dplyr::rename(rel3, LOCUS=Ind,GRP=pop, value=IR) %>% mutate(variable="IR", GRP=as.factor(GRP))
rel3 <-mutate_at(rel3,.vars=vars(GRP),.funs=forcats::fct_recode, ALL="Atl")
allstats <-full_join(allstats, rel3) %>%mutate(GRP=as.factor(GRP),LOCUS=as.factor(LOCUS), variable=as.factor(variable))

#allcolors <-c("grey","#d0d1e6","#a6bddb","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#1c9099","#016450","#016450","#016450","#016450","#016450")
allcolors <-c("grey","#d0d1e6","#a6bddb","#67a9cf","#1c9099","#016450","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#016450","#016450","#016450","#016450")

#Giant barplot
allstats <-filter(allstats,GRP %in% c("ALL","Nap","Mor","Jam","Shin","Shin_1_2016","Shin_2_2016","Shin_1_2017","Shin_2_2017","Mt","Mt_3_adults","Mt_4_adults","Mt_1_2015","Mt_2_2015","Mt_1_2016","Mt_2_2016"))
allstats$GRP <-fct_relevel(allstats$GRP, c("ALL","Jam","Mor","Mt","Nap","Shin","Mt_1_2015","Mt_2_2015","Mt_1_2016","Mt_2_2016", "Mt_3_adults","Mt_4_adults","Shin_1_2016","Shin_2_2016","Shin_1_2017","Shin_2_2017"))

allstats %>%
  mutate(GRP=as.factor(GRP)) %>%
  #filter(variable %in% c("SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis","PA","Ar", "Hdef")) %>%
  #filter(variable %in% c("SIMPSON_IDX", "EVENNESS","Ht","Fis","PA","Ar")) %>%
  filter(variable %in% c("SIMPSON_IDX", "EVENNESS","Ht","Fis","IR","Ar")) %>%
  ggplot(aes(x=fct_rev(GRP),y=value),fill=GRP)+
  geom_boxplot(aes(fill=GRP))+ 
  scale_fill_manual(name = "GRP",values = allcolors)+
  coord_flip()+ ylab("")+xlab("")+
  facet_wrap(~variable, scales="free_x", nrow=2)+ 
  theme(axis.text = element_text(size = 14, face="bold"),axis.title = element_text(size = 14,face="bold"),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),panel.spacing = unit(1, "lines"))+guides(fill = FALSE, colour = FALSE) 
ggsave('alllocstats17IR.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 14, height = 10)



### Global significance test####

#set level to test for friedmans BAY LEVEL 
fdata <- filter(loc_stats, GRP %in% c("Jam","Mor","Nap","Mt","Shin"))
#####Friedman's test for Global heterogeneity#####
#a <-friedman.test(Hs~GRP | LOCUS, data= fdata) #Hs - expected
b <-friedman.test(EVENNESS~GRP | LOCUS, data= fdata) #evenness
c<-friedman.test(Ht~GRP | LOCUS, data= fdata) #Ht - nei
#d<-friedman.test(Ho~GRP | LOCUS, data= fdata) #HO
e<-friedman.test(Fis~GRP | LOCUS, data= fdata) #Fis
f<-friedman.test(SIMPSON_IDX~GRP | LOCUS, data= fdata) 
g<-friedman.test(value~variable | LOCUS, data= meltar3) #rareified alleles #come back to this one.... 
h<-friedman.test(value~variable | LOCUS, data= meltpA) #private alleles


tests <-c("EVENNESS","Ht","Fis","SIMPSON_IDX","Ar","Pa")
pvals <-c(b$p.value, c$p.value, e$p.value, f$p.value, g$p.value, h$p.value)
glob <-as.data.frame(cbind(tests,pvals)) %>%mutate(pvals=as.numeric(as.character(pvals)))%>%
  mutate(Pr.adj=p.adjust(pvals, method="BH")) %>% mutate(sig=ifelse(Pr.adj < 0.05, "sig","not.sig"))

#set level to test for friedmans  SHINNECOCK LEVEL
fdata <- filter(loc_stats, GRP %in% c("Shin_1_2016","Shin_1_2017","Shin_2_2016","Shin_2_2017"))
#####Friedman's test for Global heterogeneity#####
a <-friedman.test(Hs~GRP | LOCUS, data= fdata) #Hs - expected
b <-friedman.test(EVENNESS~GRP | LOCUS, data= fdata) #evenness
c<-friedman.test(Ht~GRP | LOCUS, data= fdata) #Ht - nei
d<-friedman.test(Ho~GRP | LOCUS, data= fdata) #HO
e<-friedman.test(Fis~GRP | LOCUS, data= fdata) #Fis
f<-friedman.test(SIMPSON_IDX~GRP | LOCUS, data= fdata) 
g<-friedman.test(value~variable | LOCUS, data= shar) #rareified alleles #come back to this one.... 
h<-friedman.test(value~variable | LOCUS, data= meltSHpa) #private alleles

tests <-c("EVENNESS","Ht","Fis","SIMPSON_IDX","Ar","Pa")
pvals <-c(b$p.value, c$p.value, e$p.value, f$p.value, g$p.value, h$p.value)
glob <-as.data.frame(cbind(tests,pvals)) %>%mutate(pvals=as.numeric(as.character(pvals)))%>%
  mutate(Pr.adj=p.adjust(pvals, method="BH")) %>% mutate(sig=ifelse(Pr.adj < 0.05, "sig","not.sig"))

#set level to test for friedmans  MATTITUCK LEVEL
fdata <- filter(loc_stats, GRP %in% c("Mt_1_2016","Mt_1_2015","Mt_2_2016","Mt_2_2015","Mt_3_adults","Mt_4_adults"))
#####Friedman's test for Global heterogeneity#####
a <-friedman.test(Hs~GRP | LOCUS, data= fdata) #Hs - expected
b <-friedman.test(EVENNESS~GRP | LOCUS, data= fdata) #evenness
c<-friedman.test(Ht~GRP | LOCUS, data= fdata) #Ht - nei
d<-friedman.test(Ho~GRP | LOCUS, data= fdata) #HO
e<-friedman.test(Fis~GRP | LOCUS, data= fdata) #Fis
f<-friedman.test(SIMPSON_IDX~GRP | LOCUS, data= fdata) 
g<-friedman.test(value~variable | LOCUS, data= mtar) #rareified alleles #come back to this one.... 
h<-friedman.test(value~variable | LOCUS, data= meltMTpa) #private alleles

tests <-c("EVENNESS","Ht","Fis","SIMPSON_IDX","Ar","Pa")
pvals <-c(b$p.value, c$p.value, e$p.value, f$p.value, g$p.value, h$p.value)
glob <-as.data.frame(cbind(tests,pvals)) %>%mutate(pvals=as.numeric(as.character(pvals)))%>%
  mutate(Pr.adj=p.adjust(pvals, method="BH")) %>% mutate(sig=ifelse(Pr.adj < 0.05, "sig","not.sig"))



#####For the friedman's tests that came out significant#####

#Instead of loc_stats, use allstats and pivot it
allstats2 <-pivot_wider(allstats,id_cols=c("GRP","LOCUS"), names_from=variable, values_from=value)
allstats2 <-as.data.frame(allstats2)

#LI BAYS
llocstats <- filter(allstats2, GRP %in% c("Jam","Mor","Nap","Shin","Mt")) ##CHANGE HERE
#list of all the different possible pairs
comp <- as.character(unique(llocstats$GRP))  
pairs <- expand.grid(comp, comp) %>%
  filter(!Var1 == Var2) %>%
  rownames_to_column("PAIR") %>%
  split(.$PAIR) %>%
  purrr::map(function(x){
    x %>%
      dplyr::select(-PAIR) %>%
      gather(key = llocstats, value = GRP, 1:2) %>%
      dplyr::select(-llocstats)})

# empty data frame for results
results <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("pop1", "pop2", "stat", "temp", "p.value", "test"))
n <- as.numeric(length(pairs))

library("coin")
# loop over pairs
for(p in 1:n){
  pair <- pairs[[p]]$GRP
  temp <- llocstats %>%
    dplyr::filter(GRP %in% pair) %>%
    mutate(GRP = ordered(GRP, levels = pair),
           LOCUS = as.factor(LOCUS)) %>%
    droplevels()
  wilcox <- wilcoxsign_test(EVENNESS ~ GRP | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,
                            zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)),
                   "test" = "EVENNESS")  ### change this one. 
  results <- bind_rows(results, df)}


splitres <- dplyr::rename(results, Pr.exact=p.value)%>% dplyr::select(-temp)%>% split(.$test)%>%map_dfr(pradj)
write.csv(splitres,file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/wilcoxOCt29.csv" )

#we may not need to make heatmaps. I am not going to do them for now. but the code is in Microsats 2020. 

#MATTITUCK
llocstats <- filter(allstats2, GRP %in% c("Mt_1_2016","Mt_1_2015","Mt_2_2016","Mt_2_2015","Mt_3_adults","Mt_4_adults")) ##CHANGE HERE
#list of all the different possible pairs
comp <- as.character(unique(llocstats$GRP))  
pairs <- expand.grid(comp, comp) %>%
  filter(!Var1 == Var2) %>%
  rownames_to_column("PAIR") %>%
  split(.$PAIR) %>%
  purrr::map(function(x){
    x %>%
      dplyr::select(-PAIR) %>%
      gather(key = llocstats, value = GRP, 1:2) %>%
      dplyr::select(-llocstats)})

# empty data frame for results
results <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("pop1", "pop2", "stat", "temp", "p.value", "test"))
n <- as.numeric(length(pairs))

library("coin")
# loop over pairs
for(p in 1:n){
  pair <- pairs[[p]]$GRP
  temp <- llocstats %>%
    dplyr::filter(GRP %in% pair) %>%
    mutate(GRP = ordered(GRP, levels = pair),
           LOCUS = as.factor(LOCUS)) %>%
    droplevels()
  wilcox <- wilcoxsign_test(PA ~ GRP | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,
                            zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)),
                   "test" = "PA")  ### change this one. 
  results <- bind_rows(results, df)}


splitres <- dplyr::rename(results, Pr.exact=p.value)%>% dplyr::select(-temp)%>% split(.$test)%>%map_dfr(pradj)
write.csv(splitres,file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/wilcoxOCt29MT.csv" )

#SHINNECOCK
llocstats <- filter(allstats2, GRP %in% c("Shin_1_2016","Shin_2_2016","Shin_1_2017","Shin_2_2017")) ##CHANGE HERE
#list of all the different possible pairs
comp <- as.character(unique(llocstats$GRP))  
pairs <- expand.grid(comp, comp) %>%
  filter(!Var1 == Var2) %>%
  rownames_to_column("PAIR") %>%
  split(.$PAIR) %>%
  purrr::map(function(x){
    x %>%
      dplyr::select(-PAIR) %>%
      gather(key = llocstats, value = GRP, 1:2) %>%
      dplyr::select(-llocstats)})

# empty data frame for results
results <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("pop1", "pop2", "stat", "temp", "p.value", "test"))
n <- as.numeric(length(pairs))

library("coin")
# loop over pairs
for(p in 1:n){
  pair <- pairs[[p]]$GRP
  temp <- llocstats %>%
    dplyr::filter(GRP %in% pair) %>%
    mutate(GRP = ordered(GRP, levels = pair),
           LOCUS = as.factor(LOCUS)) %>%
    droplevels()
  wilcox <- wilcoxsign_test(SIMPSON_IDX ~ GRP | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,
                            zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)),
                   "test" = "SIMPSON_IDX")  ### change this one. 
  results <- bind_rows(results, df)}


splitres <- dplyr::rename(results, Pr.exact=p.value)%>% dplyr::select(-temp)%>% split(.$test)%>%map_dfr(pradj)
write.csv(splitres,file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/wilcoxOCt29SHI.csv" )


### Genetic distance ####
#### Genetic Distance ####

## Tree's using provesti's distance

#all bays?
set.seed(999)
wfpopLD %>%
  genind2genpop(pop = ~Bay) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = provesti.dist, missingno="ignore")
#This one makes a lot of sense

#Mattituck
set.seed(999)
setpop(mtco) <-~Bay/Con/Year
mtco %>%
  genind2genpop(pop = ~Bay/Con/Year) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = provesti.dist, missingno="ignore")

#Shinnecock
set.seed(999)
setPop(shinco) <-~Bay/Con/Year
shinco %>%
  genind2genpop(pop = ~Bay/Con/Year) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = provesti.dist, missingno="ignore")

### Alpha & Beta diversity ####

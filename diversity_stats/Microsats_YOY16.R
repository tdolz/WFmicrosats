### Genetic analysis Bays ###
## created 5-30-23 From Microsats_20.R
## updated to create plots in the microsat_figs_yoy16.R script on 1/24/23
## 

## NOTES ##
## I THINK THERE IS AN ISSUE HERE WITH THE MATTITUCK AND SHINNECOCK DATA BECAUSE WE REMOVED THE DATA OTHER THAN 2016 YOY EVERYWHERE... 
## SOOO THAT WILL NEED TO BE FIXED BUT CAN EASILY BE FIXED. 

#install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)
library("tidyverse")
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
library("coin")
library("reshape2")
library("strataG")
library("hierfstat")
library("viridis")
library("forcats")
library("data.table")


#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")


##### Formating the dataset #####
# We are going to use the doubl0 version. 
##### Formating the dataset #####
wfpop <- read.genalex("./data/popcorrect_17_sept20204genalex_doubl0ABC.csv")
wfpop4df <-read.csv("./data/popcorrect_17_sept2020_doubl0ABC.csv", header = TRUE) #csv

splitStrata(wfpop) <-~Ocean/Bay/Con/Year
setPop(wfpop) <-~Bay/Year/Con
wfpop <-popsub(wfpop, exclude=c("Mt_2015_5", "Mt_2016_5"))
setPop(wfpop) <-~Bay

#do data cleanup first- on full dataset. 
## Data cleanup. For details see datacleanup_comparison.R ##
wfpopLD <-genclone2genind(wfpop)
all_loci <- locNames(wfpopLD)
removeloc <- c("WF06","WF27")# we also may want to remove A441 and WF32 for missing data reasons. 
#especially WF32 because it's 70% missing in moriches and is out of HWE.
#Come back to this
keeploc <- setdiff(all_loci, removeloc)
wfpopLD <- wfpopLD[loc = keeploc]
wfpopLD <- wfpopLD %>% missingno("geno", cutoff=0.25) 
wfpopLD <-wfpopLD %>% missingno("loci", cutoff=0.25) #removes WF01
setPop(wfpopLD) <-~Bay
info_table(wfpopLD, plot = TRUE, scaled =FALSE)

#Now subset to 2016 YOY - and save a copy of the original. 
og.wfpopLD <-wfpopLD
##Subset to YOY16 only##
#remove all individuals that aren't 2016 YOY
setPop(wfpopLD) <-~Bay/Year/Con
unique(wfpopLD@pop)
wfpopLD <-popsub(wfpopLD, exclude=c("Mt_2015_1","Mt_2015_2","Mt_adults_3", "Mt_adults_4", "Shin_2017_1","Shin_2017_2"))
unique(wfpopLD@pop)


###CREATE DATASET FOR THE no 2015 COMPARISON ####
#but retain the adults and Shin 2017
#remove mt_2015. Also remove unclassified adults
setPop(og.wfpopLD)<-~Bay/Year/Con
wfpop.mtshi <-popsub(og.wfpopLD, exclude=c("Mt_2015_1","Mt_2015_2","Mt_2015_5","Mt_2016_5"))

#How many loci out of HWE? - all YOY 2016
setPop(wfpopLD) <-~Ocean
hwe.ocean <-hw.test(wfpopLD, B=1000)%>%as.data.frame()%>% 
  mutate(pop="Atl")%>%
  mutate(pval.adj =p.adjust(`Pr(chi^2 >)`, method="BH"),p.chi =ifelse(`Pr(chi^2 >)`>0.05,"in","out"))%>%
  mutate(p.chi.adj = ifelse(pval.adj >0.05,"in","out"))
# loci out of HWE - by bay 2016
setPop(wfpopLD)<-~Bay
wfhwe.pop <- seppop(wfpopLD) %>% lapply(hw.test, B=1000)
for(i in 1:length(wfhwe.pop)){
  wfhwe.pop[[i]]<-as.data.frame(wfhwe.pop[[i]])%>%
    mutate(pop=names(wfhwe.pop)[i])%>%
    mutate(pval.adj =p.adjust(`Pr(chi^2 >)`, method="BH"),p.chi =ifelse(`Pr(chi^2 >)`>0.05,"in","out"))%>%
    mutate(p.chi.adj = ifelse(pval.adj >0.05,"in","out"))
  }
wfhwe.pop <-bind_rows(wfhwe.pop)%>%bind_rows(hwe.ocean)

###HWE IN MATTITUCK CREEK & SHINNECOCK BAY ################
setPop(wfpop.mtshi)<-~Bay/Year/Con
mtshi <- seppop(wfpop.mtshi) %>% lapply(hw.test, B=1000)
for(i in 1:length(mtshi)){
  mtshi[[i]]<-as.data.frame(mtshi[[i]])%>%
    mutate(pop=names(mtshi)[i])%>%
    mutate(pval.adj =p.adjust(`Pr(chi^2 >)`, method="BH"),p.chi =ifelse(`Pr(chi^2 >)`>0.05,"in","out"))%>%
    mutate(p.chi.adj = ifelse(pval.adj >0.05,"in","out"))
}
wfhwe.pop <-bind_rows(mtshi)%>%bind_rows(wfhwe.pop)
wfhwe.pop %>% group_by(pop)%>%summarize(sig=sum(pval.adj <= 0.05))

write.csv(wfhwe.pop, file="./diversity_stats/diversity_output_files/YOY16_bay/wfhwepop16.csv")

#PopGenReport - Remember to turn this off. 
#setPop(wfpopLD) <-~Bay
#wf.gen <-genclone2genind(wfpopLD) 
#popgenreport(wf.gen,mk.counts=TRUE,mk.locihz = TRUE, mk.fst=TRUE, mk.allele.dist=TRUE, mk.null.all=TRUE,mk.allel.rich = TRUE,mk.differ.stats = TRUE,path.pgr=getwd(),mk.Rcode=TRUE,mk.pdf=TRUE )

###Summary stats from popgenreport. 
##### Shannon's summary stats ########
#https://gist.github.com/sjoleary/cdc32efbfd50c96eef446ebb7c2f2387

#dataframe with sample information
setPop(wfpopLD.og) <- ~Bay/Con/Year
df <-genind2df(wfpopLD.og, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(ncol(df),1:(ncol(df)-1))]
df[df=="NA"] <- 0 # missing data must be 0
SampleInfo <- dplyr::select(df, Ind, pop)
SampleInfo <- separate(SampleInfo, pop, c("Bay","Con","Year"))
SampleInfo <-mutate(SampleInfo, Ocean="Atl")
SampleInfo <-SampleInfo[,c(1,5,2,3,4)]

setPop(wfpopLD)<-~Ocean
gen_oce <-seppop(wfpopLD) #ALL for the bays should be only YOY2016

#update - bays with 2016 YOY only 
setPop(wfpopLD)<-~Bay
gen_bay <-seppop(wfpopLD) #here bays are all ONLY 2016 YOY

#update Mattituck only (no 2015, and no unidentified adults)
setPop(wfpop.mtshi)<-~Bay/Year/Con
genMt <-popsub(wfpop.mtshi, sublist=c("Mt_2016_1","Mt_2016_2","Mt_adults_3", "Mt_adults_4"))
gen_Mt <-seppop(genMt)

#update Shinnecock only
genShi <-popsub(wfpop.mtshi, sublist=c("Shin_2016_1", "Shin_2017_1","Shin_2016_2", "Shin_2017_2"))
gen_Shi <-seppop(genShi)

gen_grp <-c(gen_oce,
            gen_bay,
            gen_Mt,
            gen_Shi)
gen_grp[["ALL"]] <-og.wfpopLD

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
  dat <- hierfstat:::genind2hierfstat(gen_grp[[p]])
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

loc_stats[is.na(loc_stats)] <- NA
write.csv(loc_stats,file="./diversity_stats/diversity_output_files/YOY16_bay/loc_stats17.csv")
pivlocstats <-pivot_longer(loc_stats,cols = c("N_ALLELES", "SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis","SHANNON_IDX","STODD_TAYLOR_IDX"),
                             names_to="variable", values_to="value")
sumlocstats <-ddply(pivlocstats, GRP~variable, summarize, meanvar=mean(value), sdvar=sd(value))
write.csv(pivlocstats,file="./diversity_stats/diversity_output_files/YOY16_bay/pivlocstats17.csv")
write.csv(sumlocstats,file="./diversity_stats/diversity_output_files/YOY16_bay/sumlocstats17.csv")
#############################################################################################################

#create different loc stats groups for testing. 
loc_stats_bay <-filter(loc_stats, GRP %in% c("Nap","Mor","Jam","Shin","Mt"))
meltlocstats_bay <-pivot_longer(loc_stats_bay,cols = c("N_ALLELES", "SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis","SHANNON_IDX","STODD_TAYLOR_IDX"),
                                names_to="variable", values_to="value")
write.csv(meltlocstats_bay, file="./diversity_stats/diversity_output_files/YOY16_bay/meltlocstats_bay.csv")


meltlocstats2 <-filter(loc_stats, GRP %in% c("Nap","Mor","Jam","Shin","Mt","ALL"))%>%
  pivot_longer(cols = c("N_ALLELES", "SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis","SHANNON_IDX","STODD_TAYLOR_IDX"),
               names_to="variable", values_to="value")
#this is essentially by bay but also includes the "ALL" Group


#we have some reason to distrust the fact that ht and hs are the same after shannons code but not globally. so we will also do this another way.
#if you look at this, basic stats does not work at the subpopulation level and this could be why. 
g2h_all <-genind2hierfstat(wfpopLD)
g2hall <-basic.stats(g2h_all)
global <-as.data.frame(g2hall$perloc) %>% tibble::rownames_to_column("locus") %>% mutate(Bay="ALL")
write.csv(global,file="./diversity_stats/diversity_output_files/YOY16_bay/basic_stats17.csv" )


######## CALCULATE RAREFIED ALLELIC RICHNESS ----
# differences in sample size can bias the number of alleles sampled in a population
# calculate allelic richness corrected for sample size using rarefaction

# overall ====
setPop(wfpopLD) <- ~Ocean
dat <- hierfstat:::genind2hierfstat(wfpopLD)
ar <- allelic.richness(dat,diploid = TRUE)

ar <- as.data.frame(ar$Ar) %>%
  rownames_to_column("LOCUS") %>%
  dplyr::rename(ALL = 2)%>% 
  select(-dumpop)

# By Bay ===
# by region ====
setPop(wfpopLD) <- ~Bay
dat <- hierfstat:::genind2hierfstat(wfpopLD, pop = wfpopLD@pop)
df <- allelic.richness(dat,diploid = TRUE)

# I am not sure how to tell which one is which in order to label the columns, but I assume they're in the same order as WFPOPLD?
df <- as.data.frame(df$Ar) %>%
  rownames_to_column("LOCUS") #%>%dplyr::rename(Nap = V1,Mor = V2,Jam = V3,Shin = V4,Mt = V5)
ar <- left_join(ar, df)
write.csv(ar, "./diversity_stats/diversity_output_files/YOY16_bay/rarefied_allelecount.csv")

## melt AR MATTITUCK $ SHINNECOCK no yoy 2015
setPop(wfpop.mtshi) <- ~Bay/Year/Con
dat <- hierfstat:::genind2hierfstat(wfpop.mtshi, pop = wfpop.mtshi@pop)
df <- allelic.richness(dat,diploid = TRUE)

# I am not sure how to tell which one is which in order to label the columns, but I assume they're in the same order as WFPOPLD?
df <- as.data.frame(df$Ar) %>%
  rownames_to_column("LOCUS") #%>%dplyr::rename(Nap = V1,Mor = V2,Jam = V3,Shin = V4,Mt = V5)
ar_mtshi <- left_join(ar, df)
write.csv(ar_mtshi, "./diversity_stats/diversity_output_files/YOY16_bay/rarefied_allelecount_ALLGROUPS.csv")

meltarMT <-pivot_longer(ar_mtshi, cols=c("Mt_2016_1","Mt_2016_2","Mt_adults_3", "Mt_adults_4"),names_to="variable", values_to="value")
meltarSHIN <-pivot_longer(ar_mtshi, cols=c("Shin_2016_1", "Shin_2017_1", "Shin_2016_2", "Shin_2017_2"),names_to="variable", values_to="value")

#for rarified alleles, you have to make a new Allelic Richness analysis for each group you are comparing. 
meltar2 <- pivot_longer(ar, cols=c("Mt","Shin","Nap","Mor","Jam","ALL"),names_to="variable", values_to="value")
meltar3 <-filter(meltar2, variable != "ALL") #it doenst like it when you exclude the "all" category... 
meltar_all <-filter(meltar2, variable == "ALL")
write.csv(meltar3, file="./diversity_stats/diversity_output_files/YOY16_bay/meltar3.csv")

###################################Private alleles #############################################

#Between Bays
df <-genind2df(wfpopLD, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(36,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)]
df[df=="NA"] <- 0 # missing data must be 0
wf.g <-df2gtypes(df,ploidy=2)

library("data.table")
pA <-as.data.frame(privateAlleles(wf.g)) 
setDT(pA,keep.rownames=TRUE)
colnames(pA)[1] <- "LOCUS"
write.csv(pA, "./diversity_stats/diversity_output_files/YOY16_bay/private_allelecount.csv")

## melt AR MATTITUCK $ SHINNECOCK no yoy 2015 - 
#do mattituck and Shinnecock seperately, separate from each other. 
setPop(wfpop.mtshi) <- ~Bay/Year/Con
#mattituck
Mtpop <- popsub(wfpop.mtshi, sublist=c("Mt_2016_1","Mt_2016_2","Mt_adults_3", "Mt_adults_4"))
df <-genind2df(Mtpop, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(36,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)]
df[df=="NA"] <- 0 # missing data must be 0
wf.g <-df2gtypes(df,ploidy=2)
pA_mt <-as.data.frame(privateAlleles(wf.g)) 
setDT(pA_mt,keep.rownames=TRUE)
colnames(pA_mt)[1] <- "LOCUS"

#shinnecock
Shipop <- popsub(wfpop.mtshi, sublist = c("Shin_2016_1", "Shin_2017_1", "Shin_2016_2", "Shin_2017_2"))
df <-genind2df(Shipop, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(36,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)]
df[df=="NA"] <- 0 # missing data must be 0
wf.g <-df2gtypes(df,ploidy=2)
pA_shi <-as.data.frame(privateAlleles(wf.g)) 
setDT(pA_shi,keep.rownames=TRUE)
colnames(pA_shi)[1] <- "LOCUS"

pA_mtshi <-bind_rows(pA_mt, pA_shi)
write.csv(pA_mtshi, "./diversity_stats/diversity_output_files/YOY16_bay/PrivateAlleles_ALLGROUPS.csv")

meltPA_MT <-pivot_longer(pA_mtshi, cols=c("Mt_2016_1","Mt_2016_2","Mt_adults_3", "Mt_adults_4"),names_to="variable", values_to="value")
meltPA_SHIN <-pivot_longer(pA_mtshi, cols=c("Shin_2016_1", "Shin_2017_1", "Shin_2016_2", "Shin_2017_2"),names_to="variable", values_to="value")
meltpA <- pivot_longer(pA, cols=c("Mt","Shin","Nap","Mor","Jam"),names_to="variable", values_to="value")


#################### Inbreeding - Internal Relatedness ##############################
library("Rhh")
setPop(wfpopLD) <-~Bay
#wfp <-genind2df(wfpopLD)
wfp <-genind2df(wfpopLD, usepop = TRUE, oneColPerAll = TRUE) %>% tibble::rownames_to_column("Ind")
wfp <-mutate_at(wfp,vars(-pop, -Ind),as.numeric)
wfp <-mutate_at(wfp,vars(-pop, -Ind),as.numeric)
wfp[is.na(wfp)] <- 0

rel <- dplyr::select(wfp, -1, -2)
rel <-ir(rel)
#rel <-ir(wfp[,-1])

#rel
rel2 <-as.data.frame(rel)
names(rel2) <-c("IR")
rel2 <-cbind(wfp,rel2) 
rel2 <-dplyr::select(rel2,Ind,pop,IR)

##### Internal relatedness mattituck ###########
setPop(Mtpop) <-~Bay/Year/Con
wfp <-genind2df(Mtpop, usepop = TRUE, oneColPerAll = TRUE) %>% tibble::rownames_to_column("Ind")
wfp <-mutate_at(wfp,vars(-pop, -Ind),as.numeric)
wfp <-mutate_at(wfp,vars(-pop, -Ind),as.numeric)
wfp[is.na(wfp)] <- 0
rel <- dplyr::select(wfp, -1, -2)
rel <-ir(rel)
relMt <-as.data.frame(rel)
names(relMt) <-c("IR")
relMt <-cbind(wfp,relMt) 
relMt <-dplyr::select(relMt,Ind,pop,IR)

##### Internal relatedness Shinnecock ###########
setPop(Shipop) <-~Bay/Year/Con
wfp <-genind2df(Shipop, usepop = TRUE, oneColPerAll = TRUE) %>% tibble::rownames_to_column("Ind")
wfp <-mutate_at(wfp,vars(-pop, -Ind),as.numeric)
wfp <-mutate_at(wfp,vars(-pop, -Ind),as.numeric)
wfp[is.na(wfp)] <- 0
rel <- dplyr::select(wfp, -1, -2)
rel <-ir(rel)
relShi <-as.data.frame(rel)
names(relShi) <-c("IR")
relShi <-cbind(wfp,relShi) 
relShi <-dplyr::select(relShi,Ind,pop,IR)

rel2 <-bind_rows(rel2,relMt,relShi)
ddply(rel2, ~pop, summarize, meanir=mean(IR))

write.csv(rel2, "./diversity_stats/diversity_output_files/YOY16_bay/IR.csv")

######################################### TABLE 2 ###################################################################
#we will try doing this again later. 
#t2bays <-filter(sumlocstats,GRP %in% c("Atl","Jam","Mor","Mt","Nap","Shin"))%>%pivot_wider(names_from = "variable", values_from = c("meanvar","sdvar"))
#t2bays <-t(t2bays)
#t2bays <-t2bays[c("meanvar_Ho","sdvar_Ho","meanvar_Hs","sdvar_Hs","meanvar_Ht","sdvar_Ht"),]

######## TEST FOR SIGNIFICANT DIFFERENCES#####
# need to use Friedman's test for global test and Wilcoxon signed rank for pairwise tests 
# to test symmetry of numeric repeated measurements (stastic per locus) in block design
# Example using gene diversity (expected heterozygosity by estuary)

#set level to test for friedmans
fdata <- loc_stats_bay #Is the same here as loc_stats_16 because we got rid of the other individuals at the beginning. 


#####Friedman's test for Global heterogeneity#####
a<-friedman.test(Hs~GRP | LOCUS, data= fdata) #Hs - fixing this. 
b<-friedman.test(EVENNESS~GRP | LOCUS, data= fdata) #evenness
c<-friedman.test(Ht~GRP | LOCUS, data= fdata) #Ht
d<-friedman.test(Fis~GRP | LOCUS, data= fdata) #Fis
e<-friedman.test(SHANNON_IDX~GRP | LOCUS, data= fdata) 
f<-friedman.test(value~variable | LOCUS, data= meltar3) #rareified alleles #come back to this one.... 
g<-friedman.test(hobs~Bay | locus, data=h4)
h <-friedman.test(value~variable | LOCUS, data= meltpA) #private alleles
i <-friedman.test(SIMPSON_IDX~GRP| LOCUS, data= fdata)

p.vals <-c(a$p.value, b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value, h$p.value, i$p.value)
comparisons <-c(a$data.name,b$data.name,c$data.name,d$data.name,e$data.name,f$data.name,g$data.name, h$data.name, i$data.name)

friedmans_tests <-cbind(comparisons, p.vals)%>%as.data.frame%>%mutate(p.vals=as.numeric(p.vals))
friedmans_tests$chi_sq <-c(a$statistic, b$statistic,c$statistic,d$statistic,e$statistic,f$statistic,g$statistic, h$statistic, i$statistic)
friedmans_tests <-mutate(friedmans_tests, pval_adj =p.adjust(p.vals, method="BH"))
write.csv(friedmans_tests,"./diversity_stats/diversity_output_files/YOY16_bay/friedmanstests.csv")

########### Friedman test Mattituck (NO 2015 yoy) ########################
loc_stats_MT <- filter(loc_stats, GRP %in% c("Mt_2016_1","Mt_2016_2","Mt_adults_3","Mt_adults_4"))
fdata <- loc_stats_MT #Is the same here as loc_stats_16 because we got rid of the other individuals at the beginning. 

#####Friedman's test for Global heterogeneity#####
a<-friedman.test(Hs~GRP | LOCUS, data= fdata) #Hs - fixing this. 
b<-friedman.test(EVENNESS~GRP | LOCUS, data= fdata) #evenness
c<-friedman.test(Ht~GRP | LOCUS, data= fdata) #Ht
d<-friedman.test(Fis~GRP | LOCUS, data= fdata) #Fis
e<-friedman.test(SHANNON_IDX~GRP | LOCUS, data= fdata) 
i <-friedman.test(SIMPSON_IDX~GRP| LOCUS, data= fdata)

f<-friedman.test(value~variable | LOCUS, data= meltarMT) #rareified alleles 
h <-friedman.test(value~variable | LOCUS, data= meltPA_MT) #private all

p.vals <-c(a$p.value, b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,h$p.value, i$p.value)
comparisons <-c(a$data.name,b$data.name,c$data.name,d$data.name,e$data.name,f$data.name, h$data.name, i$data.name)

friedmans_tests <-cbind(comparisons, p.vals)%>%as.data.frame%>%mutate(p.vals=as.numeric(p.vals))
friedmans_tests$chi_sq <-c(a$statistic, b$statistic,c$statistic,d$statistic,e$statistic,f$statistic,h$statistic, i$statistic)
friedmans_tests <-mutate(friedmans_tests, pval_adj =p.adjust(p.vals, method="BH"))
write.csv(friedmans_tests,"./diversity_stats/diversity_output_files/YOY16_bay/friedmanstests_MTno2015.csv")


########### Friedman test Shinnecock, Should be the same as before ########################
loc_stats_SHI <- filter(loc_stats, GRP %in% c("Shin_2016_1", "Shin_2017_1", "Shin_2016_2", "Shin_2017_2"))
fdata <- loc_stats_SHI #Is the same here as loc_stats_16 because we got rid of the other individuals at the beginning. 

#####Friedman's test for Global heterogeneity#####
a<-friedman.test(Hs~GRP | LOCUS, data= fdata) #Hs - fixing this. 
b<-friedman.test(EVENNESS~GRP | LOCUS, data= fdata) #evenness
c<-friedman.test(Ht~GRP | LOCUS, data= fdata) #Ht
d<-friedman.test(Fis~GRP | LOCUS, data= fdata) #Fis
e<-friedman.test(SHANNON_IDX~GRP | LOCUS, data= fdata) 
i <-friedman.test(SIMPSON_IDX~GRP| LOCUS, data= fdata)

f<-friedman.test(value~variable | LOCUS, data= meltarSHIN) #rareified alleles 
h <-friedman.test(value~variable | LOCUS, data= meltPA_SHIN) #private all

p.vals <-c(a$p.value, b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,h$p.value, i$p.value)
comparisons <-c(a$data.name,b$data.name,c$data.name,d$data.name,e$data.name,f$data.name, h$data.name, i$data.name)

friedmans_tests <-cbind(comparisons, p.vals)%>%as.data.frame%>%mutate(p.vals=as.numeric(p.vals))
friedmans_tests$chi_sq <-c(a$statistic, b$statistic,c$statistic,d$statistic,e$statistic,f$statistic,h$statistic, i$statistic)
friedmans_tests <-mutate(friedmans_tests, pval_adj =p.adjust(p.vals, method="BH"))
write.csv(friedmans_tests,"./diversity_stats/diversity_output_files/YOY16_bay/friedmanstests_shin.csv")




####################################################################################################################
######################Wilcoxon tests #####
#tell us which groups are significantly different. 
# you have to set which level you are comparing. 
  #the different levels we have to work with are:
    #loc_stats_bay <- bays
    #loc_stats_BYC <- bay year con
    #loc_stats_MT <-bay con for mattituck only, no group 5 (unidentified adults)
    #loc_stats_shin <- bay year con shin
    #loc_stats_16 <- 2016 YOY at bay level. 

#################### WILCOX TEST YOY 2016 BAY #######################################################


llocstats <- loc_stats_bay ## This IS 2016 yoy at bay level because we previously removed the other groups. 
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

library("coin")
### I tried for days to set this up in a double loop. It seems like the issue is that
### you can't paste the column name into wilcoxsign_test to refer to it from the loop.
### even if it's the exact same characters and type. So I am just going to repeat the 
### code for each stat for now. I've wasted too much time on this. 

n <- as.numeric(length(pairs))

#NEI
results_list <-list()
for(p in 1:n){
  #tryCatch({
  pair <- pairs[[p]]$GRP
  temp <- llocstats %>%
    dplyr::filter(GRP %in% pair) %>%
    mutate(GRP = ordered(GRP, levels = pair),
           LOCUS = as.factor(LOCUS)) %>%
    droplevels()
  wilcox <- wilcoxsign_test(Ht ~ GRP | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)),
                   "test" = "Ht") 
  #if(!is.null(df)){
  results_list[[p]] <- df #}#, error=function(e){})
}
#results_list <-results_list[-which(sapply(results_list, is.null))]
results_nei<-bind_rows(results_list)%>%mutate(padj=p.adjust(p.value, method="BH"))

#Fis
results_list <-list()
for(p in 1:n){
  #tryCatch({
  pair <- pairs[[p]]$GRP
  temp <- llocstats %>%
    dplyr::filter(GRP %in% pair) %>%
    mutate(GRP = ordered(GRP, levels = pair),
           LOCUS = as.factor(LOCUS)) %>%
    droplevels()
  wilcox <- wilcoxsign_test(Fis ~ GRP | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)),
                   "test" = "Fis") 
  #if(!is.null(df)){
  results_list[[p]] <- df #}#, error=function(e){})
}
#results_list <-results_list[-which(sapply(results_list, is.null))]
results_fis<-bind_rows(results_list)%>%mutate(padj=p.adjust(p.value, method="BH"))

#EVENNESS
results_list <-list()
for(p in 1:n){
  #tryCatch({
  pair <- pairs[[p]]$GRP
  temp <- llocstats %>%
    dplyr::filter(GRP %in% pair) %>%
    mutate(GRP = ordered(GRP, levels = pair),
           LOCUS = as.factor(LOCUS)) %>%
    droplevels()
  wilcox <- wilcoxsign_test(EVENNESS ~ GRP | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)),
                   "test" = "EVENNESS") 
  #if(!is.null(df)){
  results_list[[p]] <- df #}#, error=function(e){})
}
#results_list <-results_list[-which(sapply(results_list, is.null))]
results_even<-bind_rows(results_list)%>%mutate(padj=p.adjust(p.value, method="BH"))

#SIMPSON_IDX
results_list <-list()
for(p in 1:n){
  #tryCatch({
  pair <- pairs[[p]]$GRP
  temp <- llocstats %>%
    dplyr::filter(GRP %in% pair) %>%
    mutate(GRP = ordered(GRP, levels = pair),
           LOCUS = as.factor(LOCUS)) %>%
    droplevels()
  wilcox <- wilcoxsign_test(SIMPSON_IDX ~ GRP | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)),
                   "test" = "SIMPSON_IDX") 
  #if(!is.null(df)){
  results_list[[p]] <- df #}#, error=function(e){})
}
#results_list <-results_list[-which(sapply(results_list, is.null))]
results_simpson<-bind_rows(results_list)%>%mutate(padj=p.adjust(p.value, method="BH"))


##Wilcoxon test for rarified alleles##
#list of all the different possible pairs
comp <- as.character(unique(meltar3$variable))  
pairs <- expand.grid(comp, comp) %>%
  filter(!Var1 == Var2) %>%
  rownames_to_column("PAIR") %>%
  split(.$PAIR) %>%
  purrr::map(function(x){
    x %>%
      dplyr::select(-PAIR) %>%
      gather(key = meltar3, value = variable, 1:2) %>%
      dplyr::select(-meltar3)})
meltar3 <-as.data.frame(meltar3) %>% mutate(GRP=as.factor(variable))
# empty data frame for results
results_ar <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("pop1", "pop2", "stat", "temp", "p.value"))
n <- as.numeric(length(pairs))

library("coin")
# loop over pairs
results_ar <-list()
for(p in 1:n){
  pair <- pairs[[p]]$variable
  temp <- meltar3 %>%
    dplyr::filter(variable %in% pair) %>%
    mutate(GRP = ordered(variable, levels = pair),
           LOCUS = as.factor(LOCUS)) %>%
    droplevels()
  wilcox <- wilcoxsign_test(value ~ GRP | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,
                            zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)))
  results_ar[[p]] <-df}
results_ar <- bind_rows(results_ar)%>%mutate(test="AR",padj=p.adjust(p.value, method="BH"))

#results_ar <-arrange(results_ar, test.statistic) %>% unique() %>% 
  #mutate(pair1 = pmin(pop1,pop2), pair2 =pmax(pop1,pop2)) %>% arrange(pair1, pair2)
  #

#### Wilcoxon test for private alleles
#list of all the different possible pairs
comp <- as.character(unique(meltpA$variable))  
pairs <- expand.grid(comp, comp) %>%
  filter(!Var1 == Var2) %>%
  rownames_to_column("PAIR") %>%
  split(.$PAIR) %>%
  purrr::map(function(x){
    x %>%
      dplyr::select(-PAIR) %>%
      gather(key = meltpA, value = variable, 1:2) %>%
      dplyr::select(-meltpA)})
meltpA <-as.data.frame(meltpA) %>% mutate(GRP=as.factor(variable))
# empty data frame for results
results_pA <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("pop1", "pop2", "stat", "temp", "p.value"))
n <- as.numeric(length(pairs))

library("coin")
# loop over pairs
results_pA <-list()
for(p in 1:n){
  pair <- pairs[[p]]$variable
  temp <- meltpA %>%
    dplyr::filter(variable %in% pair) %>%
    mutate(GRP = ordered(variable, levels = pair),
           LOCUS = as.factor(LOCUS)) %>%
    droplevels()
  wilcox <- wilcoxsign_test(value ~ GRP | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,
                            zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)))
  results_pA[[p]] <-df}
results_pA <- bind_rows(results_pA)%>%mutate(test="PA", padj=p.adjust(p.value, method="BH"))

#combine
results <-bind_rows(results_nei, results_fis, results_even, results_simpson, results_ar, results_pA)
results <-mutate(results, significance = ifelse(padj>=0.05,"not significant","significant"))

write.csv(results,file="./diversity_stats/diversity_output_files/YOY16_bay/wilcox.csv")

library(multcompView)


############################################ STOPPED HERE #########################################
#try this letter grouping
results_nei2 <- unite(results_nei, grp, pop1, pop2)

multcompView::multcompLetters2(x=grp,data=results_nei2$padj)

############################# WILCOX TESTS WITHIN BAYS MT AND SHIN #########################################
# Private alleles in Mattituck
#### Wilcoxon test for private alleles
#list of all the different possible pairs
comp <- as.character(unique(meltPA_MT$variable))  
pairs <- expand.grid(comp, comp) %>%
  filter(!Var1 == Var2) %>%
  rownames_to_column("PAIR") %>%
  split(.$PAIR) %>%
  purrr::map(function(x){
    x %>%
      dplyr::select(-PAIR) %>%
      gather(key = meltPA_MT, value = variable, 1:2) %>%
      dplyr::select(-meltPA_MT)})
meltPA_MT <-as.data.frame(meltPA_MT) %>% mutate(GRP=as.factor(variable))
# empty data frame for results
results_PA_MT <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("pop1", "pop2", "stat", "temp", "p.value"))
n <- as.numeric(length(pairs))

library("coin")
# loop over pairs
results_PA_MT <-list()
for(p in 1:n){
  pair <- pairs[[p]]$variable
  temp <- meltPA_MT %>%
    dplyr::filter(variable %in% pair) %>%
    mutate(GRP = ordered(variable, levels = pair),
           LOCUS = as.factor(LOCUS)) %>%
    droplevels()
  wilcox <- wilcoxsign_test(value ~ GRP | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,
                            zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)))
  results_PA_MT[[p]] <-df}
results_PA_MT <- bind_rows(results_PA_MT)%>%mutate(test="PA", padj=p.adjust(p.value, method="BH"))

# Private alleles in Shinnecock
#### Wilcoxon test for private alleles
#list of all the different possible pairs
comp <- as.character(unique(meltPA_SHIN$variable))  
pairs <- expand.grid(comp, comp) %>%
  filter(!Var1 == Var2) %>%
  rownames_to_column("PAIR") %>%
  split(.$PAIR) %>%
  purrr::map(function(x){
    x %>%
      dplyr::select(-PAIR) %>%
      gather(key = meltPA_SHIN, value = variable, 1:2) %>%
      dplyr::select(-meltPA_SHIN)})
meltPA_SHIN <-as.data.frame(meltPA_SHIN) %>% mutate(GRP=as.factor(variable))
# empty data frame for results
results_PA_SHIN <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("pop1", "pop2", "stat", "temp", "p.value"))
n <- as.numeric(length(pairs))

library("coin")
# loop over pairs
results_PA_SHIN <-list()
for(p in 1:n){
  pair <- pairs[[p]]$variable
  temp <- meltPA_SHIN %>%
    dplyr::filter(variable %in% pair) %>%
    mutate(GRP = ordered(variable, levels = pair),
           LOCUS = as.factor(LOCUS)) %>%
    droplevels()
  wilcox <- wilcoxsign_test(value ~ GRP | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,
                            zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)))
  results_PA_SHIN[[p]] <-df}
results_PA_SHIN <- bind_rows(results_PA_SHIN)%>%mutate(test="PA", padj=p.adjust(p.value, method="BH"))

##### Simpson's index Shinnecock ####
##### 
llocstats <- loc_stats_SHI ## This loc stats shinnecock 
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
#SIMPSON_IDX
results_list <-list()
for(p in 1:n){
  #tryCatch({
  pair <- pairs[[p]]$GRP
  temp <- llocstats %>%
    dplyr::filter(GRP %in% pair) %>%
    mutate(GRP = ordered(GRP, levels = pair),
           LOCUS = as.factor(LOCUS)) %>%
    droplevels()
  wilcox <- wilcoxsign_test(SIMPSON_IDX ~ GRP | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)),
                   "test" = "SIMPSON_IDX") 
  #if(!is.null(df)){
  results_list[[p]] <- df #}#, error=function(e){})
}
#results_list <-results_list[-which(sapply(results_list, is.null))]
results_simpson_shin<-bind_rows(results_list)%>%mutate(padj=p.adjust(p.value, method="BH"))

#within-bay results
results_mtshi <-bind_rows(results_simpson_shin,results_PA_MT,results_PA_SHIN)
write.csv(results_mtshi, file="./diversity_stats/diversity_output_files/YOY16_bay/wilcox_mtshi.csv")


############################### IR ##################################
# lets try testing for significant differences in mean IR for each bay pair with a wilcox sign test, idk why that's not ok... 
comp <- as.character(unique(rel2$pop))  
pairs <- expand.grid(comp, comp) %>%
  filter(!Var1 == Var2) %>%
  rownames_to_column("PAIR") %>%
  split(.$PAIR) %>%
  purrr::map(function(x){
    x %>%
      dplyr::select(-PAIR)  %>%
      gather(key = rel2, value = name, 1:2) %>%
      dplyr::select(-rel2)})
rel2 <-as.data.frame(rel2) 
# empty data frame for results
results_ir <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("pop1", "pop2", "stat", "temp", "p.value"))
n <- as.numeric(length(pairs))

library("coin")
# loop over pairs
results_ir <-list()
for(p in 1:n){
  pair <- pairs[[p]]$name
  temp <- rel2 %>%
    dplyr::filter(name %in% pair) %>%
    mutate(name = ordered(name, levels = pair),
           Ind = as.factor(Ind)) %>%
    droplevels()
  wilcox <- wilcoxsign_test(IR ~ name,   
                            data = temp,
                            zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)))
  results_ir[[p]]<-df}
  results_ir <- bind_rows(results_ir)
#these p.values are too low. 

#what about an anova followed by ls means? 
#ls means doesn't work but tukey hsd does. 
library("agricolae")
library("lsmeans")
anoIR <-lm(IR~pop, na.action=na.omit, data=rel2)
ano2 <-car::Anova(anoIR)
ano2
df<-df.residual(anoIR)
MSerror<-deviance(anoIR)/df
comparison <- HSD.test(anoIR,c("pop"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/6, group=TRUE)
comparison

## test for significant differences in mean IR for each bay-pair, with t test. Print only the p values. 
jam <-filter(rel2,pop=="Jam")
shin <-filter(rel2,pop=="Shin")
nap <-filter(rel2,pop=="Nap")
mor <-filter(rel2,pop=="Mor")
mt <-filter(rel2,pop=="Mt")

#Nap v. Shin
#var.test(nap$IR,shin$IR) 
t <-t.test(nap$IR,shin$IR,var.equal=TRUE)$p.value
a <-c(t, "Nap","Shin")
#nap v. mor
#var.test(nap$IR,mor$IR)
t <-t.test(nap$IR,mor$IR,var.equal=FALSE)$p.value
b <-c(t, "Nap","Mor")
#jam v. nap
#var.test(jam$IR,nap$IR)
t <-t.test(jam$IR,nap$IR,var.equal=TRUE)$p.value
c <-c(t, "Jam","Nap")
#mor v shin
#var.test(mor$IR,shin$IR)
t <-t.test(mor$IR,shin$IR,var.equal=TRUE)$p.value
d <-c(t, "Mor","Shin")
#jam v shin
#var.test(jam$IR,shin$IR)
t <-t.test(jam$IR,shin$IR,var.equal=TRUE)$p.value
e <-c(t, "Jam","Shin")
#jam v. mor
#var.test(jam$IR,mor$IR)
t <-t.test(jam$IR,mor$IR,var.equal=TRUE)$p.value
f <-c(t, "Jam","Mor")
#jam v. mt
#var.test(jam$IR,mt$IR)
t <-t.test(jam$IR,mt$IR,var.equal=TRUE)$p.value
g <-c(t, "Jam","Mt")
#shin v mt
#var.test(shin$IR,mt$IR)
t <-t.test(shin$IR,mt$IR,var.equal=TRUE)$p.value
h <-c(t, "Shin","Mt")
#mor v mt
#var.test(mor$IR,mt$IR)
t <-t.test(mor$IR,mt$IR,var.equal=TRUE)$p.value
i <-c(t, "Mor","Mt")
#nap v mt
#var.test(nap$IR,mt$IR)
t <-t.test(nap$IR,mt$IR,var.equal=TRUE)$p.value
j <-c(t, "Nap","Mt")

t <- rbind(a,b,c,d,e,f,g,h,i,j)
colnames(t) <- c("pr.t", "bay1","bay2")
t<-as.data.frame(t)

#significant differences between bays. 
#hmmm I think there is a more efficient way to code this. 
# why not a wilcox test? you can definitely still make a heatmap for this.  

pairwise.t.test(rel2$IR, rel2$name) #idk if this is right, but let's automate the t tests. 
t<-mutate(t,pr.t=as.numeric(as.character(pr.t)))%>% mutate(p_value = cut(pr.t, breaks=c(0,0.005, 0.01,0.05,0.1,0.5,1)), significance=ifelse(pr.t < 0.005, "sig","not.sig"))
t<- arrange(t, bay1) %>% mutate(bay1=as.character(bay1), bay2=as.character(bay2)) %>%mutate(pair1 = pmin(bay1,bay2), pair2 =pmax(bay1,bay2)) %>% arrange(pair1)

write.csv(t,"./diversity_stats/diversity_output_files/YOY16_bay/IR_ttests.csv")

######################################## DISTANCE ###############################################################
#### Genetic Distance ####

## Tree's using provesti's distance

#all bays?
set.seed(999)
jpeg(file="./diversity_stats/diversity_figs/YOY16/provesti.jpeg")
wfpopLD %>%
  genind2genpop(pop = ~Bay) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = provesti.dist, missingno="ignore")
dev.off()
#This one makes a lot of sense


#### Isolation by distance.####
library(MASS)
setPop(wfpopLD) <-~Bay
wfpop.geo <-genind2genpop(wfpopLD)
Dgen <-dist.genpop(wfpop.geo,method=5) #Prevosti, but other methods are available. like Edwards nei, etc.

#Create a coordinates list:
coord.dist =matrix(c(41.008611,72.048611,40.791667,72.696111,40.606389,73.876389,40.845278,72.500556,40.999167,72.545278),nrow=5,ncol=2,byrow=TRUE)
rownames(coord.dist) <-c("Nap","Mor","Jam","Shin","Mt")
colnames(coord.dist)<-c("x","y")

#i'm not sure if it's ok to just put in your own euclidean distance matrix like this but it seems to work. 
#previously i thought you had to attach it to the genpop objecte as $other, maybe as xy coordinates? 
#I am not sure if it understands that it corresponds to population, but let's plot anyway? 
#Dgeo <-dist(coord.dist)  # this is a euclidean distance matrix for geographic coordinates but not a shortest distnce by sea measure. 

#Experimental - manually inputting the shortest distance by sea. 
short.dist <-data.frame(site.x=c("Nap","Nap","Nap","Nap","Mor","Mor","Mor","Jam","Jam","Shin"), 
                        site.y=c("Mor","Jam","Shin","Mt","Jam","Shin","Mt","Shin","Mt","Mt"),
                        Distance=c(74.09,195.2,58.16,56.82,109.04,13.28,144.21,132.07,161.5,81.06))


Dgeo2 <- with(short.dist, Distance)
nams <- with(short.dist, unique(c(as.character(site.x), as.character(site.y))))
attributes(Dgeo2) <- with(short.dist, list(Size = length(nams),
                                           Labels = nams,
                                           Diag = FALSE,
                                           Upper = FALSE,
                                           method = "user"))
Dgeo2
class(Dgeo2) <- "dist"
Dgeo2 


ibd2 <-mantel.randtest(Dgen,Dgeo2, nrepet=10000)
ibd2
plot(ibd2)
library("MASS")
dens <- kde2d(Dgeo2,Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
jpeg(file="./diversity_stats/diversity_figs/bay/IBDbysea.jpeg")
plot(Dgeo2, Dgen, xlab="distance (KM)", ylim= c(0.12, 0.35),ylab="genetic distance (euclidean)")
image(dens, col=transp(myPal(500),0.4), add=TRUE)
abline(lm(Dgen~Dgeo2))
dev.off()
#strong looking pattern of IBD. 
library("car")
Anova(lm(Dgen~Dgeo2))
summary(lm(Dgen~Dgeo2))



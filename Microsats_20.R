### Genetic analysis Bays ###
## created 8-19-2020

## NOTES ##
## updated on Sept. 5, 2020, to include all loci but subset as needed. 
## for no half 0 genotypes, use the double0 version of the files 
## population corrected for new cohort assignment
## based on the script "skreportdnaJAN_28_2020MATTITUCK.R"
## Only including analysis that is for the publication. No exploratory analysis. 


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

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")

##### Formating the dataset #####
# We are going to use the doubl0 version. 
wfpop <- read.genalex("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_17_sept20204genalex_doubl0.csv")
wfpop4df <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_17_sept2020_doubl0.csv", header = TRUE) #csv version 

splitStrata(wfpop) <-~Ocean/Bay/Con/Year
setPop(wfpop) <-~Bay

#look at missing data
#Where missing data is greater than 10% for that locus...the locus is not informative for that bay.
info_table(wfpop, plot = TRUE, scaled =FALSE)
#ggsave("rawinfotable20.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()

#Create the LD dataset to check for Linkage Disequilibreum 
#Linkage disequilibreum  in a NONE MISSING dataset- have to do this. 
wfpop2CLEAN <- wfpop %>% missingno("geno", cutoff=0.0)
setPop(wfpop2CLEAN) <-~Bay
wfia.pair <-wfpop2CLEAN %>% clonecorrect(strata= ~Bay) %>% pair.ia(quiet=FALSE)
wfia.pair <- seppop(wfpop2CLEAN) %>% lapply(pair.ia) #by bay!
#ggsave("rawLD20.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()

#create a second dataset that removes one loci in the WF27/WF33 pair 
#also remove WF06 & WF32 due to null alllels, see null_checker script and "microsat record Sept2020.ppt" for details. 
#go back to wfpop for this. 
wfpopLD <-genclone2genind(wfpop)
all_loci <- locNames(wfpopLD)# create vector of all loci
removeloc <- c("WF06")# create vector containing loci to remove
keeploc <- setdiff(all_loci, removeloc)# create list of loci to keep
wfpopLD <- wfpopLD[loc = keeploc]# filter loci in genind object
length(locNames(wfpopLD))# check number of loci in genind obj

#now re-remove the missing data. 
setPop(wfpop2CLEAN) <-~Bay
info_table(wfpopLD, plot = TRUE, scaled =FALSE)
wfpopLD <-wfpopLD %>% missingno("loci", cutoff=0.30) # removes loci where overall missingness > 30%, so WF01
info_table(wfpopLD, plot = TRUE, scaled =FALSE)
wfpopLD <- wfpopLD %>% missingno("geno", cutoff=0.25) # remove samples that are > 25% missing
info_table(wfpopLD, plot = TRUE, scaled =FALSE)
#ggsave("25percutinfotable16LD.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()

#HWE Heatmap#
setPop(wfpopLD) <-~Ocean
hw.test(wfpopLD, B=1000) #permutation based
hw.test(wfpop2CLEAN, B=1000)
hw.test(wfpopLD, B=0) #analytical p value
wfhwe.pop <- seppop(wfpopLD) %>% lapply(hw.test)
write.csv(wfhwe.pop, file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/wfhwepop16.csv")
(wfhwe.mat <- sapply(wfhwe.pop, "[", i = TRUE, j = 3)) # Take the third column with all rows ---> output this for supplementary tables.
wfhw.mc <-sapply(wfhwe.pop, "[", i = TRUE, j = 4) #the PR exact based on the Monte carlo test! ----> p.values on the hw.test
wfhw.mc  # this is just the p values. 
alpha  <- 0.05
newmat <- wfhwe.mat
newmat[newmat > alpha] <- 1 #where the p value on the chi square is greater than 0.05, give it a 1.
# so pink means zero, which means p < 0.05, which means out of HWE. 
levelplot(t(newmat),scales=list(x=list(rot=90)))
#ggsave("HWEtest16.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()

#PopGenReport - Remember to turn this off. 
#setPop(wfpopLD) <-~Bay
#wf.gen <-genclone2genind(wfpopLD) 
#popgenreport(wf.gen,mk.counts=TRUE,mk.locihz = TRUE, mk.fst=TRUE, mk.allele.dist=TRUE, mk.null.all=TRUE,mk.allel.rich = TRUE,mk.differ.stats = TRUE,path.pgr=getwd(),mk.Rcode=TRUE,mk.pdf=TRUE )



### Summary Data ####
setPop(wfpopLD) <-~Bay
toto <-summary(wfpopLD)
barplot(toto$loc.n.all, ylab="Number of alleles",las=2,
        main="Number of alleles per locus")
barplot(toto$Hexp-toto$Hobs, main="Heterozygosity: expected-observed",ylab="Hexp - Hobs",las=2, ylim=c(-0.05, 0.25))
barplot(toto$n.by.pop, main="Sample sizes per population", ylab="Number of genotypes",las=3)
toto
hexhobs <-as.data.frame(toto$Hexp-toto$Hobs) 
hexhobs <- tibble::rownames_to_column(hexhobs,"locus")
names(hexhobs) <-c("locus","difference")

#by bay
hh <- seppop(wfpopLD) %>% lapply(summary)

hhh <-function(hexp, hobs){
  bayhet <-as.data.frame(hexp)%>% tibble::rownames_to_column("locus")
  bayhobs <-as.data.frame(hobs)%>% tibble::rownames_to_column("locus")
  bayhh <-left_join(bayhet, bayhobs, by=c("locus")) %>% mutate(hexhobs=hexp-hobs)
  bayhh
}
naph <-hhh(hh$Nap$Hexp, hh$Nap$Hobs) %>% mutate(Bay="Nap")
morh <-hhh(hh$Mor$Hexp, hh$Mor$Hobs) %>% mutate(Bay="Mor")
Mth <-hhh(hh$Mt$Hexp, hh$Mt$Hobs) %>% mutate(Bay="Mt")
jamh <-hhh(hh$Jam$Hexp, hh$Jam$Hobs) %>% mutate(Bay="Jam")
shinh <-hhh(hh$Shin$Hexp, hh$Shin$Hobs) %>% mutate(Bay="Shin")
h4 <- bind_rows(naph, morh, Mth, jamh, shinh)

#heatmap of observed minus expected. 
ggplot(h4, aes(x = locus, y = Bay, fill=hexhobs)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(hexhobs, 3)), color="grey") +
  scale_fill_viridis(option="plasma") +
  coord_fixed(ratio = 1) +
  ylab("")+
  #theme_standard() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave('hexhobsheat17.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)


###Summary stats from popgenreport. 


##### Shannon's summary stats ########
#https://gist.github.com/sjoleary/cdc32efbfd50c96eef446ebb7c2f2387

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

setPop(wfpopLD)<-~Bay/Con
gen_con <-seppop(wfpopLD)

setPop(wfpopLD)<-~Bay/Con/Year
gen_bayconyear <-seppop(wfpopLD)

setPop(wfpopLD)<-~Year
gen_year <- seppop(wfpopLD)

setPop(wfpopLD)<-~Bay/Year
gen_bayyear <- seppop(wfpopLD)

gen_grp <-c(gen_oce,gen_bay,gen_con,gen_bayconyear,gen_year)
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

loc_stats[is.na(loc_stats)] <- NA
#write.csv(loc_stats,file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/loc_stats17.csv")

#create different loc stats groups for testing. 
loc_stats_bay <-filter(loc_stats, GRP %in% c("Nap","Mor","Jam","Shin","Mt"))

#bay melt
meltlocstats_bay <-pivot_longer(loc_stats_bay,cols = c("N_ALLELES", "SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis","SHANNON_IDX","STODD_TAYLOR_IDX"),
                                names_to="variable", values_to="value")


meltlocstats2 <-filter(meltlocstats, GRP %in% c("Nap","Mor","Jam","Shin","Mt","ALL")) 
#this is essentially by bay but also includes the "ALL" Group


#we have some reason to distrust the fact that ht and hs are the same after shannons code but not globally. so we will also do this another way.
#if you look at this, basic stats does not work at the subpopulation level and this could be why. 
g2h_all <-genind2hierfstat(wfpopLD)
g2hall <-basic.stats(g2h_all)
global <-as.data.frame(g2hall$perloc) %>% tibble::rownames_to_column("locus") %>% mutate(Bay="ALL")
#write.csv(global,file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/basic_stats17.csv" )

setPop(wfpopLD) <-~Bay
g2h_bay <-genind2hierfstat(wfpopLD, pop=wfpopLD@pop)

wfshin <- popsub(wfpopLD, sublist=c("Shin"))  #%>% missingno("geno", cutoff=0.10)
wfshinh <- genind2hierfstat(wfshin) 
wfshi<-basic.stats(wfshinh)
#wfshi <-diff_stats(wfshin)
wfnap <- popsub(wfpopLD, sublist=c("Nap"))
wfmor <- popsub(wfpopLD, sublist=c("mor"))

##### Box and whisker plots for diversity stats #####
drabcolors <-c("#d0d1e6","#a6bddb", "#67a9cf", "#1c9099", "#016450")

#Nei's gene diversity
meltlocstats_bay %>%
  filter(variable == "Ht") %>%
  ggplot(aes(x=fct_rev(GRP),y=value),fill=GRP)+
  geom_boxplot(aes(fill=GRP))+ 
  scale_fill_manual(name = "Bay",values = drabcolors)+
  coord_flip()+ 
  #ylim(0.7,1.0)+
  xlab(' ')+ylab("Nei's Gene Diversity")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
#ggsave('nei_bay17.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

#Inbreeding coefficient
meltlocstats_bay %>%
  filter(variable == "Fis") %>%
  ggplot(aes(x=fct_rev(GRP),y=value),fill=GRP)+
  geom_boxplot(aes(fill=GRP))+ 
  scale_fill_manual(name = "Bay",values = drabcolors)+
  coord_flip()+ 
  #ylim(0.7,1.0)+
  xlab(' ')+ylab("Fis")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
#ggsave('FIS_bay17.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

#Evenness
meltlocstats_bay %>%
  filter(variable == "EVENNESS") %>%
  ggplot(aes(x=fct_rev(GRP),y=value),fill=GRP)+
  geom_boxplot(aes(fill=GRP))+ 
  scale_fill_manual(name = "Bay",values = drabcolors)+
  coord_flip()+ 
  #ylim(0.7,1.0)+
  xlab(' ')+ylab("Evenness")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('Evenness_bay17.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

#Shannon index
meltlocstats_bay %>%
  filter(variable == "SHANNON_IDX") %>%
  ggplot(aes(x=fct_rev(GRP),y=value),fill=GRP)+
  geom_boxplot(aes(fill=GRP))+ 
  scale_fill_manual(name = "Bay",values = drabcolors)+
  coord_flip()+ 
  #ylim(0.7,1.0)+
  xlab(' ')+ylab("Shannon Index")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('Shannon_idxbay17.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

#Expected heterozygosity
meltlocstats_bay %>%
  filter(variable == "Hs") %>%
  ggplot(aes(x=fct_rev(GRP),y=value),fill=GRP)+
  geom_boxplot(aes(fill=GRP))+ 
  scale_fill_manual(name = "Bay",values = drabcolors)+
  coord_flip()+ 
  #ylim(0.7,1.0)+
  xlab(' ')+ylab("Expected Heterozygosity")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('Hsbay17.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

#Observed heterozygosity
meltlocstats_bay %>%
  filter(variable == "Ho") %>%
  ggplot(aes(x=fct_rev(GRP),y=value),fill=GRP)+
  geom_boxplot(aes(fill=GRP))+ 
  scale_fill_manual(name = "Bay",values = drabcolors)+
  coord_flip()+ 
  #ylim(0.7,1.0)+
  xlab(' ')+ylab("Observed Heterozygosity")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('Hobay17.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)


######## CALCULATE RAREFIED ALLELIC RICHNESS ----
# differences in sample size can bias the number of alleles sampled in a population
# calculate allelic richness corrected for sample size using rarefaction

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
#df <-allelic.richness(wfpopLD, diploid=TRUE)

# I am not sure how to tell which one is which in order to label the columns, but I assume they're in the same order as WFPOPLD?
df <- as.data.frame(df$Ar) %>%
  rownames_to_column("LOCUS") %>%
  dplyr::rename(Nap = V1,
                Mor = V2,
                Jam = V3,
                Shin = V4,
                Mt = V5)

ar <- left_join(ar, df)

#write_delim(ar, "results/rarefied.allelecount", delim = "\t")

#visualize results as boxplot. 
library("forcats")
library("cowplot")

#for rarified alleles, you have to make a new Allelic Richness analysis for each group you are comparing. 
meltar2 <- pivot_longer(ar, cols=c("Mt","Shin","Nap","Mor","Jam","ALL"),names_to="variable", values_to="value")
meltar3 <-filter(meltar2, variable != "ALL") #it doenst like it when you exclude the "all" category... 

meltar3 %>%
  #filter(variable !="ALL")%>%
  mutate(GRP = as.factor(variable)) %>%
  ggplot(aes(x=fct_rev(GRP),y=value),fill=GRP)+
  geom_boxplot(aes(fill=GRP))+ 
  scale_fill_manual(name = "Bay",values = drabcolors)+
  coord_flip()+ 
  #ylim(0.7,1.0)+
  xlab(' ')+ylab("Rareified allele count")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('rareifiedallelesLD_bay17.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

##Private alleles
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

meltpA %>%
  mutate(GRP = as.factor(variable)) %>%
  ggplot(aes(x=fct_rev(GRP),y=value),fill=GRP)+
  geom_boxplot(aes(fill=GRP))+ 
  scale_fill_manual(name = "Bay",values = drabcolors)+
  coord_flip()+ 
  #ylim(0.7,1.0)+
  xlab(' ')+ylab("Private alleles")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('privateallelesLD_bay17.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)



# TEST FOR SIGNIFICANT DIFFERENCES
# need to use Friedman's test for global test and Wilcoxon signed rank for pairwise tests 
# to test symmetry of numeric repeated measurements (stastic per locus) in block design
# Example using gene diversity (expected heterozygosity by estuary)

#Shannon's way - go to the code she does it a little differently.  
# importantly she removes rows with NA from the data set, but we don't have any. 

#the different levels we have to work with are:
    #loc_stats_bay <- bays
    #loc_stats_BYC <- bay year con
    #loc_stats_MT <-bay con for mattituck only, no group 5 (unidentified adults)
    #loc_stats_shin <- bay year con shin
    #loc_stats_16 <- 2016 YOY at bay level. 

#set level to test for friedmans
fdata <- loc_stats_bay


###Friedman's test for Global heterogeneity
friedman.test(Hs~GRP | LOCUS, data= fdata) #Hs
friedman.test(EVENNESS~GRP | LOCUS, data= fdata) #evenness
friedman.test(Ht~GRP | LOCUS, data= fdata) #Ht
friedman.test(Fis~GRP | LOCUS, data= fdata) #Fis
friedman.test(SHANNON_IDX~GRP | LOCUS, data= fdata) 
friedman.test(value~variable | LOCUS, data= meltar3) #rareified alleles #come back to this one.... 



#You could automate the test result extraction similar to how you did with the lm summary,,, but not right now. 

#######Wilcoxon tests #####
#tell us which groups are significantly different. 
# you have to set which level you are comparing. 
  #the different levels we have to work with are:
    #loc_stats_bay <- bays
    #loc_stats_BYC <- bay year con
    #loc_stats_MT <-bay con for mattituck only, no group 5 (unidentified adults)
    #loc_stats_shin <- bay year con shin
    #loc_stats_16 <- 2016 YOY at bay level. 

llocstats <- loc_stats_bay ##CHANGE HERE
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
  wilcox <- wilcoxsign_test(SHANNON_IDX ~ GRP | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,
                            zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)),
                   "test" = "Shannon")  ### change this one. 
  results <- bind_rows(results, df)}

#now change test and re_run above. 
results_nei <-results
results_fis <-results
results_even <-results
results_shannon <-results

results <-bind_rows(results_nei, results_fis, results_even, results_shannon)
results <- results %>% dplyr::select(-temp)
results <-mutate(results, bonferroni=0.05/10) # number of pairwise comparisons ####CHECK THIS
results <-mutate(results, significance = ifelse(p.value>=bonferroni,"not significant","significant"))

write.csv(results,file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/wilcox.csv" )

#heatmap of results (new way)
results <- mutate(results, starsig=ifelse(significance=="significant",round(p.value,4),NA),test.statistic=abs(stat))
results<-mutate(results, p_value = cut(p.value, breaks=c(0,0.005, 0.01,0.05,0.1,0.5,1)))
cols <- c("(0,0.005]"="#034e7b", "(0.005,0.01]" = "#045a8d", "(0.01,0.05]" = "#2b8cbe", "(0.05,0.1]" = "#74a9cf", "(0.1,0.5]"  = "#a6bddb", "(0.5,1]"="#d0d1e6")

#NEI Heatmap. 
#delete duplicate pairs to form the half grid of the heatmap. 
results_nei <-filter(results, test=="Ht") %>%
  arrange(test.statistic) %>% unique() %>% 
  mutate(pair1 = pmin(pop1,pop2), pair2 =pmax(pop1,pop2)) %>% arrange(pair1)
toDelete <- seq(1, nrow(results_nei), 2)
results_nei <- results_nei[ toDelete ,]

# Filled by test statistic: NEI - shows graphical options. 
results_nei %>%
  ggplot(aes(x = pair1, y = pair2))+
  #geom_tile(aes(fill=test.statistic), show.legend = TRUE)+ # Filled by test statistic
  #geom_tile(aes(fill=p.value), show.legend = TRUE)+ #filled  by p.value (gradient)
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  #scale_fill_gradient(low ="#023858" , high = "#a6bddb", space = "Lab", na.value = "white", guide = "colourbar", aesthetics = "fill")+ #gradient fill
  scale_fill_manual(values=cols)+
  geom_text(aes(label = round(p.value,3),color=significance), size=5)+
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
#ggsave('neiBaytile17.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

#FIS heatmap
#delete duplicate pairs to form the half grid of the heatmap. 
results_Fis <-filter(results, test=="Fis") %>%
  arrange(test.statistic) %>% unique() %>% 
  mutate(pair1 = pmin(pop1,pop2), pair2 =pmax(pop1,pop2)) %>% arrange(pair1)
toDelete <- seq(1, nrow(results_Fis), 2)
results_Fis <- results_Fis[ toDelete ,]

# Filled by test statistic: FIS
results_Fis %>%
  ggplot(aes(x = pair1, y = pair2))+
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  scale_fill_manual(values=cols)+
  geom_text(aes(label = round(p.value,3),color=significance), size=5)+
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
#ggsave('FisBaytile17.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

#EVENNESS heatmap
#delete duplicate pairs to form the half grid of the heatmap. 
results_even <-filter(results, test=="Evenness") %>%
  arrange(test.statistic) %>% unique() %>% 
  mutate(pair1 = pmin(pop1,pop2), pair2 =pmax(pop1,pop2)) %>% arrange(pair1)
toDelete <- seq(1, nrow(results_even), 2)
results_even <- results_even[ toDelete ,]

# Filled by test statistic: evenness
results_even %>%
  ggplot(aes(x = pair1, y = pair2))+
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  scale_fill_manual(values=cols)+
  geom_text(aes(label = round(p.value,3),color=significance), size=5)+
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
#ggsave('evennessBaytile17.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

#shannons heatmap
#delete duplicate pairs to form the half grid of the heatmap. 
results_shannon <-filter(results, test=="Shannon") %>%
  arrange(test.statistic) %>% unique() %>% 
  mutate(pair1 = pmin(pop1,pop2), pair2 =pmax(pop1,pop2)) %>% arrange(pair1)
toDelete <- seq(1, nrow(results_shannon), 2)
results_shannon <- results_shannon[ toDelete ,]

# Filled by test statistic: evenness
results_shannon %>%
  ggplot(aes(x = pair1, y = pair2))+
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  scale_fill_manual(values=cols)+
  geom_text(aes(label = round(p.value,3),color=significance), size=5)+
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
#ggsave('shannonBaytile17.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)


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
  results_ar <- bind_rows(results_ar, df)}



results_ar <-mutate(results_ar, bonferroni=0.05/10)
results_ar <-mutate(results_ar, significance = ifelse(p.value>=bonferroni,"not significant","significant")) %>% mutate(starsig=ifelse(significance=="significant",round(p.value,4),NA),test.statistic=abs(stat))
results_ar<-mutate(results_ar, p_value = cut(p.value, breaks=c(0,0.005, 0.01,0.05,0.1,0.5,1)))
cols <- c("(0,0.005]"="#034e7b", "(0.005,0.01]" = "#045a8d", "(0.01,0.05]" = "#2b8cbe", "(0.05,0.1]" = "#74a9cf", "(0.1,0.5]"  = "#a6bddb", "(0.5,1]"="#d0d1e6")

results_ar <-arrange(results_ar, test.statistic) %>% unique() %>% 
  mutate(pair1 = pmin(pop1,pop2), pair2 =pmax(pop1,pop2)) %>% arrange(pair1, pair2)
#toDelete <- seq(1, nrow(results_ar), 2)
#results_ar <- results_ar[ toDelete ,]

# plot heatmap of results ====
results_ar %>%
  ggplot(aes(x = pair1, y = pair2))+
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  scale_fill_manual(values=cols)+
  geom_text(aes(label = round(p.value,3),color=significance), size=5)+
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
#ggsave('rariefied_allelesBaytile17.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)



##### FST ######

#Weir and Cockheram FST global & pairwise. - probably your best bet. 

wfpopLD <-wfpop17

setPop(wfpopLD) <-~Bay
wf.g2 <-genind2gtypes(wfpopLD)
statFst(wf.g2)
#popStruct <- popStructTest(wf.g2, stats = c(statFst, statFstPrime), nrep = 1000, quietly = TRUE)
popStruct <- popStructTest(wf.g2, nrep = 1000, quietly = FALSE)
popStruct

#FST heatmap
bayFST <-as.data.frame(popStruct$pairwise$result) %>% dplyr::select(strata.1, strata.2, Fst, Fst.p.val)

# Filled by test statistic: NEI - shows graphical options. 
bayFST %>%
  ggplot(aes(x = strata.1, y = strata.2))+
  geom_tile(aes(fill=Fst), show.legend = TRUE)+ # Filled by test statistic
  #geom_tile(aes(fill=p.value), show.legend = TRUE)+ #filled  by p.value (gradient)
  #geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  #scale_fill_gradient(low ="#023858" , high = "#a6bddb", space = "Lab", na.value = "white", guide = "colourbar", aesthetics = "fill")+ #gradient fill
  scale_fill_viridis_c()+
  #scale_fill_manual(values=cols)+
  geom_text(aes(label = round(Fst.p.val,3)), color="white", size=5)+
  #scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
ggsave('FSTBaytile17.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)

#see previous scripts for other kinds of FST you can calculate


#### Inbreeding - Internal Relatedness ####
library("Rhh")
setPop(wfpopLD) <-~Bay
#wfp <-genind2df(wfpopLD)
wfp <-genind2df(wfpopLD, usepop = TRUE, oneColPerAll = TRUE) %>% tibble::rownames_to_column("Ind")
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

ddply(rel2, ~pop, summarize, meanir=mean(IR))

#barplot as before
ggplot(rel2, aes(x=pop,y=IR))+
  #ggplot(aes(x=fct_inorder(GRP),y=value))+
  geom_boxplot(fill="lightgray")+ 
  coord_flip()+ 
  #ylim(0.7,1.0)+
  xlab(' ')+ylab("IR")+theme_cowplot()+guides(fill = FALSE, colour = FALSE) 
#ggsave('IRld.png', width = 7, height = 7)

rel2 <- arrange(rel2,pop)
drabcolors2 <-c("#d0d1e6","#a6bddb", "#67a9cf", "#1c9099", "#016450")
rel2 <-mutate(rel2,name=fct_relevel(pop,"Jam","Mor","Mt","Nap","Shin"))
rel2 %>%
  ggplot(aes(x=fct_rev(name),y=IR),fill=name)+
  geom_boxplot(aes(fill=name))+ 
  scale_fill_manual(name = "Bay",values = drabcolors2)+
  coord_flip()+ 
  #ylim(0.7,1.0)+
  xlab(' ')+ylab("Internal Relatedness")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
#ggsave('rel_bay17.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

# lets try testing for significant differences in mean IR for each bay pair with a wilcox sign test, idk why that's not ok... 
comp <- as.character(unique(rel2$name))  
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
  results_ir <- bind_rows(results_ir, df)}
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
cols <- c("(0,0.005]"="#034e7b", "(0.005,0.01]" = "#045a8d", "(0.01,0.05]" = "#2b8cbe", "(0.05,0.1]" = "#74a9cf", "(0.1,0.5]"  = "#a6bddb", "(0.5,1]"="#d0d1e6")
t<- arrange(t, bay1) %>% mutate(bay1=as.character(bay1), bay2=as.character(bay2)) %>%mutate(pair1 = pmin(bay1,bay2), pair2 =pmax(bay1,bay2)) %>% arrange(pair1)

t %>%
  ggplot(aes(x = pair1, y = pair2))+
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ # Filled by test statistic
  #geom_tile(aes(fill=p.value), show.legend = TRUE)+ #filled  by p.value (gradient)
  #geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  #scale_fill_gradient(low ="#023858" , high = "#a6bddb", space = "Lab", na.value = "white", guide = "colourbar", aesthetics = "fill")+ #gradient fill
  scale_fill_manual(values=cols)+
  geom_text(aes(label = round(pr.t,3),color=significance),size=5)+
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
#ggsave('IRBaytile17.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

#### Genetic Distance ####

## Tree's using provesti's distance

#all bays?
set.seed(999)
wfpopLD %>%
  genind2genpop(pop = ~Bay) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = provesti.dist, missingno="ignore")
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


ibd2 <-mantel.randtest(Dgen,Dgeo2, nrepet=1000)
ibd2
plot(ibd2)
library("MASS")
dens <- kde2d(Dgeo2,Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo2, Dgen, xlab="distance (KM)", ylim= c(0.12, 0.35),ylab="genetic distance (euclidean)")
image(dens, col=transp(myPal(500),0.4), add=TRUE)
abline(lm(Dgen~Dgeo2))
#strong looking pattern of IBD. 
library("car")
Anova(lm(Dgen~Dgeo2))
summary(lm(Dgen~Dgeo2))




##### AMOVA #####
#kevin said better to do AMOVA in Arlequin
setPop(wfpopLD) <- ~Bay
gapless <-missingno(wfpopLD) # only works with gapless data. 
strat.amo <-poppr.amova(gapless,~Bay)
strat.amo

#randomization test
library("ade4")
set.seed(1999)
stratsig  <- ade4::randtest(strat.amo, nrepet = 999)
plot(stratsig)
stratsig

#another AMOVA, this one from pegas. 
setPop(wfpopLD) <- ~Bay #idk how to change it so it doesn't care about Con/Year
bov_dist  <- dist(wfpopLD) # Euclidean distance
bov_stra  <- strata(wfpopLD)
bov_amova <- pegas::amova(bov_dist ~ Bay, data = bov_stra, nperm = 100)
bov_amova
#no statistical difference between bays i think. 

#yet another way to do an amova
library("vegan")
set.seed(20151219)
res <- adonis(bov_dist ~ Bay, data = bov_stra, permutations = 99)
res2<- adonis(bov_dist ~ Bay/Con/Year, data = bov_stra, permutations = 99)

#let's look at different bays, but 2016 only!
setPop(wfpop2) <- ~Bay/Year
wfpop2016 <- popsub(wfpop2, blacklist=c("Shin_1_2017", "Shin_2_2017"))
bov_stra2016  <- strata(wfpop2016)
bov_dist2016  <- dist(wfpop2016)
res3<- adonis(bov_dist2016 ~ Bay, data = bov_stra2016, permutations = 99)

#with gapless datasets
gapless <- missingno(wfpop2, type="geno", cutoff=0, freq=TRUE)
setPop(gapless) <- ~Bay/Year
wfpop2016 <- popsub(gapless, blacklist=c("Shin_1_2017", "Shin_2_2017"))
bov_stra2016  <- strata(wfpop2016)
bov_dist2016  <- dist(wfpop2016)
res3<- adonis(bov_dist2016 ~ Bay, data = bov_stra2016, permutations = 99)

#let's look within Shinnecock Bay
setPop(wfpop2) <- ~Bay/Con/Year
wfpopShin <- popsub(wfpop2, blacklist=c("Nap_6_2016", "Mor_6_2016", "Jam_6_2016"))
bov_straS  <- strata(wfpopShin)
bov_distS  <- dist(wfpopShin)
resS<- adonis(bov_distS ~ Con/Year, data = bov_straS, permutations = 99)


###### some more stats! #####
#private alleles
private_alleles(wfpop2)
privateAlleles(wf.g2)

#locus table
locus_table(wfpopLD)
wflt.pair <- seppop(wfpop2) %>% lapply(locus_table) #by bay!

#diversity boot
mlgtab <-mlg.table(wfpopLD)
divboot <- diversity_boot(mlgtab, n=1000, H=TRUE, G=TRUE, lambda=TRUE, E5=TRUE)
dci <-diversity_ci(wfpop2, n=10000, ci=95, plot=TRUE, center=TRUE, rarefy=FALSE, raw=TRUE)

#with rarefaction

dci <-diversity_ci(wfpop2, n=10000, ci=95, plot=FALSE, center=TRUE, rarefy=TRUE, raw=TRUE)


#other stats
library("mmod")
otherstats<-diff_stats(wfpopLD)
#plot it
per.locus <- melt(otherstats$per.locus, varnames = c("Locus", "Statistic"))
stats     <- c("Hs", "Ht", "Gst", "Gprime_st", "D", "D")
glob      <- data.frame(Statistic = stats, value = otherstats$global)
ggplot(per.locus, aes(x = Statistic, y = value)) +
  geom_boxplot() +
  geom_point() +
  geom_point(size = rel(3), color = "red", data = glob) +
  ggtitle("Estimates of population differentiation")

#bootstrap estimates
set.seed(20151219) 
bs_reps <- chao_bootstrap(wfpop2, nreps = 100) #100 for now. 
summarise_bootstrap(bs_reps, Gst_Nei) # change the function you want to use  (eg, D_Jost, Gst_Hedrick or Gst_Nei) 

#rarefaction to standardize diversity estimates
# i can't find a function to do this, but as I understand it, I should subsample randomly from the existing shinnecock basically. 

#Microsats_ bayyearcon

### Genetic analysis Mattituck and Shinnecock separately ###
## created 9-8-2020

## NOTES ##
## for no half 0 genotypes, use the double0 version of the files 
## population corrected for new cohort assignment
## based on the script "skreportdnaJAN_28_2020MATTITUCK.R" and "Microsats_20.R"
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

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")

##### Formating the dataset #####
# We are going to use the doubl0 version. 
wfpop <- read.genalex("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_20_sept20204genalex_doubl0.csv")
wfpop4df <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_20_sept2020_doubl0.csv", header = TRUE) #csv version 


splitStrata(wfpop) <-~Ocean/Bay/Con/Year
setPop(wfpop) <-~Bay

#remove the same loci as before. 
wfpopLD <-genclone2genind(wfpop)
all_loci <- locNames(wfpopLD)# create vector of all loci
removeloc <- c("WF27","WF06","WF32")# create vector containing loci to remove
keeploc <- setdiff(all_loci, removeloc)# create list of loci to keep
wfpopLD <- wfpopLD[loc = keeploc]# filter loci in genind object
length(locNames(wfpopLD))# check number of loci in genind obj
# and remove missing data
setPop(wfpop2CLEAN) <-~Bay/Year/Con
info_table(wfpopLD, plot = TRUE, scaled =FALSE)
wfpopLD <-wfpopLD %>% missingno("loci", cutoff=0.30) # removes loci where overall missingness > 30%, so WF01
info_table(wfpopLD, plot = TRUE, scaled =FALSE)
wfpopLD <- wfpopLD %>% missingno("geno", cutoff=0.25) # remove samples that are > 25% missing
info_table(wfpopLD, plot = TRUE, scaled =FALSE)

#create the mattituck and shinnecock datasets
setPop(wfpopLD) <-~Bay
wf.shin <-popsub(wfpopLD, sublist=c("Shin"))
wf.mt  <-popsub(wfpopLD, sublist=c("Mt"))
setPop(wf.shin) <-~Bay/Con/Year
setPop(wf.mt) <-~Bay/Con
wf.mt <-popsub(wf.mt, blacklist=c("Mt_5")) # we don't need the unidentified individuals. 
# an important decision here... do we want to combine the early cohorts from 2015 with those from 2016 or not?? 


#create a combined dataset -- This die not work come back to this. 
setPop(wf.mt) <- ~Bay/Con/Year
wf.mt2 <-genind2df(wf.mt, usepop=TRUE, oneColPerAll = TRUE)
wf.mt2 <-mutate(wf.mt2, newpop = ifelse(pop %in% c("Mt_3_2015","Mt_3_2016"),"Mt_3_both",ifelse(pop %in% c("Mt_4_2015","Mt_4_2016"),"Mt_4_both",as.character(pop)))) 
wf.mt2 <-wf.mt2[,c(34,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33)]
wf.mt2 <-dplyr::select(wf.mt2, -pop) %>% dplyr::rename(pop=newpop) %>%add_rownames(var="notind") ###THE IND NAME IS NOT RETAINED!!!!!
#couldn't get df2genind to work. let's write it as a csv and reimport it? 
write.csv(wf.mt2, file="mtnewpop.csv")
wf.mt3 <- read.genalex("/Users//tdolan/Documents//R-Github//WFmicrosats/mtnewpop4genalex.csv")
splitStrata(wf.mt3) <-~Bay/Con/Year

######HWE Heatmap######
#Shinnecock
shihwe.pop <- seppop(wf.shin) %>% lapply(hw.test)
write.csv(shihwe.pop, file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/wfhwepop16shin.csv")
(wfhwe.mat <- sapply(shihwe.pop, "[", i = TRUE, j = 3)) # Take the third column with all rows ---> output this for supplementary tables.
wfhw.mc <-sapply(shihwe.pop, "[", i = TRUE, j = 4) #the PR exact based on the Monte carlo test! ----> p.values on the hw.test
wfhw.mc  # this is just the p values. 
alpha  <- 0.05
newmat <- wfhwe.mat
newmat[newmat > alpha] <- 1 #where the p value on the chi square is greater than 0.05, give it a 1.
# so pink means zero, which means p < 0.05, which means out of HWE. 
levelplot(t(newmat),scales=list(x=list(rot=90)))
ggsave("HWEtest16Shin.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
dev.off()

#Mattituck
mthwe.pop <- seppop(wf.mt) %>% lapply(hw.test)
write.csv(mthwe.pop, file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/wfhwepop16mt.csv")
(wfhwe.mat <- sapply(mthwe.pop, "[", i = TRUE, j = 3)) # Take the third column with all rows ---> output this for supplementary tables.
wfhw.mc <-sapply(mthwe.pop, "[", i = TRUE, j = 4) #the PR exact based on the Monte carlo test! ----> p.values on the hw.test
wfhw.mc  # this is just the p values. 
alpha  <- 0.05
newmat <- wfhwe.mat
newmat[newmat > alpha] <- 1 #where the p value on the chi square is greater than 0.05, give it a 1.
# so pink means zero, which means p < 0.05, which means out of HWE. 
levelplot(t(newmat),scales=list(x=list(rot=90)))
ggsave("HWEtest16mt.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
dev.off()


##### Shannon's summary stats ########
#https://gist.github.com/sjoleary/cdc32efbfd50c96eef446ebb7c2f2387

library(tidyverse)
library(ggplot2)
library(adegenet)
library(hierfstat)

#dataframe with sample information
setPop(wfpopLD) <- ~Bay/Con/Year
df <-genind2df(wfpopLD, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
#reorder so that the first column has to be the sample name. 
df <-df[,c(34,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33)]
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
  #select(-Hexp) %>%
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

# combine into single data frame ====
loc_stats_2 <- plyr::ldply(loc_stats_2, data.frame) %>%
  dplyr::rename(GRP = `.id`)

loc_stats <- left_join(loc_stats, loc_stats_2) %>%
  dplyr::select(GRP, LOCUS, N_ALLELES, EVENNESS, Ho, Hs, Ht, Fis, SHANNON_IDX, SIMPSON_IDX, STODD_TAYLOR_IDX) %>%
  filter(LOCUS != "mean")

loc_stats[is.na(loc_stats)] <- NA
#write.csv(loc_stats,file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/loc_stats16.csv")

#create different loc stats groups for testing. 
loc_stats_MT <-filter(loc_stats, GRP %in% c("Mt_2", "Mt_1", "Mt_3", "Mt_4")) #does not include MT 5 because we don't know what those individuals are. 
loc_stats_shin <-filter(loc_stats, GRP %in% c("Shin_1_2016", "Shin_2_2016","Shin_1_2017","Shin_2_2017"))
#MT melt
meltlocstats_MT <-pivot_longer(loc_stats_MT,cols = c("N_ALLELES", "SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis","SHANNON_IDX","STODD_TAYLOR_IDX"),
                                names_to="variable", values_to="value")
#SHIN melt
meltlocstats_shin <-pivot_longer(loc_stats_shin,cols = c("N_ALLELES", "SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis","SHANNON_IDX","STODD_TAYLOR_IDX"),
                               names_to="variable", values_to="value")

##### Box and whisker plots for diversity stats #####
drabcolors <-c("#d0d1e6","#a6bddb", "#67a9cf", "#1c9099", "#016450")

meltlocstats_MT %>%
  filter(variable == "Ht") %>%
  ggplot(aes(x=fct_rev(GRP),y=value))+
  geom_boxplot(aes(fill=GRP),fill= "#67a9cf")+ coord_flip()+ 
  xlab(' ')+ylab("Nei's Gene Diversity")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('nei_MT.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

meltlocstats_shin %>%
  filter(variable == "Ht") %>%
  ggplot(aes(x=fct_rev(GRP),y=value))+
  geom_boxplot(aes(fill=GRP),fill= "#016450")+ coord_flip()+ 
  xlab(' ')+ylab("Nei's Gene Diversity")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('nei_shin.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

meltlocstats_MT %>%
  filter(variable == "Fis") %>%
  ggplot(aes(x=fct_rev(GRP),y=value))+
  geom_boxplot(aes(fill=GRP),fill= "#67a9cf")+ coord_flip()+xlab(' ')+ylab("Fis")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('FIS_MT.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

meltlocstats_shin %>%
  filter(variable == "Fis") %>%
  ggplot(aes(x=fct_rev(GRP),y=value))+
  geom_boxplot(aes(fill=GRP),fill= "#016450")+ coord_flip()+xlab(' ')+ylab("Fis")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('FIS_shin.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

meltlocstats_MT %>%
  filter(variable == "EVENNESS") %>%
  ggplot(aes(x=fct_rev(GRP),y=value))+
  geom_boxplot(aes(fill=GRP),fill= "#67a9cf")+ coord_flip()+xlab(' ')+ylab("Evenness")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('Even_MT.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

meltlocstats_shin %>%
  filter(variable == "EVENNESS") %>%
  ggplot(aes(x=fct_rev(GRP),y=value))+
  geom_boxplot(aes(fill=GRP),fill= "#016450")+ coord_flip()+xlab(' ')+ylab("Evenness")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('Even_shin.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

meltlocstats_MT %>%
  filter(variable == "SHANNON_IDX") %>%
  ggplot(aes(x=fct_rev(GRP),y=value))+
  geom_boxplot(aes(fill=GRP),fill= "#67a9cf")+ coord_flip()+xlab(' ')+ylab("Shannon Diversity Index")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('shannon_MT.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

meltlocstats_shin %>%
  filter(variable == "SHANNON_IDX") %>%
  ggplot(aes(x=fct_rev(GRP),y=value))+
  geom_boxplot(aes(fill=GRP),fill= "#016450")+ coord_flip()+xlab(' ')+ylab("Shannon Diversity Index")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('shannon_shin.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)





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

# By Subpopulation  ===
# by region ====
datmt <- hierfstat:::.genind2hierfstat(wf.mt, pop = wf.mt@pop)
dfmt <- allelic.richness(datmt,diploid = TRUE)

# I am not sure how to tell which one is which in order to label the columns, but I assume they're in the same order as WFPOPLD?
dfmt <- as.data.frame(dfmt$Ar) %>%
  rownames_to_column("LOCUS") %>%
  dplyr::rename(Mt_2 = V1,
                Mt_1 = V2,
                Mt_3 = V3,
                Mt_4 = V4)
ar <- left_join(ar, dfmt)

#repeat for shin
datshi <- hierfstat:::.genind2hierfstat(wf.shin, pop = wf.shin@pop)
dfshi <- allelic.richness(datshi,diploid = TRUE)
dfshi <- as.data.frame(dfshi$Ar) %>%
  rownames_to_column("LOCUS") %>%
  dplyr::rename(Shin_1_2016 = V1,
                Shin_2_2016 = V2,
                Shin_1_2017 = V3,
                Shin_2_2017 = V4)
ar <- left_join(ar, dfshi)
#write_delim(ar, "results/rarefied.allelecount", delim = "\t")

#visualize results as boxplot. 
library("forcats")
library("cowplot")

#for rarified alleles, you have to make a new Allelic Richness analysis for each group you are comparing. 
meltar2 <- pivot_longer(ar, cols=c("Mt_2","Mt_1","Mt_3","Mt_4","Shin_1_2016","Shin_2_2016","Shin_1_2017","Shin_2_2017","ALL"),names_to="variable", values_to="value")
meltar3 <-filter(meltar2, variable != "ALL") #it doenst like it when you exclude the "all" category... 

#Mattituck rareified alleles
meltar3 %>%
  filter(variable %in%c("Mt_2","Mt_1","Mt_3","Mt_4"))%>%
  mutate(GRP = as.factor(variable)) %>%
  ggplot(aes(x=fct_rev(GRP),y=value))+
  geom_boxplot(aes(fill=GRP),fill= "#67a9cf")+ coord_flip()+ 
  xlab(' ')+ylab("Rareified allele count")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('rareifiedallelesLD_MT.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

#Shinnecock rareified alleles
meltar3 %>%
  filter(variable %in%c("Shin_1_2016","Shin_2_2016","Shin_1_2017","Shin_2_2017"))%>%
  mutate(GRP = as.factor(variable)) %>%
  ggplot(aes(x=fct_rev(GRP),y=value))+
  geom_boxplot(aes(fill=GRP),fill= "#016450")+ coord_flip()+ 
  xlab(' ')+ylab("Rareified allele count")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('rareifiedallelesLD_dhin.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)


#not going to do private alleles or heatmaps by loci for these. But I will do wilcoxon test, friedmans test and heatmaps. 
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
fdata <- loc_stats_shin
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

llocstats <- loc_stats_shin ##CHANGE HERE
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
  wilcox <- wilcoxsign_test(Ht ~ GRP | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,
                            zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)),
                   "test" = "Ht")  ### change this one. 
  results <- bind_rows(results, df)}

#now change test and re_run above. 
results_nei <-results
results_fis <-results
results_even <-results
results_shannon <-results

results <-bind_rows(results_nei, results_fis, results_even, results_shannon)
results <- results %>% dplyr::select(-temp)
results <-mutate(results, bonferroni=0.05/6) # number of pairwise comparisons ####CHECK THIS
results <-mutate(results, significance = ifelse(p.value>=bonferroni,"not significant","significant"))

write.csv(results,file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/wilcoxSHI.csv" )

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
ggsave('neiSHItile.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

#FIS heatmap
#delete duplicate pairs to form the half grid of the heatmap. 
results_Fis <-filter(results, test=="Fis") %>%
  arrange(test.statistic) %>% unique() %>% 
  mutate(pair1 = pmin(pop1,pop2), pair2 =pmax(pop1,pop2)) %>% arrange(pair1)
#toDelete <- seq(1, nrow(results_Fis), 2)
#results_Fis <- results_Fis[ toDelete ,]

# Filled by test statistic: FIS
results_Fis %>%
  ggplot(aes(x = pair1, y = pair2))+
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  scale_fill_manual(values=cols)+
  geom_text(aes(label = round(p.value,3),color=significance), size=5)+
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
ggsave('FisSHItile.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

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
ggsave('evennessSHItile.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

#shannons heatmap
#delete duplicate pairs to form the half grid of the heatmap. 
results_shannon <-filter(results, test=="Shannon") %>%
  arrange(test.statistic) %>% unique() %>% 
  mutate(pair1 = pmin(pop1,pop2), pair2 =pmax(pop1,pop2)) %>% arrange(pair1)
#toDelete <- seq(1, nrow(results_shannon), 2)
#results_shannon <- results_shannon[ toDelete ,]

# Filled by test statistic: evenness
results_shannon %>%
  ggplot(aes(x = pair1, y = pair2))+
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  scale_fill_manual(values=cols)+
  geom_text(aes(label = round(p.value,3),color=significance), size=5)+
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
ggsave('shannonSHItile.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)


##Wilcoxon test for rarified alleles##
#list of all the different possible pairs
meltarMT <- filter(meltar3, variable %in% c("Mt_2","Mt_1","Mt_3","Mt_4"))
meltarSHI <- filter(meltar3, variable %in% c("Shin_1_2016", "Shin_2_2016", "Shin_1_2017", "Shin_2_2017"))

comp <- as.character(unique(meltarMT$variable))  
pairs <- expand.grid(comp, comp) %>%
  filter(!Var1 == Var2) %>%
  rownames_to_column("PAIR") %>%
  split(.$PAIR) %>%
  purrr::map(function(x){
    x %>%
      dplyr::select(-PAIR) %>%
      gather(key = meltarMT, value = variable, 1:2) %>%
      dplyr::select(-meltarMT)})
meltarMT <-as.data.frame(meltarMT) %>% mutate(GRP=as.factor(variable))
# empty data frame for results
results_ar <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("pop1", "pop2", "stat", "temp", "p.value"))
n <- as.numeric(length(pairs))

library("coin")
# loop over pairs
for(p in 1:n){
  pair <- pairs[[p]]$variable
  temp <- meltarMT %>%
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



results_ar <-mutate(results_ar, bonferroni=0.05/6)
results_ar <-mutate(results_ar, significance = ifelse(p.value>=bonferroni,"not significant","significant")) %>% mutate(starsig=ifelse(significance=="significant",round(p.value,4),NA),test.statistic=abs(stat))
results_ar<-mutate(results_ar, p_value = cut(p.value, breaks=c(0,0.005, 0.01,0.05,0.1,0.5,1)))
cols <- c("(0,0.005]"="#034e7b", "(0.005,0.01]" = "#045a8d", "(0.01,0.05]" = "#2b8cbe", "(0.05,0.1]" = "#74a9cf", "(0.1,0.5]"  = "#a6bddb", "(0.5,1]"="#d0d1e6")

results_arMT <-arrange(results_ar, test.statistic) %>% unique() %>% 
  mutate(pair1 = pmin(pop1,pop2), pair2 =pmax(pop1,pop2)) %>% arrange(pair1)
#toDelete <- seq(1, nrow(results_arMT), 2)
#results_arMT <- results_arMT[ toDelete ,]

# plot heatmap of results ====
results_arMT %>%
  ggplot(aes(x = pair1, y = pair2))+
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  scale_fill_manual(values=cols)+
  geom_text(aes(label = round(p.value,3),color=significance), size=5)+
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
ggsave('rariefied_allelesMTtile.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

### Rareified alleles heatmap shinnecock
comp <- as.character(unique(meltarSHI$variable))  
pairs <- expand.grid(comp, comp) %>%
  filter(!Var1 == Var2) %>%
  rownames_to_column("PAIR") %>%
  split(.$PAIR) %>%
  purrr::map(function(x){
    x %>%
      dplyr::select(-PAIR) %>%
      gather(key = meltarSHI, value = variable, 1:2) %>%
      dplyr::select(-meltarSHI)})
meltarSHI <-as.data.frame(meltarSHI) %>% mutate(GRP=as.factor(variable))
# empty data frame for results
results_ar <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("pop1", "pop2", "stat", "temp", "p.value"))
n <- as.numeric(length(pairs))

library("coin")
# loop over pairs
for(p in 1:n){
  pair <- pairs[[p]]$variable
  temp <- meltarSHI %>%
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



results_ar <-mutate(results_ar, bonferroni=0.05/6)
results_ar <-mutate(results_ar, significance = ifelse(p.value>=bonferroni,"not significant","significant")) %>% mutate(starsig=ifelse(significance=="significant",round(p.value,4),NA),test.statistic=abs(stat))
results_ar<-mutate(results_ar, p_value = cut(p.value, breaks=c(0,0.005, 0.01,0.05,0.1,0.5,1)))
cols <- c("(0,0.005]"="#034e7b", "(0.005,0.01]" = "#045a8d", "(0.01,0.05]" = "#2b8cbe", "(0.05,0.1]" = "#74a9cf", "(0.1,0.5]"  = "#a6bddb", "(0.5,1]"="#d0d1e6")

results_arSHI <-arrange(results_ar, test.statistic) %>% unique() %>% 
  mutate(pair1 = pmin(pop1,pop2), pair2 =pmax(pop1,pop2)) %>% arrange(pair1)
#toDelete <- seq(1, nrow(results_arMT), 2)
#results_arMT <- results_arMT[ toDelete ,]

# plot heatmap of results ====
results_arSHI %>%
  ggplot(aes(x = pair1, y = pair2))+
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  scale_fill_manual(values=cols)+
  geom_text(aes(label = round(p.value,3),color=significance), size=5)+
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
ggsave('rariefied_allelesSHItile.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

###### FST ########
#Weir and Cockheram FST global & pairwise. - probably your best bet. 
library("strataG")

wf.s2 <- genind2gtypes(wf.shin)
popStruct.S <- popStructTest(wf.s2, nrep = 1000, quietly = TRUE)
popStruct.S

#FST heatmap
ShiFST <-as.data.frame(popStruct.S$pairwise$result) %>% dplyr::select(strata.1, strata.2, Fst, Fst.p.val)%>%
mutate(p_value = cut(Fst.p.val, breaks=c(0,0.001,0.005,0.01,0.05,0.1,0.5,1)), significance=ifelse(Fst.p.val < 0.005, "sig","not.sig"))
cols <- c("(0,0.001]"="#023858","(0.001,0.005]"="#034e7b", "(0.001,0.01]" = "#045a8d", "(0.01,0.05]" = "#2b8cbe", "(0.05,0.1]" = "#74a9cf", "(0.1,0.5]"  = "#a6bddb", "(0.5,1]"="#d0d1e6")

# Filled by  p value, and reporting the fst value. 
ShiFST %>%
  ggplot(aes(x = strata.1, y = strata.2))+
  #geom_tile(aes(fill=Fst), show.legend = TRUE)+ # Filled by test statistic
  #geom_tile(aes(fill=p.value), show.legend = TRUE)+ #filled  by p.value (gradient)
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  #scale_fill_gradient(low ="#023858" , high = "#a6bddb", space = "Lab", na.value = "white", guide = "colourbar", aesthetics = "fill")+ #gradient fill
  #scale_fill_viridis_c()+
  scale_fill_manual(values=cols)+
  #geom_text(aes(label = round(Fst.p.val,3)), color="white", size=5)+  #label is p value
  geom_text(aes(label = round(Fst,3),color="signifcance"), size=5)+  #label is FST
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
ggsave('FSTShitileP.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)


#explore structure within Mattituck
wf.M2 <- genind2gtypes(wf.mt)
popStruct.M <- popStructTest(wf.M2, nrep = 1000, quietly = TRUE)
popStruct.M

#FST heatmap
MtFST <-as.data.frame(popStruct.M$pairwise$result) %>% dplyr::select(strata.1, strata.2, Fst, Fst.p.val)%>%
  mutate(p_value = cut(Fst.p.val, breaks=c(0,0.001,0.005,0.01,0.05,0.1,0.5,1)), significance=ifelse(Fst.p.val < 0.005, "sig","not.sig"))
cols <- c("(0,0.001]"="#023858","(0.001,0.005]"="#034e7b", "(0.001,0.01]" = "#045a8d", "(0.01,0.05]" = "#2b8cbe", "(0.05,0.1]" = "#74a9cf", "(0.1,0.5]"  = "#a6bddb", "(0.5,1]"="#d0d1e6")

  MtFST %>%
    ggplot(aes(x = strata.1, y = strata.2))+
    #geom_tile(aes(fill=Fst), show.legend = TRUE)+ # Filled by test statistic
    #geom_tile(aes(fill=p.value), show.legend = TRUE)+ #filled  by p.value (gradient)
    geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
    #scale_fill_gradient(low ="#023858" , high = "#a6bddb", space = "Lab", na.value = "white", guide = "colourbar", aesthetics = "fill")+ #gradient fill
    #scale_fill_viridis_c()+
    scale_fill_manual(values=cols)+
    #geom_text(aes(label = round(Fst.p.val,3)), color="white", size=5)+  #label is p value
    geom_text(aes(label = round(Fst,3),color="signifcance"), size=5)+  #label is FST
    scale_color_manual(values=c("white","red"), guide=FALSE)+
    xlab("")+ylab("")+
    theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
ggsave('FSTMttile.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)


#MATTITUCK but separate year cohorts
setPop(wf.mt) <-~Bay/Con/Year
wf.M3 <-genind2gtypes(wf.mt)
popStruct.M3 <- popStructTest(wf.M3, nrep = 1000, quietly = TRUE)
popStruct.M3
#FST heatmap
MtFST3 <-as.data.frame(popStruct.M3$pairwise$result) %>% dplyr::select(strata.1, strata.2, Fst, Fst.p.val)%>%
  mutate(p_value = cut(Fst.p.val, breaks=c(0,0.001,0.005,0.01,0.05,0.1,0.5,1)), significance=ifelse(Fst.p.val < 0.005, "sig","not.sig"))
cols <- c("(0,0.001]"="#023858","(0.001,0.005]"="#034e7b", "(0.001,0.01]" = "#045a8d", "(0.01,0.05]" = "#2b8cbe", "(0.05,0.1]" = "#74a9cf", "(0.1,0.5]"  = "#a6bddb", "(0.5,1]"="#d0d1e6")

MtFST3 %>%
  ggplot(aes(x = strata.1, y = strata.2))+
  #geom_tile(aes(fill=Fst), show.legend = TRUE)+ # Filled by test statistic
  #geom_tile(aes(fill=p.value), show.legend = TRUE)+ #filled  by p.value (gradient)
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  #scale_fill_gradient(low ="#023858" , high = "#a6bddb", space = "Lab", na.value = "white", guide = "colourbar", aesthetics = "fill")+ #gradient fill
  #scale_fill_viridis_c()+
  scale_fill_manual(values=cols)+
  #geom_text(aes(label = round(Fst.p.val,3)), color="white", size=5)+  #label is p value
  geom_text(aes(label = round(Fst,3),color="signifcance"), size=5)+  #label is FST
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
ggsave('FSTMt3tileP.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)

#MATTITUCK but adults are combined and YOY are not!
wf.mt3
wf.mt3Com <-genind2gtypes(wf.mt3)
popStruct.MCom <- popStructTest(wf.mt3Com, nrep = 1000, quietly = TRUE)
popStruct.MCom

MtFSTCom <-as.data.frame(popStruct.MCom$pairwise$result) %>% dplyr::select(strata.1, strata.2, Fst, Fst.p.val)%>%
  mutate(p_value = cut(Fst.p.val, breaks=c(0,0.001,0.005,0.01,0.05,0.1,0.5,1)), significance=ifelse(Fst.p.val < 0.005, "sig","not.sig"))
cols <- c("(0,0.001]"="#023858","(0.001,0.005]"="#034e7b", "(0.001,0.01]" = "#045a8d", "(0.01,0.05]" = "#2b8cbe", "(0.05,0.1]" = "#74a9cf", "(0.1,0.5]"  = "#a6bddb", "(0.5,1]"="#d0d1e6")

MtFSTCom %>%
  ggplot(aes(x = strata.1, y = strata.2))+
  #geom_tile(aes(fill=Fst), show.legend = TRUE)+ # Filled by test statistic
  #geom_tile(aes(fill=p.value), show.legend = TRUE)+ #filled  by p.value (gradient)
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  #scale_fill_gradient(low ="#023858" , high = "#a6bddb", space = "Lab", na.value = "white", guide = "colourbar", aesthetics = "fill")+ #gradient fill
  #scale_fill_viridis_c()+
  scale_fill_manual(values=cols)+
  #geom_text(aes(label = round(Fst.p.val,3)), color="white", size=5)+  #label is p value
  geom_text(aes(label = round(Fst,3),color="signifcance"), size=5)+  #label is FST
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
ggsave('FSTMtCombinedAdultstileP.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)

#just look at MT YOY
wf.mtyoy <-popsub(wf.mt3, blacklist=c("Mt_3_both", "Mt_4_both"))
wf.mty <-genind2gtypes(wf.mtyoy)
popStruct.y <- popStructTest(wf.mty, nrep = 1000, quietly = TRUE)
popStruct.y

MtFSTy <-as.data.frame(popStruct.y$pairwise$result) %>% dplyr::select(strata.1, strata.2, Fst, Fst.p.val)%>%
  mutate(p_value = cut(Fst.p.val, breaks=c(0,0.001,0.005,0.01,0.05,0.1,0.5,1)), significance=ifelse(Fst.p.val < 0.005, "sig","not.sig"))
cols <- c("(0,0.001]"="#023858","(0.001,0.005]"="#034e7b", "(0.001,0.01]" = "#045a8d", "(0.01,0.05]" = "#2b8cbe", "(0.05,0.1]" = "#74a9cf", "(0.1,0.5]"  = "#a6bddb", "(0.5,1]"="#d0d1e6")

MtFSTy %>%
  ggplot(aes(x = strata.1, y = strata.2))+
  #geom_tile(aes(fill=Fst), show.legend = TRUE)+ # Filled by test statistic
  #geom_tile(aes(fill=p.value), show.legend = TRUE)+ #filled  by p.value (gradient)
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  #scale_fill_gradient(low ="#023858" , high = "#a6bddb", space = "Lab", na.value = "white", guide = "colourbar", aesthetics = "fill")+ #gradient fill
  #scale_fill_viridis_c()+
  scale_fill_manual(values=cols)+
  #geom_text(aes(label = round(Fst.p.val,3)), color="white", size=5)+  #label is p value
  geom_text(aes(label = round(Fst,3),color="signifcance"), size=5)+  #label is FST
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
ggsave('FSTMtyoyonlyP.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)
 #nothing here. 




####### IR ############### We haven't done this yet but will come back to this. 

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
ggsave('rel_bay.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

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
ggsave('IRBaytile.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

### Genetic Distance ###

## Tree's using provesti's distance

#Within Shinnecock
set.seed(999)
wf.s %>%
  genind2genpop(pop = ~Bay/Con/Year) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = provesti.dist, missingno="ignore")



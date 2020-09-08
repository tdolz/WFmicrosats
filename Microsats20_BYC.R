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
#wfia.pair <- seppop(wfpop2CLEAN) %>% lapply(pair.ia) #by bay!
#ggsave("rawLD20.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()

#create a second dataset that removes one loci in the WF27/WF33 pair 
#also remove WF06 & WF32 due to null alllels, see null_checker script and "microsat record Sept2020.ppt" for details. 
#go back to wfpop for this. 
wfpopLD <-genclone2genind(wfpop)
all_loci <- locNames(wfpopLD)# create vector of all loci
removeloc <- c("WF27","WF06","WF32")# create vector containing loci to remove
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
setPop(wfpopLD) <-~Bay
wfhwe.pop <- seppop(wfpopLD) %>% lapply(hw.test)
#write.csv(wfhwe.pop, file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/wfhwepop16.csv")
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

#Test HWE over Con/Year
#Test HWE over Year population
setPop(wfpopLD) <-~Bay/Con/Year
wfhwe.pop <- seppop(wfpopLD) %>% lapply(hw.test)
write.csv(wfhwe.pop, file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/wfhwepop16BAYCONYEAR.csv")
(wfhwe.mat <- sapply(wfhwe.pop, "[", i = TRUE, j = 3)) # Take the third column with all rows
wfhw.mc <-sapply(wfhwe.pop, "[", i = TRUE, j = 4) #the PR exact based on the Monte carlo test! 
wfhw.mc
alpha  <- 0.05
newmat <- wfhwe.mat
newmat[newmat > alpha] <- 1
levelplot(t(newmat),scales=list(x=list(rot=90)))#slightly different if you do it with genind object (which it calls for) instead of genclone.
#ggsave("HWEtest20Bayconyear.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()


### Summary Data ####
setPop(wfpopLD) <-~Bay
toto <-summary(wfpopLD)
barplot(toto$loc.n.all, ylab="Number of alleles",las=2,
        main="Number of alleles per locus")
barplot(toto$Hexp-toto$Hobs, main="Heterozygosity: expected-observed",ylab="Hexp - Hobs",las=2, ylim=c(-0.05, 0.06))
barplot(toto$n.by.pop, main="Sample sizes per population", ylab="Number of genotypes",las=3)
toto

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

# write file with genetic diversity stats by locus to file
#write_delim(loc_stats, "results/gendiv.locstats", delim = "\t")


#create different loc stats groups for testing. 
loc_stats_bay <-filter(loc_stats, GRP %in% c("Nap","Mor","Jam","Shin","Mt"))
loc_stats_BYC <-filter(loc_stats, GRP %in% c("Nap_6_2016","Mor_6_2016","Jam_6_2016","Shin_1_2016","Shin_2_2016","Mt_2_2015","Mt_1_2015","Mt_2_2016",  
                                             "Mt_1_2016","Mt_3_2015","Mt_4_2015","Mt_3_2016","Mt_4_2016","Shin_1_2017","Shin_2_2017")) #does not include MT 5 because we don't know what those individuals are. 
loc_stats_MT <-filter(loc_stats, GRP %in% c("Mt_2", "Mt_1", "Mt_3", "Mt_4")) #does not include MT 5 because we don't know what those individuals are. 
loc_stats_shin <-filter(loc_stats, GRP %in% c("Shin_1_2016", "Shin_2_2016","Shin_1_2017","Shin_2_2017"))
loc_stats_16 <-filter(loc_stats, GRP %in% c("Nap_2016","Mor_2016","Shin_2016","Mt_2016","Jam_2016"))

#generic melt
meltlocstats <-pivot_longer(loc_stats,cols = c("N_ALLELES", "SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis","SHANNON_IDX","STODD_TAYLOR_IDX"),
                            names_to="variable", values_to="value")
#bay melt
meltlocstats_bay <-pivot_longer(loc_stats_bay,cols = c("N_ALLELES", "SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis","SHANNON_IDX","STODD_TAYLOR_IDX"),
                                names_to="variable", values_to="value")


meltlocstats2 <-filter(meltlocstats, GRP %in% c("Nap","Mor","Jam","Shin","Mt","ALL")) 
#this is essentially by bay but also includes the "ALL" Group

##### Box and whisker plots for diversity stats #####
drabcolors <-c("#d0d1e6","#a6bddb", "#67a9cf", "#1c9099", "#016450")

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
ggsave('nei_bay.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

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
ggsave('FIS_bay.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

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
ggsave('Evenness_bay.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

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
ggsave('Shannon_idxbay.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)


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
#ggsave('rareifiedallelesLD_bay.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)

##Private alleles
df <-genind2df(wfpopLD, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(34,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33)]
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
#ggsave('privateallelesLD_bay.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 5)


##### Heatmaps #####
#heatmap
ggplot(meltpA, aes(x = LOCUS, y = variable, fill=value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_viridis() +
  coord_fixed(ratio = 1) +
  ylab("")+
  #theme_standard() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
#ggsave('privateallelesLDheat.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)

#heatmap of FIS
meltlocstats %>%
  filter(variable == "Fis") %>%
  ggplot(aes(x = LOCUS, y = GRP, fill=value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(value, 2)),colour="lightblue") +
  scale_fill_viridis(option="A") +
  coord_fixed(ratio = 1) +
  ylab("")+
  #theme_standard() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  
ggsave('Fis.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

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
ggsave('neiBaytile.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

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
ggsave('FisBaytile.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

#EVENNESS heatmap
#delete duplicate pairs to form the half grid of the heatmap. 
results_even <-filter(results, test=="EVENNESS") %>%
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
ggsave('evennessBaytile.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)

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
ggsave('shannonBaytile.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)


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
  mutate(pair1 = pmin(pop1,pop2), pair2 =pmax(pop1,pop2)) %>% arrange(pair1)
toDelete <- seq(1, nrow(results_ar), 2)
results_ar <- results_ar[ toDelete ,]

# plot heatmap of results ====
results_ar %>%
  ggplot(aes(x = pair1, y = pair2))+
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  scale_fill_manual(values=cols)+
  geom_text(aes(label = round(p.value,3),color=significance), size=5)+
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
ggsave('rariefied_allelesBaytile.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 10, height = 7)



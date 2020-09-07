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
meltlocstats2 %>%
  filter(variable == "Ht") %>%
  #mutate(GRP = fct_rev(as.factor(GRP))) %>%
  ggplot(aes(x=GRP,y=value))+
  #ggplot(aes(x=fct_inorder(GRP),y=value))+
  geom_boxplot(fill="lightgray")+ 
  coord_flip()+ 
  #ylim(0.7,1.0)+
  xlab(' ')+ylab("Nei's Gene Diversity")+theme_cowplot()+guides(fill = FALSE, colour = FALSE) 
ggsave('nei.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)

meltlocstats2 %>%
  filter(variable == "Fis") %>%
  #mutate(GRP = fct_rev(as.factor(GRP))) %>%
  ggplot(aes(x=GRP,y=value))+
  #ggplot(aes(x=fct_inorder(GRP),y=value))+
  geom_boxplot(fill="lightgray")+ 
  coord_flip()+ 
  #ylim(0.7,1.0)+
  xlab(' ')+ylab("FIS")+theme_cowplot()+guides(fill = FALSE, colour = FALSE) 
ggsave('FIS.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)

meltlocstats2 %>%
  filter(variable == "EVENNESS") %>%
  #mutate(GRP = fct_rev(as.factor(GRP))) %>%
  ggplot(aes(x=GRP,y=value))+
  #ggplot(aes(x=fct_inorder(GRP),y=value))+
  geom_boxplot(fill="lightgray")+ 
  coord_flip()+ 
  #ylim(0.7,1.0)+
  xlab(' ')+ylab("Evenness")+theme_cowplot()+guides(fill = FALSE, colour = FALSE) 
ggsave('Evenness.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)

meltlocstats2 %>%
  filter(variable == "SHANNON_IDX") %>%
  #mutate(GRP = fct_rev(as.factor(GRP))) %>%
  ggplot(aes(x=GRP,y=value))+
  #ggplot(aes(x=fct_inorder(GRP),y=value))+
  geom_boxplot(fill="lightgray")+ 
  coord_flip()+ 
  #ylim(0.7,1.0)+
  xlab(' ')+ylab("Shannon Index")+theme_cowplot()+guides(fill = FALSE, colour = FALSE) 
ggsave('Shannon_idx.png', path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)


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
meltar <- melt(ar, id.vars=c("LOCUS"), measure.vars=c("Mt","Shin","Nap","Mor","Jam","ALL"))

#meltar <-filter(meltar, variable != "ALL")

meltar %>%
  filter(variable !="ALL")%>%
  #mutate(variable = fct_rev(as.factor(variable))) %>%
  ggplot(aes(x=variable,y=value))+
  geom_boxplot(fill="lightgray")+ 
  coord_flip()+ 
  xlab(' ')+ylab('rareified allele count')+theme_cowplot()+guides(fill = FALSE, colour = FALSE)
#ggsave('rareifiedallelesLD.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)

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

meltpA <- melt(pA, id.vars=c("LOCUS"), measure.vars=c("Mt","Shin","Nap","Mor","Jam"))

meltpA %>%
  mutate(variable = fct_rev(as.factor(variable))) %>%
  ggplot(aes(x=variable,y=value))+
  geom_boxplot(fill="lightgray")+ 
  coord_flip()+ 
  xlab(' ')+ylab('private alleles')+theme_cowplot()+guides(fill = FALSE, colour = FALSE)
#ggsave('privateallelesLD.png',path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs", width = 7, height = 7)

library("viridis")
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

#for rarified alleles, you have to make a new Allelic Richness analysis for each group you are comparing. 
meltar2 <- pivot_longer(ar, cols=c("Mt","Shin","Nap","Mor","Jam","ALL"),names_to="variable", values_to="value")
meltar3 <-filter(meltar2, variable != "ALL") #it doenst like it when you exclude the "all" category... 
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
results <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("pop1", "pop2", "stat", "temp", "p.value"))

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
                   "p-value" = as.numeric(pvalue(wilcox)))
  results <- bind_rows(results, df)}

results <- results %>%
  dplyr::select(-temp)
#write_delim(results, "results/exphet_est.wilcox")

results <-mutate(results, bonferroni=0.05/10) # number of pairwise comparisons ####CHECK THIS
results <-mutate(results, significance = ifelse(p.value>=bonferroni,"not significant","significant"))

#heatmap of results (new way)
results <- mutate(results, starsig=ifelse(significance=="significant",round(p.value,4),NA),newstat=abs(stat))
#write.csv(results,file="/Users/tdolan/Documents/WIP research/microsats/microsats_results/wilcox.csv" )

#delete duplicate pairs
results <-arrange(results, newstat)
toDelete <- seq(1, nrow(results), 2)
results2 <- results[ toDelete ,]

#we're not interested in the wilcoxon stat, we're interested in the p value
results2 %>%
  ggplot(aes(x = pop1, y = pop2))+
  geom_tile(aes(fill=newstat), show.legend = TRUE)+
  #scale_fill_viridis(option="D") +
  scale_fill_gradient(low = "#034e7b", high = "#f1eef6",
                      space = "Lab", na.value = "white", guide = "colourbar",
                      aesthetics = "fill")+
  #geom_text(aes(label = round(p.value,4),color=significance), size=3)+
  geom_text(aes(label = starsig), color="black", size=3)+
  xlab("")+ylab("")+
  theme(axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))

#though I don't think we need this heatmap, I think we do letter groups on the bar charts? Idk how you would do this. 
#you could say groups with different letters are in different groups. 

##Wilcoxon test for rarified alleles##
meltarr <-filter(meltar, variable %in% c("Nap","Mor","Jam","Shin", "Mt"))

#list of all the different possible pairs
comp <- as.character(unique(meltarr$variable))  
pairs <- expand.grid(comp, comp) %>%
  filter(!Var1 == Var2) %>%
  rownames_to_column("PAIR") %>%
  split(.$PAIR) %>%
  purrr::map(function(x){
    x %>%
      dplyr::select(-PAIR) %>%
      gather(key = meltarr, value = variable, 1:2) %>%
      dplyr::select(-meltarr)})

# empty data frame for results
results <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("pop1", "pop2", "stat", "temp", "p.value"))
n <- as.numeric(length(pairs))

library("coin")
# loop over pairs
for(p in 1:n){
  pair <- pairs[[p]]$variable
  temp <- meltarr %>%
    dplyr::filter(variable %in% pair) %>%
    mutate(GRP = ordered(variable, levels = pair),
           LOCUS = as.factor(LOCUS)) %>%
    droplevels()
  wilcox <- wilcoxsign_test(value ~ variable | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,
                            zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)))
  results <- bind_rows(results, df)}

results <- results %>%
  dplyr::select(-temp)
write_delim(results, "results/exphet_est.wilcox")

results <-mutate(results, bonferroni=0.05/9)
results <-mutate(results, significance = ifelse(p.value>=bonferroni,"not significant","significant"))

# plot heatmap of results ====
ggplot(results, aes(x = pop1, y = pop2, fill = stat)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(stat,2),color=significance)) +
  scale_fill_viridis(option="D") +
  coord_fixed(ratio = 1) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Wilcox Signed-Rank Test: Rareified Alleles ")  ### CHANGE THE TITLE
#ggsave('rareallelesLD.png', width = 7, height = 7)  #### CHANGE THE FILE NAME





##### FST ######

#Weir and Cockheram FST global & pairwise. - probably your best bet. 
library("strataG")
setPop(wfpopLD) <-~Bay/Con
wf.g2 <-genind2gtypes(wfpopLD)
#statFst(wf.g)
popStruct <- popStructTest(wf.g2, stats = c(statFst, statFstPrime), nrep = 1000, quietly = TRUE)
popStruct

#some custom for mike
setPop(wfpopLD) <- ~Bay/Con

#explore structure within a bay using strataG
setPop(wfpopLD) <-~Bay/Con/Year
wf.s <-popsub(wfpopLD, c("Shin_1_2016","Shin_2_2016","Shin_1_2017", "Shin_2_2017"))
wf.s2 <- genind2gtypes(wf.s)
popStruct.S <- popStructTest(wf.s2, stats = c(statFst, statFstPrime), nrep = 1000, quietly = TRUE)
popStruct.S

bayyrcon <-as.data.frame(popStruct.S$pairwise$pair.mat$Fst)
write.csv(bayyrcon,file="bayyerfst.csv")

#explore structure within Mattituck
setPop(wfpopLD) <-~Bay/Con/Year
wf.M <-popsub(wfpopLD, c("Mt_1_2015","Mt_2_2015","Mt_1_2016", "Mt_2_2016", "Mt_3_2015", "Mt_4_2015", "Mt_5_2015", "Mt_3_2016", "Mt_4_2016", "Mt_5_2016"))
wf.M2 <- genind2gtypes(wf.M)
popStruct.M <- popStructTest(wf.M2, stats = c(statFst, statFstPrime), nrep = 1000, quietly = TRUE)
popStruct.M



#with the pruned dataset from shin
popStruct.nosibs <- popStructTest(wfpopLDnosibs, stats = c(statFst, statFstPrime), nrep = 1000, quietly = TRUE)
popStruct.nosibs

#see previous scripts for other kinds of FST you can calculate

# ultimately, you will need help interpreting these. 

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

#violin plot
p <- ggplot(rel2, aes(x=pop, y=IR)) + 
  geom_violin(trim=FALSE,fill='#A4A4A4', color="darkred")+
  stat_summary(fun.data="mean_sdl", geom="pointrange")

#barplot as before
ggplot(rel2, aes(x=pop,y=IR))+
  #ggplot(aes(x=fct_inorder(GRP),y=value))+
  geom_boxplot(fill="lightgray")+ 
  coord_flip()+ 
  #ylim(0.7,1.0)+
  xlab(' ')+ylab("IR")+theme_cowplot()+guides(fill = FALSE, colour = FALSE) 
#ggsave('IRld.png', width = 7, height = 7)

## test for significant differences in mean IR for each bay-pair, with t test. Print only the p values. 
jam <-filter(rel2,pop=="Jam")
shin <-filter(rel2,pop=="Shin")
nap <-filter(rel2,pop=="Nap")
mor <-filter(rel2,pop=="Mor")
mt <-filter(rel2,pop=="Mt")

print("Nap vs. Shin")
var.test(nap$IR,shin$IR)
t.test(nap$IR,shin$IR,var.equal=TRUE)$p.value
print("Nap vs. Mor")
var.test(nap$IR,mor$IR)
t.test(nap$IR,mor$IR,var.equal=TRUE)$p.value
print("Jam vs. Nap")
var.test(jam$IR,nap$IR)
t.test(jam$IR,nap$IR,var.equal=TRUE)$p.value
print("Mor vs. Shin")
var.test(mor$IR,shin$IR)
t.test(mor$IR,shin$IR,var.equal=TRUE)$p.value
print("Jam vs. Shin")
var.test(jam$IR,shin$IR)
t.test(jam$IR,shin$IR,var.equal=TRUE)$p.value
print("Jam vs. Mor")
var.test(jam$IR,mor$IR)
t.test(jam$IR,mor$IR,var.equal=TRUE)$p.value

print("Jam vs. Mt")
var.test(jam$IR,mt$IR)
t.test(jam$IR,mt$IR,var.equal=TRUE)$p.value
print("Shin vs. Mt")
var.test(shin$IR,mt$IR)
t.test(shin$IR,mt$IR,var.equal=TRUE)$p.value
print("Mor vs. Mt")
var.test(mor$IR,mt$IR)
t.test(mor$IR,mt$IR,var.equal=TRUE)$p.value
print("Nap vs. Mt")
var.test(nap$IR,mt$IR)
t.test(nap$IR,mt$IR,var.equal=TRUE)$p.value
#no significant differences between bays. 
#hmmm I think there is a more efficient way to code this. 

#### Genetic Distance ####

## Tree's using provesti's distance

#Within Shinnecock
set.seed(999)
wf.s %>%
  genind2genpop(pop = ~Bay/Con/Year) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = provesti.dist, missingno="ignore")

#all bays?
set.seed(999)
wfpopLD %>%
  genind2genpop(pop = ~Bay) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = provesti.dist, missingno="ignore")
#This one makes a lot of sense

#all bays?
set.seed(999)
wfpopLD %>%
  genind2genpop(pop = ~Bay/Con/Year) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = provesti.dist, missingno="ignore")
#This one makes a lot of sense


#### Isolation by distance.####
library(MASS)
setPop(wfpopLD) <-~Bay
wfpop.geo <-genind2genpop(wfpop2)
Dgen <-dist.genpop(wfpop.geo,method=2) #edward's euclidean distance, but other methods are available. like provesti, nei, etc.

#Create a coordinates list:
coord.dist =matrix(c(41.008611,72.048611,40.791667,72.696111,40.606389,73.876389,40.845278,72.500556),nrow=4,ncol=2,byrow=TRUE)
rownames(coord.dist) <-c("Nap","Mor","Jam","Shin")
colnames(coord.dist)<-c("x","y")

#i'm not sure if it's ok to just put in your own euclidean distance matrix like this but it seems to work. #previously i thought you had to attach it to the genpop objecte as $other, maybe as xy coordinates? I am not sure if it understands that it corresponds to population, but let's plot anyway? 
Dgeo <-dist(coord.dist)

ibd2 <-mantel.randtest(Dgen,Dgeo)
ibd2
plot(ibd2)


dens <- kde2d(Dgeo,Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, Dgen, xlab="distance (KM)", ylim= c(0.1, 0.5),ylab="genetic distance (euclidean)")
image(dens, col=transp(myPal(300),0.7), add=TRUE)
abline(lm(Dgen~Dgeo))
title("Isolation by distance plot")
#strong looking pattern of IBD. 
library("car")
Anova(lm(Dgen~Dgeo))

#try with ggplot


##### AMOVA #####
setPop(wfpopLD) <- ~Bay
gapless <-missingno(wfpop2) # only works with gapless data. 
strat.amo <-poppr.amova(gapless,~Bay)
strat.amo

#randomization test
library("ade4")
set.seed(1999)
stratsig  <- ade4::randtest(strat.amo, nrepet = 999)
plot(stratsig)
stratsig

#another AMOVA, this one from pegas. 
setPop(wfpop2) <- ~Bay #idk how to change it so it doesn't care about Con/Year
bov_dist  <- dist(wfpop2) # Euclidean distance
bov_stra  <- strata(wfpop2)
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
locus_table(wfpop2)
wflt.pair <- seppop(wfpop2) %>% lapply(locus_table) #by bay!

#diversity boot
mlgtab <-mlg.table(wfpop2)
divboot <- diversity_boot(mlgtab, n=1000, H=TRUE, G=TRUE, lambda=TRUE, E5=TRUE)
dci <-diversity_ci(wfpop2, n=10000, ci=95, plot=TRUE, center=TRUE, rarefy=FALSE, raw=TRUE)

#with rarefaction

dci <-diversity_ci(wfpop2, n=10000, ci=95, plot=FALSE, center=TRUE, rarefy=TRUE, raw=TRUE)


#other stats
otherstats<-diff_stats(wfpop2)
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

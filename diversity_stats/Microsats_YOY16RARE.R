### Genetic analysis Bays ###
## created 5-30-23 From Microsats_20.R
## The same as MicrosatsYOY16.R but with rareifaction. 

## NOTES ##

#install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)
library("tidyverse")
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
library("poolr") #for fisher test, might get rid. 


#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")


##### Formating the dataset #####
# We are going to use the doubl0 version. 
##### Formating the dataset #####
wfpop <- read.genalex("./data/popcorrect_17_sept20204genalex_doubl0ABC.csv")
#wfpop4df <-read.csv("./data/popcorrect_17_sept2020_doubl0ABC.csv", header = TRUE) #csv

splitStrata(wfpop) <-~Ocean/Bay/Con/Year
setPop(wfpop) <-~Bay

## Data cleanup. For details see datacleanup_comparison.R ##
wfpopLD <-genclone2genind(wfpop)
all_loci <- locNames(wfpopLD)
removeloc <- c("WF06","WF27", "WF32")# we also may want to remove A441 and WF32 for missing data reasons. 
#especially WF32 because it's 70% missing in moriches and is out of HWE.
#Come back to this
keeploc <- setdiff(all_loci, removeloc)
wfpopLD <- wfpopLD[loc = keeploc]
wfpopLD <- wfpopLD %>% missingno("geno", cutoff=0.25) 
wfpopLD <-wfpopLD %>% missingno("loci", cutoff=0.25) #removes WF01
setPop(wfpopLD) <-~Bay
info_table(wfpopLD, plot = TRUE, scaled =FALSE)
ggsave("missing_by_loci.png", path="./diversity_stats/diversity_figs/YOY16_rare")

##Subset to YOY16 only##
#remove all individuals that aren't 2016 YOY
#first save an original version of the dataset
wfpopLD.og <-wfpopLD
setPop(wfpopLD) <-~Bay/Year
wfpop <-popsub(wfpopLD, exclude=c("Mt_2015","Mt_adults","Shin_2017"))
setPop(wfpopLD) <-~Bay

#HWE Heatmap#
setPop(wfpopLD) <-~Ocean
hw.test(wfpopLD, B=1000) #permutation based
hw.test(wfpop2CLEAN, B=1000)
hw.test(wfpopLD, B=0) #analytical p value
wfhwe.pop <- seppop(wfpopLD) %>% lapply(hw.test)
write.csv(wfhwe.pop, file="./diversity_stats/diversity_output_files/YOY16_rare/wfhwepop16.csv")
(wfhwe.mat <- sapply(wfhwe.pop, "[", i = TRUE, j = 3)) # Take the third column with all rows ---> output this for supplementary tables.
wfhw.mc <-sapply(wfhwe.pop, "[", i = TRUE, j = 4) #the PR exact based on the Monte carlo test! ----> p.values on the hw.test
wfhw.mc  # this is just the p values. 
alpha  <- 0.05
newmat <- wfhwe.mat
newmat[newmat > alpha] <- 1 #where the p value on the chi square is greater than 0.05, give it a 1.
# so pink means zero, which means p < 0.05, which means out of HWE. 
levelplot(t(newmat),scales=list(x=list(rot=90)))
ggsave("HWEtest16.png", path="./diversity_stats/diversity_figs/YOY16_rare")
dev.off()

#PopGenReport - Remember to turn this off. 
#setPop(wfpopLD) <-~Bay
#wf.gen <-genclone2genind(wfpopLD) 
#popgenreport(wf.gen,mk.counts=TRUE,mk.locihz = TRUE, mk.fst=TRUE, mk.allele.dist=TRUE, mk.null.all=TRUE,mk.allel.rich = TRUE,mk.differ.stats = TRUE,path.pgr=getwd(),mk.Rcode=TRUE,mk.pdf=TRUE )

### Summary Data ####

### Expected-observed heterozygosity with rarefaction ####
setPop(wfpopLD) <-~Bay
#Set up the loop 
df <-genind2df(wfpopLD, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(ncol(df),1:(ncol(df)-1))]
df[df=="NA"] <- 0 # missing data must be 0

########### RAREIFICATION LOOP 1 ###########################
#randomly sample from within shinnecock and mattituck 100x
rare.out <-list()
for(i in 1:100){
  df.split <-split(df, df$pop)
  new.shin <-sample_n(df.split$Shin,30)
  new.mt <-sample_n(df.split$Mt,30)
  rare.bays <-bind_rows(new.shin,new.mt,df.split$Jam,df.split$Mor,df.split$Nap)
  wf.rare <-df2gtypes(rare.bays,ploidy=2)
  wf.rare <-gtypes2genind(wf.rare)
  #wfpopLD.rare <-df2genind(rare.bays, ind.names = rare.bays$Ind, pop=rare.bays$pop, ploidy=2, sep=",")
#perform calculation
  toto <-summary(wf.rare)
  hexhobs <-as.data.frame(toto$Hexp-toto$Hobs) 
  hexhobs <- tibble::rownames_to_column(hexhobs,"locus")
  names(hexhobs) <-c("locus","difference")
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
  h4 <- bind_rows(naph, morh, Mth, jamh, shinh)%>%mutate(count=i)
  rare.out[[i]]<-h4
}
########### END RAREIFICATION LOOP 1 ###########################

mean.h4 <-bind_rows(rare.out)%>%
  mutate(locus=as.factor(locus), Bay=as.factor(Bay))%>%
  group_by(locus,Bay)%>%
  dplyr::summarize(mean.hexp=mean(hexp), sd.hexp=sd(hexp),
            mean.hobs=mean(hobs), sd.hobs=sd(hobs),
            mean.hexhobs=mean(hexhobs), sd.hexhobs=sd(hexhobs),.groups="drop")
write.csv(mean.h4, file="./diversity_stats/diversity_output_files/YOY16_bay/exp_obs_heterozygosity_noWF32.csv")
#write.csv(mean.h4, file="./diversity_stats/diversity_output_files/YOY16_bay/exp_obs_heterozygosity.csv")


#heatmap of observed minus expected. 
ggplot(mean.h4, aes(x = locus, y = Bay, fill=mean.hexhobs)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(mean.hexhobs, 3)), color="grey") +
  scale_fill_viridis(option="plasma") +
  coord_fixed(ratio = 1) +
  ylab("")+
  #theme_standard() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
#ggsave('hexhobsheat17.png',path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 7)
ggsave('hexhobsheat17.png',path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 7)


##### Shannon's summary stats ########
#https://gist.github.com/sjoleary/cdc32efbfd50c96eef446ebb7c2f2387

# genind object of genotypes
wfpopLD

#sample metadata
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


########### RAREIFICATION LOOP 2 ###########################

grp_sumlocstats <-list()
grp_locstats <-list()
gen_grp <-list()
for (j in 1:100){

setPop(wfpopLD.og)<-~Ocean #the ocean should be all individuals - do not rareify?. 
gen_oce <-seppop(wfpopLD.og)

# BAY COMPARISON WITH 2016 YOY ONLY - need to rareify
setPop(wfpopLD)<-~Bay
df <-genind2df(wfpopLD, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(ncol(df),1:(ncol(df)-1))]
df[df=="NA"] <- 0 
df.split <-split(df, df$pop)
new.shin <-sample_n(df.split$Shin,30)
new.mt <-sample_n(df.split$Mt,30)
rare.bays <-bind_rows(new.shin,new.mt,df.split$Jam,df.split$Mor,df.split$Nap)
#wfpopLD <-df2genind(rare.bays, pop=pop, ploidy=2)
wf.rare16 <-df2gtypes(rare.bays,ploidy=2)
wf.rare16 <-gtypes2genind(wf.rare16)
#setPop(wf.rare16)<-~Bay
gen_bay <-seppop(wf.rare16)

#Bay/con The meaningful comparison here is between migrant and resident adults in mattituck
# so we will rareify migrant adults. 
# The shinnecock & Mattituck bay/cohort matters, but the sample sizes are roughly equal already. 
setPop(wfpopLD.og)<-~Bay/Con
df <-genind2df(wfpopLD.og, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(ncol(df),1:(ncol(df)-1))]
df[df=="NA"] <- 0 
new.mt_3 <-filter(df, pop =="Mt_3")%>% sample_n(15)
rare.mt3 <-filter(df, pop != "Mt_3")%>%bind_rows(new.mt_3)
rare.mt3 <-df2gtypes(rare.mt3, ploidy=2)%>%gtypes2genind()
gen_con <-seppop(rare.mt3)

#The meaningful comparisons here are between cohorts within years.
#The year factor of the MT adults is not useful here. 
setPop(wfpopLD.og)<-~Bay/Con/Year
#df <-genind2df(wfpopLD.og, usepop = TRUE, oneColPerAll = TRUE) 
#df$Ind <- rownames(df)
#df <-df[,c(ncol(df),1:(ncol(df)-1))]
#df[df=="NA"] <- 0 
#df.split <-split(df,df$pop)
# the mattituck cohorts in 2015 have very low sample sizes. 
# Maybe we shouldn't even include them?
# otherwise these seem pretty equal? 
# We could rareify the Shinnecock cohorts down to 20? but it's not a huge discrepancy. 
gen_bayconyear <-seppop(wfpopLD)

#This is not a comparison that we use in the paper, as of now. 
#so we will not bother to rareify it.
setPop(wfpopLD.og)<-~Year
gen_year <- seppop(wfpopLD)

#We are interested in whether the year effect trumps the cohort effect. 
setPop(wfpopLD)<-~Bay/Year
gen_bayyear <- seppop(wfpopLD)

### On this version of the sheet we are using ONLY 2016 YOY FOR THE BAY LEVEL#####

gen_grp <-c(gen_oce,
            gen_bay,
            gen_con,
            gen_bayconyear,
            gen_year)
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
  filter(LOCUS != "mean")%>%
  mutate(iter=j)

loc_stats[is.na(loc_stats)] <- NA
pivlocstats <-pivot_longer(loc_stats,cols = c("N_ALLELES", "SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis","SHANNON_IDX","STODD_TAYLOR_IDX"),
                             names_to="variable", values_to="value")
sumlocstats <-pivlocstats%>% group_by(GRP,variable)%>%summarize(meanvar=mean(value), sdvar=sd(value),.groups="drop")%>%
  mutate(iter=j, var=sdvar^2)

grp_sumlocstats[[j]] <-sumlocstats
grp_locstats[[j]] <-loc_stats
}
########### END RAREIFICATION LOOP 2 #########################################################
#our output from this loop is a list of dataframes. Now we need to combine the dataframes by averaging the averages. 

grp_sumlocstats <-bind_rows(grp_sumlocstats)
grp_locstats <-bind_rows(grp_locstats)

sumlocstats_mean <-grp_sumlocstats %>% group_by(GRP, variable) %>% 
  summarize(av_meanvar = mean(meanvar), av_sd =sqrt(sum(var)/n_distinct(grp_sumlocstats$iter)), .groups="drop")

locstats_mean <-grp_locstats%>%group_by(GRP,LOCUS)%>% 
  summarize(N_ALLELES=mean(N_ALLELES),EVENNESS=mean(EVENNESS), Ho=mean(Ho), Hs=mean(Hs),Ht=mean(Ht),Fis=mean(Fis),
            SHANNON_IDX=mean(SHANNON_IDX), SIMPSON_IDX=mean(SIMPSON_IDX), STODD_TAYLOR_IDX=mean(STODD_TAYLOR_IDX),.groups="drop")

write.csv(locstats_mean,file="./diversity_stats/diversity_output_files/YOY16_rare/mean_loc_stats17.csv")
write.csv(sumlocstats_mean,file="./diversity_stats/diversity_output_files/YOY16_rare/mean_sumlocstats17.csv")


###plot it
allcolors <-c("ALL"="grey","Atl"="grey","Jam"="#d0d1e6","Jam_6"="#d0d1e6","Mor"="#a6bddb","Mor_6"="#a6bddb",
              "Mt"="#67a9cf","Mt_1"="#67a9cf","Mt_2"="#67a9cf","Mt_3"="#67a9cf","Mt_4"="#67a9cf", "Mt_5"="#67a9cf",
              "Nap"="#1c9099","Nap_6"="#1c9099","Shin"="#016450","Shin_1"="#016450","Shin_2"="#016450")
            

#new piv
pivlocstats <-pivot_longer(locstats_mean,cols = c("N_ALLELES", "SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis","SHANNON_IDX","STODD_TAYLOR_IDX"),
                           names_to="variable", values_to="value")

#Giant barplot
pivlocstats %>%
  mutate(GRP=as.factor(GRP)) %>%
  filter(variable %in% c("N_ALLELES", "SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis")) %>%
  #filter(GRP %in% c("Atl","Nap","Mor","Jam","Shin","Shin_1_2016","Shin_2_2016","Shin_1_2017","Shin_2_2017","Mt","Mt_3","Mt_4","Mt_1_2015","Mt_2_2015","Mt_1_2016","Mt_2_2016"))%>%
  ggplot(aes(x=fct_rev(GRP),y=value),fill=GRP)+
  geom_boxplot(aes(fill=GRP))+ 
  scale_fill_manual(name = "GRP",values = allcolors)+
  coord_flip()+ 
  facet_wrap(~variable, scales="free_x", nrow=2)+ 
  theme(axis.text = element_text(size = 10),axis.title = element_text(size = 12),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('pivlocstats17.png', path="./diversity_stats/diversity_figs/YOY16_rare", width = 12, height = 5)

#create different loc stats groups for testing. 
loc_stats_bay <-filter(locstats_mean, GRP %in% c("Nap","Mor","Jam","Shin","Mt"))

#bay melt
meltlocstats_bay <-pivot_longer(loc_stats_bay,cols = c("N_ALLELES", "SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis","SHANNON_IDX","STODD_TAYLOR_IDX"),
                                names_to="variable", values_to="value")


meltlocstats2 <-filter(locstats_mean, GRP %in% c("Nap","Mor","Jam","Shin","Mt","ALL"))%>%
  pivot_longer(cols = c("N_ALLELES", "SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis","SHANNON_IDX","STODD_TAYLOR_IDX"),
               names_to="variable", values_to="value")
#this is essentially by bay but also includes the "ALL" Group

##### Box and whisker plots for diversity stats #####
drabcolors <-c("#d0d1e6","#a6bddb", "#67a9cf", "#1c9099", "#016450")
#drabcolors <-c("gray","#d0d1e6","#a6bddb", "#67a9cf", "#1c9099", "#016450")

#Nei's gene diversity
#meltlocstats2 %>%
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
ggsave('nei_bay17.png', path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 5)

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
ggsave('FIS_bay17.png', path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 5)

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
ggsave('Evenness_bay17.png', path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 5)

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
ggsave('Shannon_idxbay17.png', path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 5)

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
ggsave('Hsbay17.png', path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 5)

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
ggsave('Hobay17.png', path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 5)


######## CALCULATE RAREFIED ALLELIC RICHNESS ----
# differences in sample size can bias the number of alleles sampled in a population
# calculate allelic richness corrected for sample size using rarefaction


################################ START RAREIFACTION LOOP 3 ###############################################

grp_ar <-list()

for (j in 1:100){
# overall ====

#rareify
setPop(wfpopLD)<-~Bay
  df <-genind2df(wfpopLD, usepop = TRUE, oneColPerAll = TRUE) 
  df$Ind <- rownames(df)
  df <-df[,c(ncol(df),1:(ncol(df)-1))]
  df[df=="NA"] <- 0 
  df.split <-split(df, df$pop)
  new.shin <-sample_n(df.split$Shin,30)
  new.mt <-sample_n(df.split$Mt,30)
  rare.bays <-bind_rows(new.shin,new.mt,df.split$Jam,df.split$Mor,df.split$Nap)
  #wfpopLD <-df2genind(rare.bays, pop=pop, ploidy=2)
  wf.rare16 <-df2gtypes(rare.bays,ploidy=2)
  wf.rare16 <-gtypes2genind(wf.rare16)
#reattach the strata
  str16 <-data.frame(rep("Atl",146), wf.rare16@pop)
  names(str16)<-c("Ocean","Bay")
  strata(wf.rare16) <-str16
  
#allelic richness ocean level  
setPop(wf.rare16) <- ~Ocean
dat <- hierfstat:::genind2hierfstat(wf.rare16)
ar <- allelic.richness(dat,diploid = TRUE)

ar <- as.data.frame(ar$Ar) %>%
  rownames_to_column("LOCUS") %>%
  dplyr::rename(ALL = 2)%>% 
  select(-dumpop)

# By Bay ===
setPop(wf.rare16) <- ~Bay
dat <- hierfstat:::genind2hierfstat(wf.rare16, pop = wf.rare16@pop)
df <- allelic.richness(dat,diploid = TRUE)

df <- as.data.frame(df$Ar) %>%
  rownames_to_column("LOCUS") #%>%dplyr::rename(Nap = V1,Mor = V2,Jam = V3,Shin = V4,Mt = V5)
ar <- left_join(ar, df)%>%mutate(iter=j)
grp_ar[[j]]<-ar
}


################################# END RAREIFACTION LOOP 3 #########################################
ar <-bind_rows(grp_ar)

ar_mean <-ar%>%group_by(LOCUS)%>% 
  summarize(ALL=mean(ALL),Jam=mean(Jam),Mor=mean(Mor),Mt=mean(Mt),Nap=mean(Nap),Shin=mean(Shin))

ar_sd <-ar%>%group_by(LOCUS)%>% 
  summarize(ALL=sqrt(sum(ALL)/n_distinct(ar$iter)),Jam=sqrt(sum(Jam)/n_distinct(ar$iter)),
            Mor=sqrt(sum(Mor)/n_distinct(ar$iter)),Mt=sqrt(sum(Mt)/n_distinct(ar$iter)),
            Nap=sqrt(sum(Nap)/n_distinct(ar$iter)),Shin=sqrt(sum(Shin)/n_distinct(ar$iter)))

write.csv(ar_mean, "./diversity_stats/diversity_output_files/YOY16_rare/rarefied_allelecount_mean.csv")
write.csv(ar_sd, "./diversity_stats/diversity_output_files/YOY16_rare/rarefied_allelecount_SD.csv")


#visualize results as boxplot. 

#for rarified alleles, you have to make a new Allelic Richness analysis for each group you are comparing. 
meltar2 <- pivot_longer(ar_mean, cols=c("Mt","Shin","Nap","Mor","Jam","ALL"),names_to="variable", values_to="value")
meltar3 <-filter(meltar2, variable != "ALL") #it doenst like it when you exclude the "all" category... 
meltar_all <-filter(meltar2, variable == "ALL")

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
ggsave('rareifiedallelesLD_bay17.png',path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 5)

############################# START RAREIFACTION LOOP 4 ############################################
##Private alleles

grp_pa <-list()

for (j in 1:100){
#rareify
  setPop(wfpopLD)<-~Bay
  df <-genind2df(wfpopLD, usepop = TRUE, oneColPerAll = TRUE) 
  df$Ind <- rownames(df)
  df <-df[,c(ncol(df),1:(ncol(df)-1))]
  df[df=="NA"] <- 0 
  df.split <-split(df, df$pop)
  new.shin <-sample_n(df.split$Shin,30)
  new.mt <-sample_n(df.split$Mt,30)
  rare.bays <-bind_rows(new.shin,new.mt,df.split$Jam,df.split$Mor,df.split$Nap)
  #wfpopLD <-df2genind(rare.bays, pop=pop, ploidy=2)
  wf.rare16 <-df2gtypes(rare.bays,ploidy=2)
  wf.rare16 <-gtypes2genind(wf.rare16) 
  #setPop(wf.rare16)<-~Bay

df <-genind2df(wf.rare16, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df<-df[,c(ncol(df),1:(ncol(df)-1))]
df[df=="NA"] <- 0 # missing data must be 0
wf.g <-df2gtypes(df,ploidy=2)

pA <-as.data.frame(privateAlleles(wf.g))%>%mutate(iter=j)
setDT(pA,keep.rownames=TRUE)
colnames(pA)[1] <- "LOCUS"

grp_pa[[j]] <-pA
}
############################# END RAREIFACTION LOOP 4 ############################################

pA <-bind_rows(grp_pa)

pA_roundedmean <-pA%>%group_by(LOCUS)%>% 
  summarize(Jam=round(mean(Jam)),Mor=round(mean(Mor)),
            Mt=round(mean(Mt)),Nap=round(mean(Nap)),Shin=round(mean(Shin)))

pA_sd <-pA%>%group_by(LOCUS)%>% 
  summarize(Jam=sqrt(sum(Jam)/n_distinct(pA$iter)),
            Mor=sqrt(sum(Mor)/n_distinct(pA$iter)),Mt=sqrt(sum(Mt)/n_distinct(pA$iter)),
            Nap=sqrt(sum(Nap)/n_distinct(pA$iter)),Shin=sqrt(sum(Shin)/n_distinct(pA$iter)))

write.csv(pA_roundedmean, "./diversity_stats/diversity_output_files/YOY16_rare/privateAlleles_roundedmean.csv")
write.csv(pA_sd, "./diversity_stats/diversity_output_files/YOY16_rare/privateAlleles_SD.csv")


meltpA <- pivot_longer(pA_roundedmean, cols=c("Mt","Shin","Nap","Mor","Jam"),names_to="variable", values_to="value")

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
ggsave('privateallelesLD_bay17.png',path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 5)

##heatmap of private alleles by locus
ggplot(meltpA, aes(x = LOCUS, y =factor(variable, level=c('Jam', 'Mor', 'Shin', 'Nap',"Mt")), fill=value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = value), color="grey") +
  scale_fill_viridis(option="plasma") +
  coord_fixed(ratio = 1) +
  ylab("BAY")+
  #theme_standard() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave('privateAllelesHeat.png',path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 7)


######## TEST FOR SIGNIFICANT DIFFERENCES#####
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

#To do this the rareifaction way, we perform the test on each iteration of the rareification 
# then we perform fisher's test to adjust the p value. 
# I think this is more appropriate than performing the test on the mean values. 


#set level to test for friedmans
fdata <- filter(grp_locstats, GRP %in% c("Nap","Mor","Jam","Shin","Mt")) #Bay
#return to list form
fdata <-split(fdata, fdata$iter)

#turn the friedman's tests into a function.
#####Start function########### 

fri_test <-function(fdata){
  #df = unname(fdata)
#####Friedman's test for Global heterogeneity#####
a<-friedman.test(Hs~GRP | LOCUS, data= fdata) #Hs - fixing this. 
b<-friedman.test(EVENNESS~GRP | LOCUS, data= fdata) #evenness
c<-friedman.test(Ht~GRP | LOCUS, data= fdata) #Ht
d<-friedman.test(Fis~GRP | LOCUS, data= fdata) #Fis
e<-friedman.test(SHANNON_IDX~GRP | LOCUS, data= fdata) 
f<-friedman.test(value~variable | LOCUS, data= meltar3) #rareified alleles #come back to this one.... 
#h4 from before. 
g<-friedman.test(hobs~Bay | locus, data=h4)

p.vals <-c(a$p.value, b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value)
comparisons <-c(a$data.name,b$data.name,c$data.name,d$data.name,e$data.name,f$data.name,g$data.name)

friedmans_tests <-cbind(comparisons, p.vals)%>%as.data.frame
}

#### END FUNCTION #####

#Apply the function to df in the list#
fri_grp <- fdata %>% lapply(fri_test)
fri_grp <-bind_rows(fri_grp, .id="iter")
#save raw values
write.csv(fri_grp,"./diversity_stats/diversity_output_files/YOY16_rare/friedmanstests_raw.csv")

#apply the fisher test to the comparisons.
# we will do bonferroni as well just to be safe.  
fisher <-mutate(fri_grp, p.vals=as.numeric(p.vals))%>%
                  split(fri_grp$comparisons)
fisher_p <-list()
fisher_pvals <-c()
bonferroni_pvals <-c()
stouffer_pvals <-c()
simes_pvals <-c()
for (i in 1:length(fisher)){
  df <-fisher[[i]]
  #Empirical Browns test
  #The browns method is better because it accounts for nonindependence 
  #https://www.bioconductor.org/packages/devel/bioc/html/EmpiricalBrownsMethod.html 
 # test_stat <-word(df$comparisons, start=1, end=1)#get the comparison
 # fd <-dplyr::select(fdata[[i]],c(GRP, LOCUS,test_stat[1]))%>%
 #   pivot_wider(names_from = LOCUS, values_from = test_stat[1])%>%as.data.frame() #create the matrix
 # brown_p[[i]] <-empiricalBrownsMethod(data_matrix = fd, p_values = )
  
  #fisher test
  fisher_p[[i]] <-fisher(df$p.vals, adjust="none")
  fisher_pvals[i]<-fisher_p[[i]]$p
  #bonferroni test
  bp <-bonferroni(df$p.vals, adjust="none")
  bonferroni_pvals[i]<-bp$p
  #stouffer test
  sp <-stouffer(df$p.vals, adjust="none")
  stouffer_pvals[i]<-sp$p
  #simes test)
  simes_p <-metapod::combineParallelPValues(as.list(df$p.vals), method="simes")
  simes_pvals[i]<-simes_p$p.value
}
names(fisher_p)<-names(fisher)
fisher_pvals<-data.table(fisher_pvals,bonferroni_pvals,stouffer_pvals, simes_pvals, comparison=names(fisher))
#However, this fisher test is for independent tests which these are not. 

write.csv(fisher_pvals,"./diversity_stats/diversity_output_files/YOY16_rare/friedmanstests_adjustpval.csv")


#######Wilcoxon tests #####
#tell us which groups are significantly different. 
# you have to set which level you are comparing. 
  #the different levels we have to work with are:
    #loc_stats_bay <- bays
    #loc_stats_BYC <- bay year con
    #loc_stats_MT <-bay con for mattituck only, no group 5 (unidentified adults)
    #loc_stats_shin <- bay year con shin
    #loc_stats_16 <- 2016 YOY at bay level. 

#from before, fdata is the list of df of loc_stats_bay
#llocstats <- loc_stats_bay ## This IS 2016 yoy at bay level because we previously removed the other groups. 
#list of all the different possible pairs


##################### START WILCOXON FUNCTION #######################################
w_test <-function(llocstats){
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
results_nei<-bind_rows(results_list)

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
results_fis<-bind_rows(results_list)

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
results_even<-bind_rows(results_list)

#SHANNON_IDX
results_list <-list()
for(p in 1:n){
  #tryCatch({
  pair <- pairs[[p]]$GRP
  temp <- llocstats %>%
    dplyr::filter(GRP %in% pair) %>%
    mutate(GRP = ordered(GRP, levels = pair),
           LOCUS = as.factor(LOCUS)) %>%
    droplevels()
  wilcox <- wilcoxsign_test(SHANNON_IDX ~ GRP | LOCUS,    ###### MUST REMEMBER TO CHANGE WHAT TEST YOU ARE DOING. 
                            data = temp,zero.method = "Pratt")
  df <- data.frame("pop1" = pair[1], 
                   "pop2" = pair[2], 
                   "stat" = as.numeric(wilcox@statistic@teststatistic), 
                   "p-value" = as.numeric(pvalue(wilcox)),
                   "test" = "SHANNON_IDX") 
  #if(!is.null(df)){
  results_list[[p]] <- df #}#, error=function(e){})
}
#results_list <-results_list[-which(sapply(results_list, is.null))]
results_shannon<-bind_rows(results_list)

#combine
results <-bind_rows(results_nei, results_fis, results_even, results_shannon)
}
############# END WILCOXON FUNCTION ####################################

##Apply wilcoxon function
wil_test <- fdata%>% lapply(w_test)
wil_df <-bind_rows(wil_test, .id="iter")
#save raw values
write.csv(wil_df,"./diversity_stats/diversity_output_files/YOY16_rare/wilcoxontests_raw.csv")

#apply the fisher test to the comparisons.
# we will do bonferroni as well just to be safe.  
wt <-mutate(wil_df, p.value=as.numeric(p.value))%>%
  split(wil_df$test)
fisher_p <-list()
fisher_pvals <-c()
bonferroni_pvals <-c()
stouffer_pvals <-c()
simes_pvals <-c()
for (i in 1:length(wt)){
  df <-wt[[i]]
  #fisher
  fisher_p[[i]] <-fisher(df$p.value, adjust="none")
  fisher_pvals[i]<-fisher_p[[i]]$p
  #bonferroni
  bp <-bonferroni(df$p.value, adjust="none")
  bonferroni_pvals[i]<-bp$p
  #stouffer
  sp <-stouffer(df$p.value, adjust="none")
  stouffer_pvals[i]<-sp$p
  #simes
  simes_p <-metapod::combineParallelPValues(as.list(df$p.value), method="simes")
  simes_pvals[i]<-simes_p$p.value
}
wilcoxon_fisher_pvals<-data.table(fisher_pvals,bonferroni_pvals,stouffer_pvals,simes_pvals,comparison=names(wt))
#However, this fisher test is for independent tests which these are not. 
write.csv(wilcoxon_fisher_pvals,"./diversity_stats/diversity_output_files/YOY16_rare/wilcoxon_adjustpval.csv")


################## WE STOPPED HERE 8/18/23 WITH THE RAREFACTION CONVERSION #########################



############ MAKE PLOTS #############################################
#results <- results %>% dplyr::select(-temp)
#results <-mutate(results, bonferroni=0.05/10) # number of pairwise comparisons ####CHECK THIS
#results <-mutate(results, significance = ifelse(p.value>=bonferroni,"not significant","significant"))
#write.csv(results,file="./diversity_stats/diversity_output_files/YOY16_rare/wilcox.csv")

#heatmap of results (new way)
#results <- mutate(results, starsig=ifelse(significance=="significant",round(p.value,4),NA),test.statistic=abs(stat))
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
ggsave('neiBaytile17.png',path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 7)

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
ggsave('FisBaytile17.png',path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 7)

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
ggsave('evennessBaytile17.png',path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 7)

#shannons heatmap
#delete duplicate pairs to form the half grid of the heatmap. 
results_shannon <-filter(results, test=="SHANNON_IDX") %>%
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
ggsave('shannonBaytile17.png',path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 7)


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
  results_ar <- bind_rows(results_ar)

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
ggsave('rariefied_allelesBaytile17.png',path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 7)



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
ggsave('IRld.png',path="./diversity_stats/diversity_figs/YOY16_rare", width = 7, height = 7)

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
ggsave('rel_bay17.png', path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 5)



##### FINAL GIANT BARPLOT ########

allcolors <-c("grey","#d0d1e6","#a6bddb","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#67a9cf","#1c9099","#016450","#016450","#016450","#016450","#016450")
#Giant barplot
pivlocstats %>%
  mutate(GRP=as.factor(GRP)) %>%
  filter(variable %in% c("N_ALLELES", "SIMPSON_IDX", "EVENNESS","Ho","Hs","Ht","Fis")) %>%
  filter(GRP %in% c("Atl","Nap","Mor","Jam","Shin","Shin_1_2016","Shin_2_2016","Shin_1_2017","Shin_2_2017","Mt","Mt_3","Mt_4","Mt_1_2015","Mt_2_2015","Mt_1_2016","Mt_2_2016"))%>%
  ggplot(aes(x=fct_rev(GRP),y=value),fill=GRP)+
  geom_boxplot(aes(fill=GRP))+ 
  scale_fill_manual(name = "GRP",values = allcolors)+
  coord_flip()+ 
  facet_wrap(~variable, scales="free_x", nrow=2)+ 
  theme(axis.text = element_text(size = 10),axis.title = element_text(size = 12),panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),plot.margin=margin(0.5,1,0.5,0.5,"cm"))+guides(fill = FALSE, colour = FALSE) 
ggsave('pivlocstats17.png', path="./diversity_stats/diversity_figs/YOY16_rare", width = 12, height = 5)


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
ggsave('IRBaytile17.png',path="./diversity_stats/diversity_figs/YOY16_rare", width = 10, height = 7)

#### Genetic Distance ####

## Tree's using provesti's distance

#all bays?
set.seed(999)
jpeg(file="./diversity_stats/diversity_figs/YOY16_rare/provesti.jpeg")
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


ibd2 <-mantel.randtest(Dgen,Dgeo2, nrepet=1000)
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



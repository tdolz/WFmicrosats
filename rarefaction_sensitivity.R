# Rarefaction sensitivity
# 
# #install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)
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
wfpop4df <-read.csv("./data/popcorrect_17_sept2020_doubl0ABC.csv", header = TRUE) #csv

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

##Subset to YOY16 only##
#remove all individuals that aren't 2016 YOY
#first save an original version of the dataset
wfpopLD.og <-wfpopLD
setPop(wfpopLD) <-~Bay/Year
wfpop <-popsub(wfpopLD, exclude=c("Mt_2015","Mt_adults","Shin_2017"))
setPop(wfpopLD) <-~Bay

########### RAREIFICATION LOOP 30 ###########################

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
  
  ### On this version of the sheet we are using ONLY 2016 YOY FOR THE BAY LEVEL#####
  
  gen_grp <-c(gen_oce,
              gen_bay)
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
  sumlocstats <-pivlocstats%>% group_by(GRP,variable)%>%dplyr::summarize(meanvar=mean(value), sdvar=sd(value),.groups="drop")%>%
    mutate(iter=j, var=sdvar^2)
  
  grp_sumlocstats[[j]] <-sumlocstats
  grp_locstats[[j]] <-loc_stats
}
########### END RAREIFICATION LOOP 30 #########################################################
#our output from this loop is a list of dataframes. Now we need to combine the dataframes by averaging the averages. 

grp_sumlocstats30 <-bind_rows(grp_sumlocstats)
grp_locstats30 <-bind_rows(grp_locstats)

sumlocstats_mean30 <-grp_sumlocstats30 %>% group_by(GRP, variable) %>% 
  dplyr::summarize(av_meanvar = mean(meanvar), av_sd =sqrt(sum(var)/n_distinct(grp_sumlocstats$iter)), .groups="drop")

locstats_mean30 <-grp_locstats30%>%group_by(GRP,LOCUS)%>% 
  dplyr::summarize(N_ALLELES=mean(N_ALLELES),EVENNESS=mean(EVENNESS), Ho=mean(Ho), Hs=mean(Hs),Ht=mean(Ht),Fis=mean(Fis),
            SHANNON_IDX=mean(SHANNON_IDX), SIMPSON_IDX=mean(SIMPSON_IDX), STODD_TAYLOR_IDX=mean(STODD_TAYLOR_IDX),.groups="drop")


########################RAREFICATION 25##############

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
  new.shin <-sample_n(df.split$Shin,25)
  new.mt <-sample_n(df.split$Mt,25)
  rare.bays <-bind_rows(new.shin,new.mt,df.split$Jam,df.split$Mor,df.split$Nap)
  #wfpopLD <-df2genind(rare.bays, pop=pop, ploidy=2)
  wf.rare16 <-df2gtypes(rare.bays,ploidy=2)
  wf.rare16 <-gtypes2genind(wf.rare16)
  #setPop(wf.rare16)<-~Bay
  gen_bay <-seppop(wf.rare16)
  
  ### On this version of the sheet we are using ONLY 2016 YOY FOR THE BAY LEVEL#####
  
  gen_grp <-c(gen_oce,
              gen_bay)
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
  sumlocstats <-pivlocstats%>% group_by(GRP,variable)%>%dplyr::summarize(meanvar=mean(value), sdvar=sd(value),.groups="drop")%>%
    mutate(iter=j, var=sdvar^2)
  
  grp_sumlocstats[[j]] <-sumlocstats
  grp_locstats[[j]] <-loc_stats
}
########### END RAREIFICATION LOOP 25 #########################################################
#our output from this loop is a list of dataframes. Now we need to combine the dataframes by averaging the averages. 

grp_sumlocstats25 <-bind_rows(grp_sumlocstats)
grp_locstats25 <-bind_rows(grp_locstats)

sumlocstats_mean25 <-grp_sumlocstats25 %>% group_by(GRP, variable) %>% 
  dplyr::summarize(av_meanvar = mean(meanvar), av_sd =sqrt(sum(var)/n_distinct(grp_sumlocstats$iter)), .groups="drop")

locstats_mean25 <-grp_locstats25%>%group_by(GRP,LOCUS)%>% 
  dplyr::summarize(N_ALLELES=mean(N_ALLELES),EVENNESS=mean(EVENNESS), Ho=mean(Ho), Hs=mean(Hs),Ht=mean(Ht),Fis=mean(Fis),
            SHANNON_IDX=mean(SHANNON_IDX), SIMPSON_IDX=mean(SIMPSON_IDX), STODD_TAYLOR_IDX=mean(STODD_TAYLOR_IDX),.groups="drop")

#####################################RAREFICATION LOOP 20
#####################################
########### RAREIFICATION LOOP 20 ###########################

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
  new.shin <-sample_n(df.split$Shin,20)
  new.mt <-sample_n(df.split$Mt,20)
  rare.bays <-bind_rows(new.shin,new.mt,df.split$Jam,df.split$Mor,df.split$Nap)
  #wfpopLD <-df2genind(rare.bays, pop=pop, ploidy=2)
  wf.rare16 <-df2gtypes(rare.bays,ploidy=2)
  wf.rare16 <-gtypes2genind(wf.rare16)
  #setPop(wf.rare16)<-~Bay
  gen_bay <-seppop(wf.rare16)
  
  ### On this version of the sheet we are using ONLY 2016 YOY FOR THE BAY LEVEL#####
  
  gen_grp <-c(gen_oce,
              gen_bay)
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
  sumlocstats <-pivlocstats%>% group_by(GRP,variable)%>%dplyr::summarize(meanvar=mean(value), sdvar=sd(value),.groups="drop")%>%
    mutate(iter=j, var=sdvar^2)
  
  grp_sumlocstats[[j]] <-sumlocstats
  grp_locstats[[j]] <-loc_stats
}
########### END RAREIFICATION LOOP 20 #########################################################
#our output from this loop is a list of dataframes. Now we need to combine the dataframes by averaging the averages. 

grp_sumlocstats20 <-bind_rows(grp_sumlocstats)
grp_locstats20 <-bind_rows(grp_locstats)

sumlocstats_mean20 <-grp_sumlocstats20 %>% group_by(GRP, variable) %>% 
  dplyr::summarize(av_meanvar = mean(meanvar), av_sd =sqrt(sum(var)/n_distinct(grp_sumlocstats$iter)), .groups="drop")

locstats_mean20 <-grp_locstats20%>%group_by(GRP,LOCUS)%>% 
  dplyr::summarize(N_ALLELES=mean(N_ALLELES),EVENNESS=mean(EVENNESS), Ho=mean(Ho), Hs=mean(Hs),Ht=mean(Ht),Fis=mean(Fis),
            SHANNON_IDX=mean(SHANNON_IDX), SIMPSON_IDX=mean(SIMPSON_IDX), STODD_TAYLOR_IDX=mean(STODD_TAYLOR_IDX),.groups="drop")


########### RAREIFICATION LOOP 15 ###########################

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
  new.shin <-sample_n(df.split$Shin,15)
  new.mt <-sample_n(df.split$Mt,15)
  rare.bays <-bind_rows(new.shin,new.mt,df.split$Jam,df.split$Mor,df.split$Nap)
  #wfpopLD <-df2genind(rare.bays, pop=pop, ploidy=2)
  wf.rare16 <-df2gtypes(rare.bays,ploidy=2)
  wf.rare16 <-gtypes2genind(wf.rare16)
  #setPop(wf.rare16)<-~Bay
  gen_bay <-seppop(wf.rare16)
  
  ### On this version of the sheet we are using ONLY 2016 YOY FOR THE BAY LEVEL#####
  
  gen_grp <-c(gen_oce,
              gen_bay)
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
  sumlocstats <-pivlocstats%>% group_by(GRP,variable)%>%dplyr::summarize(meanvar=mean(value), sdvar=sd(value),.groups="drop")%>%
    mutate(iter=j, var=sdvar^2)
  
  grp_sumlocstats[[j]] <-sumlocstats
  grp_locstats[[j]] <-loc_stats
}
########### END RAREIFICATION LOOP 15 #########################################################
#our output from this loop is a list of dataframes. Now we need to combine the dataframes by averaging the averages. 

grp_sumlocstats15 <-bind_rows(grp_sumlocstats)
grp_locstats15 <-bind_rows(grp_locstats)

sumlocstats_mean15 <-grp_sumlocstats15 %>% group_by(GRP, variable) %>% 
  dplyr::summarize(av_meanvar = mean(meanvar), av_sd =sqrt(sum(var)/n_distinct(grp_sumlocstats$iter)), .groups="drop")

locstats_mean15 <-grp_locstats15%>%group_by(GRP,LOCUS)%>% 
  dplyr::summarize(N_ALLELES=mean(N_ALLELES),EVENNESS=mean(EVENNESS), Ho=mean(Ho), Hs=mean(Hs),Ht=mean(Ht),Fis=mean(Fis),
            SHANNON_IDX=mean(SHANNON_IDX), SIMPSON_IDX=mean(SIMPSON_IDX), STODD_TAYLOR_IDX=mean(STODD_TAYLOR_IDX),.groups="drop")

########### RAREIFICATION LOOP 10 ###########################

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
  new.shin <-sample_n(df.split$Shin,10)
  new.mt <-sample_n(df.split$Mt,10)
  rare.bays <-bind_rows(new.shin,new.mt,df.split$Jam,df.split$Mor,df.split$Nap)
  #wfpopLD <-df2genind(rare.bays, pop=pop, ploidy=2)
  wf.rare16 <-df2gtypes(rare.bays,ploidy=2)
  wf.rare16 <-gtypes2genind(wf.rare16)
  #setPop(wf.rare16)<-~Bay
  gen_bay <-seppop(wf.rare16)
  
  ### On this version of the sheet we are using ONLY 2016 YOY FOR THE BAY LEVEL#####
  
  gen_grp <-c(gen_oce,
              gen_bay)
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
  sumlocstats <-pivlocstats%>% group_by(GRP,variable)%>%dplyr::summarize(meanvar=mean(value), sdvar=sd(value),.groups="drop")%>%
    mutate(iter=j, var=sdvar^2)
  
  grp_sumlocstats[[j]] <-sumlocstats
  grp_locstats[[j]] <-loc_stats
}
########### END RAREIFICATION LOOP 10 #########################################################
#our output from this loop is a list of dataframes. Now we need to combine the dataframes by averaging the averages. 

grp_sumlocstats10 <-bind_rows(grp_sumlocstats)
grp_locstats10 <-bind_rows(grp_locstats)

sumlocstats_mean10 <-grp_sumlocstats10 %>% group_by(GRP, variable) %>% 
  dplyr::summarize(av_meanvar = mean(meanvar), av_sd =sqrt(sum(var)/n_distinct(grp_sumlocstats$iter)), .groups="drop")

locstats_mean10 <-grp_locstats10%>%group_by(GRP,LOCUS)%>% 
  dplyr::summarize(N_ALLELES=mean(N_ALLELES),EVENNESS=mean(EVENNESS), Ho=mean(Ho), Hs=mean(Hs),Ht=mean(Ht),Fis=mean(Fis),
            SHANNON_IDX=mean(SHANNON_IDX), SIMPSON_IDX=mean(SIMPSON_IDX), STODD_TAYLOR_IDX=mean(STODD_TAYLOR_IDX),.groups="drop")

########### RAREIFICATION LOOP 5 ###########################

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
  new.shin <-sample_n(df.split$Shin,5)
  new.mt <-sample_n(df.split$Mt,5)
  rare.bays <-bind_rows(new.shin,new.mt,df.split$Jam,df.split$Mor,df.split$Nap)
  #wfpopLD <-df2genind(rare.bays, pop=pop, ploidy=2)
  wf.rare16 <-df2gtypes(rare.bays,ploidy=2)
  wf.rare16 <-gtypes2genind(wf.rare16)
  #setPop(wf.rare16)<-~Bay
  gen_bay <-seppop(wf.rare16)
  
  ### On this version of the sheet we are using ONLY 2016 YOY FOR THE BAY LEVEL#####
  
  gen_grp <-c(gen_oce,
              gen_bay)
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
  sumlocstats <-pivlocstats%>% group_by(GRP,variable)%>%dplyr::summarize(meanvar=mean(value), sdvar=sd(value),.groups="drop")%>%
    mutate(iter=j, var=sdvar^2)
  
  grp_sumlocstats[[j]] <-sumlocstats
  grp_locstats[[j]] <-loc_stats
}
########### END RAREIFICATION LOOP 5 #########################################################
#our output from this loop is a list of dataframes. Now we need to combine the dataframes by averaging the averages. 

grp_sumlocstats5 <-bind_rows(grp_sumlocstats)
grp_locstats5 <-bind_rows(grp_locstats)

sumlocstats_mean5 <-grp_sumlocstats5 %>% group_by(GRP, variable) %>% 
  dplyr::summarize(av_meanvar = mean(meanvar), av_sd =sqrt(sum(var)/n_distinct(grp_sumlocstats$iter)), .groups="drop")

locstats_mean5 <-grp_locstats5%>%group_by(GRP,LOCUS)%>% 
  dplyr::dplyr::summarize(N_ALLELES=mean(N_ALLELES),EVENNESS=mean(EVENNESS), Ho=mean(Ho), Hs=mean(Hs),Ht=mean(Ht),Fis=mean(Fis),
            SHANNON_IDX=mean(SHANNON_IDX), SIMPSON_IDX=mean(SIMPSON_IDX), STODD_TAYLOR_IDX=mean(STODD_TAYLOR_IDX),.groups="drop")


####################### COMBINE ##########################

locstats_mean30 <-mutate(locstats_mean30,nsamp=30)
locstats_mean25 <-mutate(locstats_mean25,nsamp=25)
locstats_mean20 <-mutate(locstats_mean20,nsamp=20)
locstats_mean15 <-mutate(locstats_mean15,nsamp=15)
locstats_mean10 <-mutate(locstats_mean10,nsamp=10)
locstats_mean5 <-mutate(locstats_mean4,nsamp=5)

allloc <-bind_rows(locstats_mean30,locstats_mean25,locstats_mean20,
                   locstats_mean15,locstats_mean10,locstats_mean5)
allocpiv <-pivot_longer(allloc,cols=3:11,names_to = "stat")%>%#mutate(nsamp=as.factor(nsamp))%>%
  filter(nsamp !="NA")

#plot
allocpiv%>%
  filter(stat=="EVENNESS")%>%
ggplot(aes(x = nsamp, y = value, color=GRP)) +
  geom_point() +
  facet_grid(~LOCUS, scales="free")+
  theme_classic()

sumalloc <- allocpiv %>%group_by(GRP,nsamp,stat)%>%dplyr::summarize(value=mean(value),.groups="drop")

sumalloc%>%
  filter(GRP %in% c("Shin","Mt"))%>%
  ggplot(aes(x = nsamp, y = value, color=GRP)) +
  geom_line() + geom_point()+
  facet_wrap(~stat, scales="free", ncol=4)+
  theme_classic()


########## YOU FORGOT TO RARIFY THE OTHER BAYS YA DUM DUM!
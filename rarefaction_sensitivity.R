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

########### RAREIFICATION DOUBLE LOOP ###########################

se <-seq(5,25,5) #input
grp_grp_sumlocstats <-list()
grp_grp_locstats <-list()
#outer loop
for (i in 1:length(se)){
grp_sumlocstats <-list()
grp_locstats <-list()
gen_grp <-list()
#inner loop
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
  new.shin <-sample_n(df.split$Shin,se[i])
  new.mt <-sample_n(df.split$Mt,se[i])
  new.jam <-sample_n(df.split$Jam,se[i])
  new.mor <-sample_n(df.split$Mor,se[i])
  new.nap <-sample_n(df.split$Nap,se[i])
  rare.bays <-bind_rows(new.shin,new.mt, new.mor,new.jam, new.nap)
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

grp_sumlocstats <-bind_rows(grp_sumlocstats)%>%mutate(nsamp=se[i])
grp_locstats <-bind_rows(grp_locstats)%>%mutate(nsamp=se[i])

grp_grp_sumlocstats[[i]] <-grp_sumlocstats
grp_grp_locstats[[i]] <-grp_locstats
}

########### END RAREIFICATION DOUBLE LOOP #########################################################

#bind the list into a dataframe
grp_sumlocstats <-bind_rows(grp_grp_sumlocstats)
grp_locstats <-bind_rows(grp_grp_locstats)

#downscale by averaging the averages. 
sumlocstats_mean <-grp_sumlocstats %>% group_by(GRP, variable, nsamp) %>% 
  dplyr::summarize(av_meanvar = mean(meanvar), av_sd =sqrt(sum(var)/n_distinct(grp_sumlocstats$iter)), .groups="drop")

locstats_mean <-grp_locstats%>%group_by(GRP,LOCUS,nsamp)%>% 
  dplyr::summarize(N_ALLELES=mean(N_ALLELES),EVENNESS=mean(EVENNESS), Ho=mean(Ho), Hs=mean(Hs),Ht=mean(Ht),Fis=mean(Fis),
            SHANNON_IDX=mean(SHANNON_IDX), SIMPSON_IDX=mean(SIMPSON_IDX), STODD_TAYLOR_IDX=mean(STODD_TAYLOR_IDX),.groups="drop")

locstats_mean_piv <-pivot_longer(locstats_mean,cols=4:12,names_to = "stat")%>%#mutate(nsamp=as.factor(nsamp))%>%
  filter(nsamp !="NA")

#plot
locstats_mean_piv%>%
  filter(stat=="EVENNESS"& GRP != "Atl" & GRP != "ALL")%>%
ggplot(aes(x = nsamp, y = value, color=GRP)) +
  geom_point() + geom_line()+
  facet_wrap(~LOCUS, scales="free", nrow=3)+
  theme_classic()

sumalloc <- locstats_mean_piv %>%group_by(GRP,nsamp,stat)%>%dplyr::summarize(value=mean(value),.groups="drop")

sumalloc%>%
  #filter(GRP %in% c("Shin","Mt"))%>%
  filter(GRP != "ALL" & GRP != "Atl")%>%
  ggplot(aes(x = nsamp, y = value, color=GRP)) +
  geom_line() + geom_point()+
  facet_wrap(~stat, scales="free", ncol=4)+
  theme_classic()


########## RAREIFY SHINNECOCK AND MT ONLY ###########################

se <-data.table(se.s=seq(25,50,5),se.m=append(seq(25,35,5),rep(39,3)))

grp_grp_sumlocstats <-list()
grp_grp_locstats <-list()
#outer loop
for (i in 1:length(se$se.s)){
  grp_sumlocstats <-list()
  grp_locstats <-list()
  gen_grp <-list()
  #inner loop
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
    new.shin <-sample_n(df.split$Shin,se[[i,1]])
    new.mt <-sample_n(df.split$Mt,se[[i,2]])
    rare.bays <-bind_rows(new.shin,new.mt, df.split$Mor,df.split$Jam, df.split$Nap)
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
  
  grp_sumlocstats <-bind_rows(grp_sumlocstats)%>%mutate(nsamp=se[[i,1]])
  grp_locstats <-bind_rows(grp_locstats)%>%mutate(nsamp=se[[i,1]])
  
  grp_grp_sumlocstats[[i]] <-grp_sumlocstats
  grp_grp_locstats[[i]] <-grp_locstats
}
################### END SHIN/MT RAREFACTION DOUBLE LOOP ################

#bind the list into a dataframe
grp_sumlocstats <-bind_rows(grp_grp_sumlocstats)
grp_locstats <-bind_rows(grp_grp_locstats)

#downscale by averaging the averages. 
sumlocstats_mean <-grp_sumlocstats %>% group_by(GRP, variable, nsamp) %>% 
  dplyr::summarize(av_meanvar = mean(meanvar), av_sd =sqrt(sum(var)/n_distinct(grp_sumlocstats$iter)), .groups="drop")

locstats_mean <-grp_locstats%>%group_by(GRP,LOCUS,nsamp)%>% 
  dplyr::summarize(N_ALLELES=mean(N_ALLELES),EVENNESS=mean(EVENNESS), Ho=mean(Ho), Hs=mean(Hs),Ht=mean(Ht),Fis=mean(Fis),
                   SHANNON_IDX=mean(SHANNON_IDX), SIMPSON_IDX=mean(SIMPSON_IDX), STODD_TAYLOR_IDX=mean(STODD_TAYLOR_IDX),.groups="drop")

locstats_mean_piv <-pivot_longer(locstats_mean,cols=4:12,names_to = "stat")%>%#mutate(nsamp=as.factor(nsamp))%>%
  filter(nsamp !="NA")

#plot
locstats_mean_piv%>%
  filter(GRP %in% c("Mt","Shin"))%>%
  filter(stat=="EVENNESS"& GRP != "Atl" & GRP != "ALL")%>%
  ggplot(aes(x = nsamp, y = value, color=GRP)) +
  geom_point() + geom_line()+
  facet_wrap(~LOCUS, scales="free", nrow=3)+
  theme_classic()

sumalloc <- locstats_mean_piv %>%group_by(GRP,nsamp,stat)%>%dplyr::summarize(value=mean(value),.groups="drop")

sumalloc%>%
  #filter(GRP %in% c("Shin","Mt"))%>%
  filter(GRP != "ALL" & GRP != "Atl")%>%
  ggplot(aes(x = nsamp, y = value, color=GRP)) +
  geom_line() + geom_point()+
  facet_wrap(~stat, scales="free", ncol=4)+
  theme_classic()


################# RAREFACTION SENSITIVITY FST YOY 2016 ONLY ##########################

setPop(wfpopLD)<-~Bay
df <-genind2df(wfpopLD, usepop = TRUE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(ncol(df),1:(ncol(df)-1))]
df[df=="NA"] <- 0 

#randomly sample from within shinnecock and mattituck 100x

se <-data.table(se.s=seq(25,50,5),se.m=append(seq(25,35,5),rep(39,3)))
grp_FST <-list()
#outer loop
for (j in 1:length(se$se.s)){
rare.out <-list()
for(i in 1:100){
  df.split <-split(df, df$pop)
  new.shin <-sample_n(df.split$Shin,se[[j,1]])
  new.mt <-sample_n(df.split$Mt,se[[j,2]])
  rare.bays <-bind_rows(new.shin,new.mt,df.split$Jam,df.split$Mor,df.split$Nap)
  wf.rare <-df2gtypes(rare.bays,ploidy=2)
  popStruct.b <- popStructTest(wf.rare, nrep = 1000, quietly = TRUE) 
  mtres <- as.data.frame(popStruct.b$pairwise$result)
  mtres <-dplyr::select(mtres, strata.1,strata.2, Fst, Fst.p.val) %>% dplyr::rename(Pr.exact=Fst.p.val)%>%
    mutate(count=i)
  rare.out[[i]]<-mtres
}
grp_FST[[j]] <-bind_rows(rare.out)%>%mutate(nsamp=se[[j,1]])
}
## end FST calc loop 

#to pool p values and Fst, combine across pairwise comparisons. 
mean.fst <-bind_rows(grp_FST)%>%group_by(strata.1,strata.2,nsamp)%>%
  dplyr::summarize(mean.fst=mean(Fst),.groups="drop")
#pool pvals using simes
rare.out2 <-bind_rows(rare.out)%>%unite(pairs,strata.1, strata.2)
rare.out2 <-split(rare.out2,rare.out2$pairs)

simes_pvals <-c()
for (i in 1:length(rare.out2)){
  simes_p <-metapod::combineParallelPValues(as.list(rare.out2[[i]]$Pr.exact),method="simes")
  simes_pvals[i]<-simes_p$p.value
}
mean.fst$comparison <-names(rare.out2); mean.fst$p.adj <-simes_pvals
mean.fst <-mutate(mean.fst, p.adj.BH=p.adjust(p.adj, method="BH"))

mean.fst%>%
  #filter(GRP %in% c("Shin","Mt"))%>%
  ggplot(aes(x = nsamp, y = mean.fst)) +
  geom_line() + geom_point()+
  facet_grid(strata.1~strata.2)

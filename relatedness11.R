### Relatedness #####

## NOTES ##
## includes only 12 loci for now.
## for no half 0 genotypes, use the double0 version of the files 
## population corrected for new cohort assignment
## based on the script "skreportdnaJAN_28_2020MATTITUCK.R"
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

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")

##### Formating the dataset #####
wfpop <- read.genalex("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrected12satsAug1920204genalex_doubl0.csv")
wfpop4df <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrected_12_satsAug2020_doubl0.csv", header = TRUE) #csv version 

#create a second dataset that removes one loci in the WF27/WF33 pair 
#go back to wfpop for this. 
wfpopLD <-genclone2genind(wfpop)
all_loci <- locNames(wfpopLD)# create vector of all loci
removeloc <- c("WF27")# create vector containing loci to remove
keeploc <- setdiff(all_loci, removeloc)# create list of loci to keep
wfpopLD <- wfpopLD[loc = keeploc]# filter loci in genind object
length(locNames(wfpopLD))# check number of loci in genind obj

#now re-remove the missing data. 
info_table(wfpopLD, plot = TRUE, scaled =FALSE)
wfpopLD <- wfpopLD %>% missingno("geno", cutoff=0.20) # remove samples that are > 20% missing
#ggsave("20percutinfotable11LD.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()


######## Pairwise Relatedness ######

#There are many different estimators for relatedness. In previous scripts we have tried more robust ways. 
# for now, Shannon's way. https://gist.github.com/sjoleary/3efd4a7d56b115fad319781298765a31

#first convert to the right data format
setPop(wfpopLD) <-~Bay
df <-genind2df(wfpopLD, usepop = FALSE, oneColPerAll = TRUE) 
df$Ind <- rownames(df)
df <-df[,c(23,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)]
df[df=="NA"] <- 0 # missing data must be 0
write.table(df, "scratch",row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) #write it.
genotypedata <- readgenotypedata("scratch")# import input file as list (gdata, nloci, nalleles, ninds, freqs)

# pairwise relatedness (Lynch & Ritland 1999)
relatedness_lynchrd <- coancestry(genotypedata$gdata,
                                  lynchrd = 1)

# write relatedness to file
relatedn <- relatedness_lynchrd$relatedness %>%
  dplyr::select(pair.no, ind1.id, ind2.id, lynchrd)

#not sure what's going on here. 
library("readr")
write_delim(relatedn, "pairwise_relatedness")
# write inbreeding to file
inbreed <- relatedness_lynchrd$inbreeding %>%
  dplyr::select(ind.id, L3, LH) %>%
  dplyr::rename(INDV = ind.id)
write_delim(inbreed, "inbreeding")

#the simulation
#change bootstrap back to 100 because 1000 just takes too long. 
cosim <-compareestimators(relatedness_lynchrd, ninds=100)
cosim
#ggsave("relatendess_simulation.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()

#write_delim(cosim, "cosim") #wrote this after doing 1000 runs

#extract mean, median & CI values for relatedness from the simulation. 
relsim <- as.data.frame(cosim$data)
relply<- ddply(relsim, Estimator~relationship,summarize, mean.rel = mean(Relatedness_Value), LIrel = quantile(Relatedness_Value,0.05), HIrel = quantile(Relatedness_Value, 0.95), medianrel= quantile(Relatedness_Value, 0.5))
relply

# relatedness - remember to change these values depending on what dataset you are using. 
ggplot(relatedn, aes(x = lynchrd)) +
  geom_histogram(binwidth = 0.001, color = "black", fill = "darkorange") +
  geom_vline(aes(xintercept = mean(lynchrd, na.rm = TRUE)),
             color = "darkblue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = quantile(lynchrd, 0.95)),
             color = "darkred", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = 0.23), color = "green", linetype = "dashed", size = 1) + #the estimated mean value for half sibs in the simulation
  geom_vline(aes(xintercept = 0.48), color = "green", linetype = "dashed", size = 1) + #the estimated mean value for full sibs in the simulation
  geom_vline(aes(xintercept = 0.47), color = "green", linetype = "dashed", size = 1) + #the estimated mean value for parent offspring in the simulation
  geom_vline(aes(xintercept = 0.03),color = "purple", linetype = "dashed", size = 0.5) + #the Lynchrd 0.05 quantile on the simulated half sibs relationship
  geom_vline(aes(xintercept = 0.53), color = "darkred", linetype = "dashed", size = 1) + #red is the 0.95 quantile
  labs(x = "relatedness", y = "number of pairs")+
  theme_cowplot()
#ggsave("pairwiserelatenessALL_11.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()

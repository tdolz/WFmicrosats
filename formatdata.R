#format data for microdrop. 

###September 19, 2020

#install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)
library('plyr')
library("dplyr")
library('readr')

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")
wfpop4df <-read.csv("/Users//tdolan/Documents//R-Github//WFmicrosats/popcorrect_20_sept2020_doubl0.csv", header = TRUE) #csv version 

allele1 <- dplyr::select(wfpop4df, Ind,Pop, "J42.1", "WF22.1","PAM27.1","PAM79.1","WF27.1","WF16.1", "WF33.1","WF3.1","WF517.1","WF196.1","WF223.1","WF421.1",
                      "A441.1","PAM21.1", "Psy087.1", "Psy022.1", "WF06.1","WF12.1","WF01.1","WF32.1")

names(allele1) <-c("Ind","Pop","J42", "WF22","PAM27","PAM79","WF27","WF16", "WF33","WF3","WF517","WF196","WF223","WF421",
                   "A441","PAM21", "Psy087", "Psy022", "WF06","WF12","WF01","WF32")


allele2 <-dplyr::select(wfpop4df, Ind,Pop, "J42.2", "WF22.2","PAM27.2","PAM79.2","WF27.2","WF16.2", "WF33.2","WF3.2","WF517.2","WF196.2","WF223.2","WF421.2",
                         "A441.2","PAM21.2", "Psy087.2", "Psy022.2", "WF06.2","WF12.2","WF01.2","WF32.2")

names(allele2) <-c("Ind","Pop","J42", "WF22","PAM27","PAM79","WF27","WF16", "WF33","WF3","WF517","WF196","WF223","WF421",
                   "A441","PAM21", "Psy087", "Psy022", "WF06","WF12","WF01","WF32")

microdrop <-full_join(allele1,allele2)
microdrop <-arrange(microdrop, Ind)
microdrop[microdrop=="0"] <- -9 # missing data must be -9

write_delim(microdrop, "microdrop_datafile")

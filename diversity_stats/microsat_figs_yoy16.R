#Microsats YOY2016 plots 
#created 1/3/24
#
#install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)
library("tidyverse")
library("poppr")
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
library("data.table")

#################### giant barplot 1 #####################################################
##### FINAL GIANT BARPLOT ########
pivlocstats <-read.csv("./diversity_stats/diversity_output_files/YOY16_bay/pivlocstats17.csv", header=T)
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
ggsave('pivlocstats17.png', path="./diversity_stats/diversity_figs/YOY16", width = 12, height = 5)



################################## BAY diversity barplots #############################################################
################################## 
meltlocstats_bay <- read.csv("./diversity_stats/diversity_output_files/YOY16_bay/meltlocstats_bay.csv",header=T)

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
ggsave('nei_bay17.png', path="./diversity_stats/diversity_figs/YOY16", width = 10, height = 5)

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
ggsave('FIS_bay17.png', path="./diversity_stats/diversity_figs/YOY16", width = 10, height = 5)

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
ggsave('Evenness_bay17.png', path="./diversity_stats/diversity_figs/YOY16", width = 10, height = 5)

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
ggsave('Shannon_idxbay17.png', path="./diversity_stats/diversity_figs/YOY16", width = 10, height = 5)

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
ggsave('Hsbay17.png', path="./diversity_stats/diversity_figs/YOY16", width = 10, height = 5)

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
ggsave('Hobay17.png', path="./diversity_stats/diversity_figs/YOY16", width = 10, height = 5)

#rareiied allelic richness
meltar3 <- read.csv(file="./diversity_stats/diversity_output_files/YOY16_bay/meltar3.csv",header=T)
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
ggsave('rareifiedallelesLD_bay17.png',path="./diversity_stats/diversity_figs/YOY16", width = 10, height = 5)

### private alleles ##
pA <-read.csv("./diversity_stats/diversity_output_files/YOY16_bay/private_allelecount.csv",header=T)
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
ggsave('privateallelesLD_bay17.png',path="./diversity_stats/diversity_figs/YOY16", width = 10, height = 5)

## internal relatedness ##
rel2 <-read.csv("./diversity_stats/diversity_output_files/YOY16_bay/IR.csv",header=T)
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
ggsave('rel_bay17.png', path="./diversity_stats/diversity_figs/YOY16", width = 10, height = 5)



############################################ MATTITUCK DIVERSITY BARPLOTS NO YOY 2015 ###################################
############################################ 

#allelic richness
ar_mtshi <-read.csv("./diversity_stats/diversity_output_files/YOY16_bay/rarefied_allelecount_ALLGROUPS.csv",header=T)
meltarMT <-pivot_longer(ar_mtshi, cols=c("Mt_2016_1","Mt_2016_2","Mt_adults_3", "Mt_adults_4"),names_to="variable", values_to="value")

## private alleles
pA_mtshi <- read.csv("./diversity_stats/diversity_output_files/YOY16_bay/rarefied_allelecount_ALLGROUPS.csv",header=T)
meltPA_MT <-pivot_longer(pA_mtshi, cols=c("Mt_2016_1","Mt_2016_2","Mt_adults_3", "Mt_adults_4"),names_to="variable", values_to="value")

#########################################SHINNECOCK DIVERSITY BARPLOTS ###########################################################

#rareified allelic richness
meltarSHIN <-pivot_longer(ar_mtshi, cols=c("Shin_2016_1", "Shin_2017_1", "Shin_2016_2", "Shin_2017_2"),names_to="variable", values_to="value")

#private alleles 
meltPA_SHIN <-pivot_longer(pA_mtshi, cols=c("Shin_2016_1", "Shin_2017_1", "Shin_2016_2", "Shin_2017_2"),names_to="variable", values_to="value")

#internal relatedness

################################# WILCOXON TESTS HEATMAPS BAY #######################################################
################################# 

results <-read.csv(file="./diversity_stats/diversity_output_files/YOY16_bay/wilcox.csv",header=T)
results <- mutate(results, starsig=ifelse(significance=="significant",round(p.value,4),NA),test.statistic=abs(stat))
results<-mutate(results, p_value = cut(p.value, breaks=c(0,0.005, 0.01,0.05,0.1,0.5,1)))

cols <- c("(0,0.005]"="#034e7b", "(0.005,0.01]" = "#045a8d", "(0.01,0.05]" = "#2b8cbe", "(0.05,0.1]" = "#74a9cf", "(0.1,0.5]"  = "#a6bddb", "(0.5,1]"="#d0d1e6")
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
ggsave('neiBaytile17.png',path="./diversity_stats/diversity_figs/YOY16", width = 10, height = 7)

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
ggsave('FisBaytile17.png',path="./diversity_stats/diversity_figs/YOY16", width = 10, height = 7)

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
ggsave('evennessBaytile17.png',path="./diversity_stats/diversity_figs/YOY16", width = 10, height = 7)

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
ggsave('shannonBaytile17.png',path="./diversity_stats/diversity_figs/YOY16", width = 10, height = 7)

# allelic richness ====
results_ar %>%
  ggplot(aes(x = pair1, y = pair2))+
  geom_tile(aes(fill=p_value), show.legend = TRUE)+ #filled  by p_value (discrete)
  scale_fill_manual(values=cols)+
  geom_text(aes(label = round(p.value,3),color=significance), size=5)+
  scale_color_manual(values=c("white","red"), guide=FALSE)+
  xlab("")+ylab("")+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "white", colour = "black"))
ggsave('rariefied_allelesBaytile17.png',path="./diversity_stats/diversity_figs/YOY16", width = 10, height = 7)

### INTERNAL RELATEDNESS
t <-read.csv("./diversity_stats/diversity_output_files/YOY16_bay/IR_ttests.csv",header=T)
cols <- c("(0,0.005]"="#034e7b", "(0.005,0.01]" = "#045a8d", "(0.01,0.05]" = "#2b8cbe", "(0.05,0.1]" = "#74a9cf", "(0.1,0.5]"  = "#a6bddb", "(0.5,1]"="#d0d1e6")

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
ggsave('IRBaytile17.png',path="./diversity_stats/diversity_figs/YOY16", width = 10, height = 7)



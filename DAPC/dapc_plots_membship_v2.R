### DAPC graphs ###
### march 14, 2024
### attempting to make the membership plots prettier. 

library('plyr')
library("poppr")
library("tidyr")
library('purrr')
library("ggplot2")
library("adegenet")
library("pegas")
library("cowplot")
library("vcfR")
library("adegenet")
library("forcats")
library("dplyr")

################ LI BAYS ####################
grp_membership <-read.csv("./DAPC/DAPC_figs/yoy16dapc.csv")

baycolors <-c("#332288","#88CCEE","#44AA99","#999933","#DDCC77","#CC6677","#882255","#AA4499")
drabcolors <-c("#d0d1e6","#a6bddb", "#67a9cf", "#1c9099", "#016450")

grp_membership <- mutate(grp_membership, groupnum=as.character(groupnum), Ind=as.factor(Ind))

row1 <-filter(grp_membership, nclust==2)%>%arrange(pop,major.group,Ind)

row1.x <-row1 %>%arrange(pop,major.group)
level_order <-row1.x$Ind
level_order_df <-as.data.frame(level_order)
level_order_df$index <-seq(1,nrow(level_order_df),1)
row1$level_order <-level_order_df$index

test_plot <-row1 %>% ggplot(aes(x=fct_reorder(Ind,level_order), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  #scale_x_discrete(limits=level_order)+
  scale_fill_manual(values = baycolors)+
  theme_cowplot(); test_plot
#works!
k2 <-filter(grp_membership, nclust==2)%>%mutate(Ind=as.factor(Ind))
#k2$level_order <-level_order

#for (i in 1:length(grp_membership$Ind)){
#  p <-filter(k2, level_order == grp_membership$Ind[i])%>%
 #   arrange(major.group)
 # grp_membership$level_order[i] <-p$level_order[1]
#}

grp_membership <-grp_membership %>% arrange(nclust,level_order)
#grp_membership$level_order <-level_order_df$index

#row1.x <-grp_membership %>%arrange(nclust,pop,major.group,Ind)
#level_order <-row1.x$Ind
#level_order_df <-as.data.frame(level_order)
#level_order_df$index <-seq(1,nrow(level_order_df),1)
#grp_membership$level_order <-level_order_df$index
#grp_membership <-mutate(grp_membership, groupnum=as.factor(groupnum), nclust=as.factor(nclust),level_order=as.numeric(level_order))

kmeans_plot <- grp_membership%>%
  filter(nclust !="ap")%>%
  arrange(nclust,pop,MEMBSHIP)%>%
  #ggplot(aes(x = fct_reorder(Ind,level_order), y = MEMBSHIP, fill = groupnum)) +
  ggplot(aes(x = fct_reorder(Ind,MEMBSHIP), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  #scale_x_discrete(limits=level_order)+
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(
    axis.text.x = element_blank(),
    axis.title.x=element_blank(),
    #axis.text.y = element_text(size=9),
    strip.background =element_rect(fill="transparent"),
    #strip.text = element_text(size=9),
    #legend.position = "none",
    legend.position = "right",
    panel.spacing = unit(1.5, "mm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = "dark grey"),
    plot.background = element_rect(fill = NA,colour = "black", size=0.5),
    panel.border = element_rect(colour="black",fill=NA, size=0.3)
  ); kmeans_plot

ap_plot <- grp_membership%>%
  filter(nclust =="ap")%>%
  arrange(major.group,Ind,MEMBSHIP)%>%
  ggplot(aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = drabcolors)+
  theme_cowplot()+
  theme(
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size=9),
    #axis.title.y=element_blank(),
    strip.background =element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    #legend.position = "right",
    panel.spacing = unit(1.5, "mm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = "dark grey"),
    plot.background = element_rect(fill = NA,colour = "black", size=0.5),
    panel.border = element_rect(colour="black",fill=NA, size=0.3)
  ); ap_plot

plot_grid(kmeans_plot, ap_plot, ncol=1, rel_heights = c(1, 1/4))
ggsave("bays16_membship_facet.png", width=8, height=10, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")


#### plot grid separately. but can we order it by the first one? 
grp_membership <-read.csv("./DAPC/DAPC_figs/yoy16dapc.csv")
grp_membership <- mutate(grp_membership, groupnum=as.character(groupnum), Ind=as.factor(Ind))

#k2
k2 <-filter(grp_membership, nclust==2)%>%arrange(pop,major.group,Ind)
k2.x <-k2 %>%arrange(pop,major.group)
level_order <-k2.x$Ind
level_order_df <-as.data.frame(level_order)
names(level_order_df)<-c("Ind")
level_order_df$index <-seq(1,nrow(level_order_df),1)
k2$level_order <-level_order_df$index

k2p <-k2 %>% ggplot(aes(x=fct_reorder(Ind,level_order), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        #axis.text.y = element_text(size=9),
        strip.background =element_rect(fill="transparent"),
        #strip.text = element_text(size=9),
        legend.position = "none",
        panel.spacing = unit(1.5, "mm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "dark grey"),
        plot.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=0.3)); k2p
#works!

# just do the panel plot. 
#list the middle 5 plots with a for loop
plot_list <-list()
for (i in 3:8){
  p <- grp_membership%>%
    filter(nclust==i)%>%arrange(pop,major.group,Ind)%>%
    ggplot(aes(x=fct_reorder(Ind,major.group), y = MEMBSHIP, fill = groupnum)) +
    geom_bar(stat = "identity") +
    labs(x = "", y = "") +
    facet_grid(nclust ~ pop, scales = "free", space = "free") +
    scale_fill_manual(values = baycolors)+
    theme_cowplot()+
    theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
          strip.background =element_blank(),strip.text.x = element_blank(),
          legend.position = "none", panel.spacing = unit(1.5, "mm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = "dark grey"),
          plot.background = element_blank(),
          panel.border = element_rect(colour="black",fill=NA, size=0.3))
  plot_list[[i]]<-p
}

#kap
kap <- grp_membership%>%
  filter(nclust=="ap")%>%arrange(pop,major.group,Ind)%>%
  ggplot(aes(x=fct_reorder(Ind,MEMBSHIP), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = drabcolors)+
  theme_cowplot()+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        strip.background =element_blank(),strip.text.x = element_blank(),
        legend.position = "none", panel.spacing = unit(1.5, "mm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "dark grey"),
        plot.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=0.3)); kap

plot_grid(k2p,plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],kap, ncol=1)
ggsave("bays16_membship_plotgrid_selforder.png", width=8, height=10, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

###################### MATTITUCK CREEK (including 2015) ########################
#################################################################################################################################################################
################################################################################
#### plot grid separately. but can we order it by the first one? 
grp_membership <-read.csv("./DAPC/DAPC_figs/mt2015_dapc.csv")
grp_membership <- mutate(grp_membership, groupnum=as.character(groupnum), Ind=as.factor(Ind))

#k2
k2 <-filter(grp_membership, nclust==2)%>%arrange(pop,major.group,Ind)
k2.x <-k2 %>%arrange(pop,major.group)
level_order <-k2.x$Ind
level_order_df <-as.data.frame(level_order)
names(level_order_df)<-c("Ind")
level_order_df$index <-seq(1,nrow(level_order_df),1)
k2$level_order <-level_order_df$index

k2p <-k2 %>% ggplot(aes(x=fct_reorder(Ind,level_order), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        #axis.text.y = element_text(size=9),
        strip.background =element_rect(fill="transparent"),
        #strip.text = element_text(size=9),
        legend.position = "none",
        panel.spacing = unit(1.5, "mm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "dark grey"),
        plot.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=0.3)); k2p
#works!

# just do the panel plot. 
k3p <- grp_membership%>%
  filter(nclust==3)%>%arrange(pop,major.group,Ind)%>%
  ggplot(aes(x=fct_reorder(Ind,major.group), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        strip.background =element_blank(),strip.text.x = element_blank(),
        legend.position = "none", panel.spacing = unit(1.5, "mm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "dark grey"),
        plot.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=0.3)); k3p

# k4 
k4p <- grp_membership%>%
  filter(nclust==4)%>%arrange(pop,major.group,Ind)%>%
  ggplot(aes(x=fct_reorder(Ind,major.group), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        strip.background =element_blank(),strip.text.x = element_blank(),
        legend.position = "none", panel.spacing = unit(1.5, "mm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "dark grey"),
        plot.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=0.3)); k4p

#k5
k5p <- grp_membership%>%
  filter(nclust==5)%>%arrange(pop,major.group,Ind)%>%
  ggplot(aes(x=fct_reorder(Ind,major.group), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        strip.background =element_blank(),strip.text.x = element_blank(),
        legend.position = "none", panel.spacing = unit(1.5, "mm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "dark grey"),
        plot.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=0.3)); k5p

# k6
k6p <- grp_membership%>%
  filter(nclust==6)%>%arrange(pop,major.group,Ind)%>%
  ggplot(aes(x=fct_reorder(Ind,major.group), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        strip.background =element_blank(),strip.text.x = element_blank(),
        legend.position = "none", panel.spacing = unit(1.5, "mm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "dark grey"),
        plot.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=0.3)); k6p

# k7 
k7p <- grp_membership%>%
  filter(nclust==7)%>%arrange(pop,MEMBSHIP,Ind)%>%
  ggplot(aes(x=fct_reorder(Ind,major.group), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        strip.background =element_blank(),strip.text.x = element_blank(),
        legend.position = "none", panel.spacing = unit(1.5, "mm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "dark grey"),
        plot.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=0.3)); k7p

# k8. 
k8p <- grp_membership%>%
  filter(nclust==8)%>%arrange(pop,MEMBSHIP,Ind)%>%
  ggplot(aes(x=fct_reorder(Ind,major.group), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        strip.background =element_blank(),strip.text.x = element_blank(),
        legend.position = "none", panel.spacing = unit(1.5, "mm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "dark grey"),
        plot.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=0.3)); k8p

#kap
kap <- grp_membership%>%
  filter(nclust=="ap")%>%arrange(pop,major.group,Ind)%>%
  ggplot(aes(x=fct_reorder(Ind,MEMBSHIP), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = drabcolors)+
  theme_cowplot()+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        strip.background =element_blank(),strip.text.x = element_blank(),
        legend.position = "none", panel.spacing = unit(1.5, "mm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "dark grey"),
        plot.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=0.3)); kap

plot_grid(k2p,k3p,k4p,k5p,k6p,k7p,k8p,kap, ncol=1)
ggsave("bays16_membship_plotgrid_selforder.png", width=8, height=10, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")


########################################## MATTITUCK 2015 #######################################################################################
grp_membership <-read.csv("./DAPC/DAPC_figs/mt2015_dapc.csv")

baycolors <-c("#332288","#88CCEE","#44AA99","#999933","#DDCC77","#CC6677","#882255","#AA4499")
Mtcolors4 <-c("#fad2e1","#d8bbff","#ff006e","#8338ec","#4361ee","#4cc9f0")

grp_membership <- mutate(grp_membership, groupnum=as.character(groupnum), Ind=as.factor(Ind))

row1 <-filter(grp_membership, nclust==2)%>%arrange(pop,major.group,Ind)

row1.x <-row1 %>%arrange(pop,major.group)
level_order <-row1.x$Ind
level_order_df <-as.data.frame(level_order)
level_order_df$index <-seq(1,nrow(level_order_df),1)
row1$level_order <-level_order_df$index

kmeans_plot <- grp_membership%>%
  filter(nclust !="ap")%>%
  arrange(nclust,pop,MEMBSHIP)%>%
  #ggplot(aes(x = fct_reorder(Ind,level_order), y = MEMBSHIP, fill = groupnum)) +
  ggplot(aes(x = fct_reorder(Ind,MEMBSHIP), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  #scale_x_discrete(limits=level_order)+
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(
    axis.text.x = element_blank(),
    axis.title.x=element_blank(),
    #axis.text.y = element_text(size=9),
    strip.background =element_rect(fill="transparent"),
    #strip.text = element_text(size=9),
    #legend.position = "none",
    legend.position = "right",
    panel.spacing = unit(1.5, "mm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = "dark grey"),
    plot.background = element_rect(fill = NA,colour = "black", size=0.5),
    panel.border = element_rect(colour="black",fill=NA, size=0.3)
  ); kmeans_plot

ap_plot <- grp_membership%>%
  filter(nclust =="ap")%>%
  arrange(major.group,Ind,MEMBSHIP)%>%
  ggplot(aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = Mtcolors4)+
  theme_cowplot()+
  theme(
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size=9),
    #axis.title.y=element_blank(),
    strip.background =element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    #legend.position = "right",
    panel.spacing = unit(1.5, "mm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = "dark grey"),
    plot.background = element_rect(fill = NA,colour = "black", size=0.5),
    panel.border = element_rect(colour="black",fill=NA, size=0.3)
  ); ap_plot

plot_grid(kmeans_plot, ap_plot, ncol=1, rel_heights = c(1, 1/4))
ggsave("mt15_membship_facet.png", width=8, height=10, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

###############################

#k2
k2 <-filter(grp_membership, nclust==2)%>%arrange(pop,major.group,Ind)
k2.x <-k2 %>%arrange(pop,major.group)
level_order <-k2.x$Ind
level_order_df <-as.data.frame(level_order)
names(level_order_df)<-c("Ind")
level_order_df$index <-seq(1,nrow(level_order_df),1)
k2$level_order <-level_order_df$index

k2p <-k2 %>% ggplot(aes(x=fct_reorder(Ind,level_order), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        #axis.text.y = element_text(size=9),
        strip.background =element_rect(fill="transparent"),
        #strip.text = element_text(size=9),
        legend.position = "none",
        panel.spacing = unit(1.5, "mm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "dark grey"),
        plot.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=0.3)); k2p

#list the middle 5 plots 
plot_list <-list()
for (i in 3:8){
  p <- grp_membership%>%
    filter(nclust==i)%>%arrange(pop,major.group,Ind)%>%
    ggplot(aes(x=fct_reorder(Ind,major.group), y = MEMBSHIP, fill = groupnum)) +
    geom_bar(stat = "identity") +
    labs(x = "", y = "") +
    facet_grid(nclust ~ pop, scales = "free", space = "free") +
    scale_fill_manual(values = baycolors)+
    theme_cowplot()+
    theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
          strip.background =element_blank(),strip.text.x = element_blank(),
          legend.position = "none", panel.spacing = unit(1.5, "mm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = "dark grey"),
          plot.background = element_blank(),
          panel.border = element_rect(colour="black",fill=NA, size=0.3))
  plot_list[[i]]<-p
}

#kap
kap <- grp_membership%>%
  filter(nclust=="ap")%>%arrange(pop,major.group,Ind)%>%
  ggplot(aes(x=fct_reorder(Ind,MEMBSHIP), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = Mtcolors4)+
  theme_cowplot()+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        strip.background =element_blank(),strip.text.x = element_blank(),
        legend.position = "none", panel.spacing = unit(1.5, "mm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "dark grey"),
        plot.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=0.3)); kap

plot_grid(k2p,plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],kap, ncol=1)
ggsave("MT15_membship_plotgrid_selforder.png", width=8, height=10, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")


########################################## SHINNECOCK #######################################################################################
grp_membership <-read.csv("./DAPC/DAPC_figs/shin_dapc.csv")

baycolors <-c("#332288","#88CCEE","#44AA99","#999933","#DDCC77","#CC6677","#882255","#AA4499")
shincolors2 <-c("#29bf12","#abff4f","#3f8efc","#3b28cc")

grp_membership <- mutate(grp_membership, groupnum=as.character(groupnum), Ind=as.factor(Ind))

row1 <-filter(grp_membership, nclust==2)%>%arrange(pop,major.group,Ind)

row1.x <-row1 %>%arrange(pop,major.group)
level_order <-row1.x$Ind
level_order_df <-as.data.frame(level_order)
level_order_df$index <-seq(1,nrow(level_order_df),1)
row1$level_order <-level_order_df$index

kmeans_plot <- grp_membership%>%
  filter(nclust !="ap")%>%
  arrange(nclust,pop,MEMBSHIP)%>%
  #ggplot(aes(x = fct_reorder(Ind,level_order), y = MEMBSHIP, fill = groupnum)) +
  ggplot(aes(x = fct_reorder(Ind,MEMBSHIP), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "membership probability") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  #scale_x_discrete(limits=level_order)+
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(
    axis.text.x = element_blank(),
    axis.title.x=element_blank(),
    #axis.text.y = element_text(size=9),
    strip.background =element_rect(fill="transparent"),
    #strip.text = element_text(size=9),
    #legend.position = "none",
    legend.position = "right",
    panel.spacing = unit(1.5, "mm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = "dark grey"),
    plot.background = element_rect(fill = NA,colour = "black", size=0.5),
    panel.border = element_rect(colour="black",fill=NA, size=0.3)
  ); kmeans_plot

ap_plot <- grp_membership%>%
  filter(nclust =="ap")%>%
  arrange(major.group,Ind,MEMBSHIP)%>%
  ggplot(aes(x = fct_reorder(Ind,major.group), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = shincolors2)+
  theme_cowplot()+
  theme(
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size=9),
    #axis.title.y=element_blank(),
    strip.background =element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    #legend.position = "right",
    panel.spacing = unit(1.5, "mm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = "dark grey"),
    plot.background = element_rect(fill = NA,colour = "black", size=0.5),
    panel.border = element_rect(colour="black",fill=NA, size=0.3)
  ); ap_plot

plot_grid(kmeans_plot, ap_plot, ncol=1, rel_heights = c(1, 1/4))
ggsave("shin_membship_facet.png", width=8, height=10, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")

###############################

#k2
k2 <-filter(grp_membership, nclust==2)%>%arrange(pop,major.group,Ind)
k2.x <-k2 %>%arrange(pop,major.group)
level_order <-k2.x$Ind
level_order_df <-as.data.frame(level_order)
names(level_order_df)<-c("Ind")
level_order_df$index <-seq(1,nrow(level_order_df),1)
k2$level_order <-level_order_df$index

k2p <-k2 %>% ggplot(aes(x=fct_reorder(Ind,level_order), y = MEMBSHIP, fill = groupnum)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = baycolors)+
  theme_cowplot()+
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        #axis.text.y = element_text(size=9),
        strip.background =element_rect(fill="transparent"),
        #strip.text = element_text(size=9),
        legend.position = "none",
        panel.spacing = unit(1.5, "mm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "dark grey"),
        plot.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=0.3)); k2p

#list the middle 5 plots 
plot_list <-list()
for (i in 3:8){
  p <- grp_membership%>%
    filter(nclust==i)%>%arrange(pop,major.group,Ind)%>%
    ggplot(aes(x=fct_reorder(Ind,major.group), y = MEMBSHIP, fill = groupnum)) +
    geom_bar(stat = "identity") +
    labs(x = "", y = "") +
    facet_grid(nclust ~ pop, scales = "free", space = "free") +
    scale_fill_manual(values = baycolors)+
    theme_cowplot()+
    theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
          strip.background =element_blank(),strip.text.x = element_blank(),
          legend.position = "none", panel.spacing = unit(1.5, "mm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = "dark grey"),
          plot.background = element_blank(),
          panel.border = element_rect(colour="black",fill=NA, size=0.3))
  plot_list[[i]]<-p
}

#kap
kap <- grp_membership%>%
  filter(nclust=="ap")%>%arrange(pop,major.group,Ind)%>%
  ggplot(aes(x=fct_reorder(Ind,MEMBSHIP), y = MEMBSHIP, fill = GRP)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "") +
  facet_grid(nclust ~ pop, scales = "free", space = "free") +
  scale_fill_manual(values = shincolors2)+
  theme_cowplot()+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        strip.background =element_blank(),strip.text.x = element_blank(),
        legend.position = "none", panel.spacing = unit(1.5, "mm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = "dark grey"),
        plot.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=0.3)); kap

plot_grid(k2p,plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],kap, ncol=1)
ggsave("SHIN_membship_plotgrid_selforder.png", width=8, height=10, path="/Users/tdolan/documents/R-Github/WFmicrosats/DAPC/DAPC_figs")








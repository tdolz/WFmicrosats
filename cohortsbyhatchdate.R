### WF Microsats ###
## Mattituck Cohort assignment by hatch date ###
## August 18, 2020 ###

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github//WFmicrosats")

library("dplyr")
library("tidyr")
library('purrr')
library("lubridate")
library("chron")
library("ggplot2")
library('ggrepel')

mtage <-read.csv(file="mtcohortreassignAug182020.csv", header=TRUE)
mtage <-mutate(mtage, bday=as.Date(hatch.date), cohort=as.factor(new.cohort)) %>%slice(1:49)

mtage %>%
  #ggplot(aes(x=bday, y=n_distinct(samp)))+ #turn off cohort sorting
 ggplot(aes(x=bday, y=cohort,group=cohort))+
 geom_point(aes(color=cohort,size=TL))+
  #geom_point(aes(size=TL))+ #turn off cohort sorting
  #geom_text(aes(label=samp),hjust=0,vjust=0)+ #another way to do text. 
  geom_label_repel(aes(label = samp),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_x_date(date_breaks="1 week",date_labels = ("%Y-%m-%d"),)+
  xlab("settlement (maybe hatch) date")+
  ggtitle("Mattituck Cohort Assignment")+
  facet_grid(~Year, scales="free")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = 'white', colour = 'light gray'),
        panel.grid.major = element_line(colour = "white"))
#ggsave("Mtcohortassign.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()

#Ridges graphs for Mattituck

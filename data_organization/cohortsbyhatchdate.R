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
allmelt <-read.csv(file="allmeltdec2019.csv", header=TRUE)
spawndate <-read.csv(file="mtspawndate.csv", header=TRUE)

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
  xlab("settlement date")+
  ggtitle("Mattituck Cohort Assignment")+
  facet_grid(~Year, scales="free")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = 'white', colour = 'light gray'),
        panel.grid.major = element_line(colour = "white"))
#ggsave("Mtcohortassign.png", path="/Users/tdolan/Documents/WIP research/microsats/microsat_figs")
#dev.off()

#merge spawn date and fish information. 
spawndate <-dplyr::rename(spawndate, fishID=Tag.ID)
mtage2 <-left_join(mtage,spawndate, by="fishID")
#now we have spawning date and settlement date. 
filter(mtage2, cohort=="2") %>% mean(Spawning.Date)

#histograms where catherine's cohorts are overlaid on the regular cohort data. 
## original script is called "field graphs 9_10_19.R"
mtage <-dplyr::select(mtage,-cohort,-old.pop)%>%mutate(source="cc")

mtmelt <-filter(allmelt, Bay=="Mattituck") %>%dplyr::rename(catch.date=Date, Location=Bay, TL=value) %>%mutate(source="tows")
weeks <-dplyr::select(mtmelt,catch.date,Week)
mtage <-left_join(mtage,weeks)
mtfnew <-full_join(mtage, mtmelt, by=c("catch.date","TL","Location","Year","source","Week")) %>% mutate(Year=as.factor(Year),catch.date=as.Date(catch.date))

split <-c(40,	47,	49,	43,	57,	70)
splits15 <- as.data.frame(split) %>% mutate(Year="2015")
splits15$catch.date <- as.Date(c("2015-06-24",	"2015-07-07",	"2015-07-21",	"2015-08-03",	"2015-08-17",	"2015-09-01"))
split <-c(49,	52,	63,	84)
splits16 <- as.data.frame(split) %>% mutate(Year="2016")
splits16$catch.date <- as.Date(c("2016-06-30",	"2016-07-11",	"2016-07-26",	"2016-08-08"))
splits <-bind_rows(splits15,splits16) %>% mutate(Year=as.factor(Year),catch.date=as.Date(catch.date))

mtfnew <-full_join(mtfnew, splits, by=c("Year","catch.date"))
mtfnew <-unique(mtfnew) %>% mutate(day = format(as.Date(catch.date), format="%m/%d"))

#both years in same plot
ggplot(mtfnew,aes(x=TL)) + 
  geom_histogram(data = subset(mtfnew,source=="cc"& new.cohort=="1"), fill="#e66101",color="black", position="identity",   show.legend=TRUE)+ 
  geom_histogram(data = subset(mtfnew,source=="cc"& new.cohort=="2"), fill="#018571",color="black",position="identity",   show.legend=TRUE)+
  geom_histogram(data = subset(mtfnew,source=="tows"), fill="#5e3c99",color="#5e3c99",position="identity",   alpha = 0.2, show.legend=TRUE)+
  geom_vline(aes(xintercept=split), color="black", linetype="dashed", size=0.5) +
  facet_grid(day~Year,scales="free_y")+
  theme_classic()+
  labs(x = "length (mm)",
       title = "YOY cohorts by length",
       subtitle = "Mattituck") 

#separate plots for each year - Remember to change out year!
mtfnew2 <-filter(mtfnew, Year=="2015")
ggplot(mtfnew2,aes(x=TL)) + 
  geom_histogram(data = subset(mtfnew2,source=="cc"& new.cohort=="1"), fill="#e66101",color="black", position="identity",   show.legend=TRUE)+ 
  geom_histogram(data = subset(mtfnew2,source=="cc"& new.cohort=="2"), fill="#018571",color="black",position="identity",   show.legend=TRUE)+
  geom_histogram(data = subset(mtfnew2,source=="tows"), fill="#5e3c99",color="#5e3c99",position="identity",   alpha = 0.2, show.legend=TRUE)+
  geom_vline(aes(xintercept=split), color="black", linetype="dashed", size=0.5) +
  facet_wrap(~day,scales="free_y",nrow=6)+
  theme_classic()+
  labs(x = "length (mm)",
       title = "YOY cohorts by length",
       subtitle = "Mattituck 2015") 


##Shinnecock!
shimelt <-filter(allmelt, Bay=="Shinnecock") %>%dplyr::rename(catch.date=Date, TL=value) %>%mutate(source="tows")
shico <-read.csv(file="shinextracted.csv", header=TRUE)
shico <-dplyr::select(shico,-Age,-old.pop)%>%mutate(source="cc")

weeks <-dplyr::select(shimelt,catch.date,Week)
shico <-left_join(shico,weeks)
shifnew <-full_join(shico, shimelt, by=c("catch.date","TL","Bay","Year","source","Week")) %>% mutate(Year=as.factor(Year),catch.date=as.Date(catch.date))

split <-c(36,	40,	53,	55,	63,	65,	80,	82)
splits16 <- as.data.frame(split) %>% mutate(Year="2016")
splits16$catch.date <- as.Date(c("2016-06-17",	"2016-06-24",	"2016-07-01",	"2016-07-08",	"2016-07-13",	"2016-07-21",	"2016-07-28",	"2016-08-04"))
split <-c(33,	35,	40,	47,	51,	68,	72,	85)
splits17 <- as.data.frame(split) %>% mutate(Year="2017")
splits17$catch.date <- as.Date(c("2017-06-02",	"2017-06-09",	"2017-06-14",	"2017-06-22",	"2017-06-28",	"2017-07-05",	"2017-07-13",	"2017-07-19"))
splits <-bind_rows(splits16,splits17) %>% mutate(Year=as.factor(Year),catch.date=as.Date(catch.date))

shifnew <-full_join(shifnew, splits, by=c("Year","catch.date"))
shifnew <-unique(shifnew) %>% mutate(day = format(as.Date(catch.date), format="%m/%d")) %>%filter(!is.na(Week))

#both years in same plot
ggplot(shifnew,aes(x=TL)) + 
  geom_histogram(data = subset(shifnew,source=="cc"& new.cohort=="1"), fill="#e66101",color="black", position="identity",   show.legend=TRUE)+ 
  geom_histogram(data = subset(shifnew,source=="cc"& new.cohort=="2"), fill="#018571",color="black",position="identity",   show.legend=TRUE)+
  geom_histogram(data = subset(shifnew,source=="tows"), fill="#5e3c99",color="#5e3c99",position="identity",   alpha = 0.2, show.legend=TRUE)+
  geom_vline(aes(xintercept=split), color="black", linetype="dashed", size=0.5) +
  facet_grid(day~Year,scales="free_y")+
  theme_classic()+
  labs(x = "length (mm)",
       title = "YOY cohorts by length",
       subtitle = "Shinnecock")


#separate plots for each year - Remember to change out year!
shifnew2 <-filter(shifnew, Year=="2017")
ggplot(shifnew2,aes(x=TL)) + 
  geom_histogram(data = subset(shifnew2,source=="cc"& new.cohort=="1"), fill="#e66101",color="black", position="identity",   show.legend=TRUE)+ 
  geom_histogram(data = subset(shifnew2,source=="cc"& new.cohort=="2"), fill="#018571",color="black",position="identity",   show.legend=TRUE)+
  geom_histogram(data = subset(shifnew2,source=="tows"), fill="#5e3c99",color="#5e3c99",position="identity",   alpha = 0.2, show.legend=TRUE)+
  geom_vline(aes(xintercept=split), color="black", linetype="dashed", size=0.5) +
  facet_wrap(~day,scales="free_y",nrow=6)+
  theme_classic()+
  labs(x = "length (mm)",
       title = "YOY cohorts by length",
       subtitle = "Shinnecock 2017") 

## Stuff to consider ###


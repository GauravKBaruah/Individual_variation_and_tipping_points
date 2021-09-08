rm(list=ls())
source('~/02_functions_tipping_point.R')
require(deSolve) ## for integrating ordinary differential equations
require(tidyverse) ## for efficient data manipulation & plotting
require(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(Hmisc)
library(popbio)
library(effects)
library(logihist)
library(regclass)
mydir = 'Network_data'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:101]
#newfiles<-newfiles[-17]

adj.mat<-function(data){
  #dat <- paste('network.csv',sep='')
  d <- read.csv(file=data,header=FALSE )
  dat<-as.matrix(d)
  dat[dat > 0] = 1
  dat<-apply(dat,2,as.numeric)
  return(dat)}

# code to evaluate the max number of species
data_species<- expand.grid( `web` =newfiles) %>%
  as_tibble %>%
  mutate(`No.species`=0)

for( r in 1:101) {
  g<-adj.mat(myfiles[which(myfiles == data_species$web[r])]) #network web names
  #g<-g[-1,-1] 
  data_species$No.species[r]<-nrow(g)+ncol(g)  
  
}

max(data_species$No.species)
min(data_species$No.species)


# all the data types

#load("01_species_level_tipping_point_trade_off_competition_exp_distribution.RData")
load("01_species_level_tipping_point_trade_off_no_trade_off_no_noise_1.RData") #original data
#load("01_species_level_tipping_point_trade_off_strong_competition.RData")
#load("01_species_level_tipping_point_trade_off_additive_noise.RData")
#load("01_species_level_tipping_point_trade_off_multiplicative_noise_1.RData")

#new_ddf<-new_ddf %>% filter(Web != "plant_pollinator/M_PL_046.csv")

for(i in 1:nrow(new_ddf)){
  
  g<-adj.mat(myfiles[which(myfiles == new_ddf$Web[i])]) #network web names
  #g<-g[-1,-1] 
  new_ddf$network.size[i]<- nrow(g)+ncol(g)
  
}


new_ddf$Web<-as.factor(new_ddf$Web)
new_ddf$Species<-as.numeric(as.character(new_ddf$Species))
new_ddf$Tipping_points<-as.numeric(as.character(new_ddf$Tipping_points))
new_ddf$Nestedness<-as.numeric(as.character(new_ddf$Nestedness))
new_ddf$Connectance<-as.numeric(as.character(new_ddf$Connectance))
new_ddf$Abrupt_collapse<-as.numeric(as.character(new_ddf$Abrupt_collapse))
str(new_ddf)

myfiles<-myfiles[-102]
#myfiles<-myfiles[-17]
no.of.species_trade_off_h<-no.of.species_trade_off_l<-no.of.species_ntrade_off_h<-no.of.species_ntrade_off_l<-numeric()
for(i in 1:101){
  
  temp_trade_off_h<-new_ddf %>% filter(Web ==myfiles[i], Individual_variation == "high",
                           Interaction_type == "trade_off")
  
  temp_trade_off_l<-new_ddf %>% filter(Web ==myfiles[i], Individual_variation == "low",
                           Interaction_type == "trade_off")
  
  temp_n_trade_off_h<-new_ddf %>% filter(Web ==myfiles[i], Individual_variation == "high",
                           Interaction_type == "trade_off")
  
  temp_n_trade_off_l<-new_ddf %>% filter(Web ==myfiles[i], Individual_variation == "low",
                           Interaction_type == "no_trade_off")
  
  which(temp_trade_off_h$Tipping_points>= 0 )
  no.of.species_trade_off_h[i]<- length(which(temp_trade_off_h$Tipping_points >= 0))/length(temp_trade_off_h$Species)
  no.of.species_trade_off_l[i]<- length(which(temp_trade_off_l$Tipping_points >=0))/length(temp_trade_off_l$Species)
  no.of.species_ntrade_off_h[i]<- length(which(temp_n_trade_off_h$Tipping_points>=0))/length(temp_n_trade_off_h$Species)
  no.of.species_ntrade_off_l[i]<- length(which(temp_n_trade_off_l$Tipping_points>=0))/length(temp_n_trade_off_l$Species)
  
  
  
  
}

#data frame
species_tipping_point_dat<-data.frame(fraction_of_species = c(no.of.species_trade_off_h,no.of.species_trade_off_l,
                                   no.of.species_ntrade_off_h,no.of.species_ntrade_off_l),
  Individual_variation =c(rep("High",length(no.of.species_trade_off_h)),rep("Low",length(no.of.species_trade_off_h)), 
                          rep("High",length(no.of.species_trade_off_h)), rep("Low",length(no.of.species_trade_off_h))),
  Interaction_type =c(rep("trade_off",(length(no.of.species_trade_off_h)+length(no.of.species_trade_off_h))),
                      rep("no_trade_off", (length(no.of.species_trade_off_h)+length(no.of.species_trade_off_h)))))



species_tipping_point_dat_1<- species_tipping_point_dat %>% filter(Interaction_type == "trade_off")
ddd<-species_tipping_point_dat_1 %>% group_by(Individual_variation) %>% summarise(mean_=mean(fraction_of_species,na.rm=T),
                                                                             n=n(),
                                                                             SD_=sd(fraction_of_species,na.rm=T)/sqrt(n))


ddd
#17% +- 2.3% <- trade-off for high variation
#17.1% 1- 2.31%<-no trade-off for high variation
#0% for trade-off low variation
#2.96% for 0.7% for no-trade-off low variation

m2_species<-ggplot(species_tipping_point_dat_1, aes(x=Individual_variation, y=(fraction_of_species), 
                                            color=Individual_variation))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2),size=3.5, alpha=0.5)+
  scale_fill_manual(values=c( "#E69F00", "#56B4E9"))+ggtitle("C")+
  theme_cowplot()+ 
  ylab("Fraction of species in a network \n exhibiting abrupt collapse") +
  xlab("")+ 
  scale_color_manual(values = c( "#E69F00", "#56B4E9"))+
  ylim(c(0,1))+
  xlab("")+ 
  theme(legend.position = "none")#+labs(color='individual variation') 

m2_species



#### ANALYSIS of TIPPING POINTS AT THE NETWORK LEVEL #####

#### ANALYSIS of TIPPING POINTS AT THE NETWORK LEVEL #####


load("01_network_level_tipping_point_trade_off_no_trade_off_no_noise_1.RData") #original data
#load("01_Appendix_network_level_tipping_point_trade_off_competition_exp_distribution.RData")
#load("01_Appendix_network_level_tipping_point_trade_off_strong_competition.RData")
#load("01_Appendix_network_level_tipping_point_trade_off_additive_noise.RData")
#load("01_Appendix_network_level_tipping_point_trade_off_multiplicative_noise_1.RData")


new_gamma_df$Web<-as.factor(new_gamma_df$Web)
new_gamma_df$Total_network_abundance<-as.numeric(as.character(new_gamma_df$Total_network_abundance))
new_gamma_df$Mutualistic_strength<-as.numeric(as.character(new_gamma_df$Mutualistic_strength))
new_gamma_df$Nestedness<-as.numeric(as.character(new_gamma_df$Nestedness))
new_gamma_df$Connectance<-as.numeric(as.character(new_gamma_df$Connectance))
new_gamma_df$Trait_matching<-as.numeric(as.character(new_gamma_df$Trait_matching))


#new_gamma_df<-new_gamma_df %>% filter(Web != "plant_pollinator/M_PL_046.csv")
#myfiles<-myfiles[-102]
#myfiles<-myfiles[-17]
for(i in 1:nrow(new_gamma_df)){
  
  g<-adj.mat(myfiles[which(myfiles == new_gamma_df$Web[i])]) #network web names
  #g<-g[-1,-1] 
  new_gamma_df$network.size[i]<- nrow(g)+ncol(g)
  
}





# mutualistic webs going through a tipping point shown for trade-off only criteria
new_gamma_df_trade_off<-new_gamma_df %>% filter(Interaction_type == "trade_off" )


m1<-ggplot(new_gamma_df_trade_off, aes(y = (Total_network_abundance),
                             x = Mutualistic_strength,
                             colour = Individual_variation))+
  geom_point(size = 1.25)+ggtitle("A")+
  scale_color_manual(values=c( "#E69F00", "#56B4E9"))+
  theme_cowplot()+theme(legend.position="bottom")+
  labs(x = expression(paste("Mutualistic strength, ",gamma[0])))+
  ylab("Equilibrium network \n abundance")+
  geom_smooth(method="loess",se = T)#+facet_wrap(.~Interaction_type)

m1


#summarizing mean trait matching 
fact_m_t_matching<-new_gamma_df_trade_off %>% group_by(Individual_variation, Mutualistic_strength) %>% 
  dplyr::summarize(mean_matching = mean(Trait_matching, na.rm = TRUE),
                   sd_matching=1.96*sd(Trait_matching, na.rm = TRUE)/sqrt(101))



t1<-ggplot(fact_m_t_matching, aes(y = (mean_matching),
                                  x = Mutualistic_strength,
                                  colour = Individual_variation))+
  geom_pointrange(aes(ymin=mean_matching-sd_matching, ymax=mean_matching+sd_matching))+
  scale_color_manual(values=c( "#E69F00", "#56B4E9"))+ggtitle("A")+
  theme_cowplot()+ ylab("Mean trait matching")+ labs(x = expression(paste("Mutualistic strength, ",gamma[0])))+
  labs(color='individual variation') #+facet_wrap(.~Interaction_type)
  

t1


## function for evaluating the point of collapses, frequency of collapse
collapse_point<-function(webs,data){
  # filtering web names and then finding the point at which the web of network collapses
  collapse.point.highvar.trade_off<-collapse.point.lowvar.trade_off<-collapse.point.highvar.ntrade_off<-collapse.point.lowvar.ntrade_off<-
    connc.h<-connc.l<-nestd.h.trade_off<-nestd.l.trade_off<-nestd.h.ntrade_off<-nestd.l.ntrade_off<-connc.h.trade_off<-connc.l.trade_off<-connc.h.ntrade_off<-connc.l.ntrade_off<-
    collapse.net.high.trade_off<-collapse.net.low.trade_off<-collapse.net.high.ntrade_off<-collapse.net.low.ntrade_off<-numeric()
  network.size.h.ntrade_off<-network.size.l.ntrade_off<-network.size.h.trade_off<-network.size.l.trade_off<-numeric()
 
    
    if(all(data$Interaction_type == "trade_off") == TRUE){
      
      
      for(i in 1:length(webs)){
      t1.high.trade_off<- filter(data, Web == webs[i] & Individual_variation == "high" & Interaction_type == "trade_off")
      t2.low.trade_off<-filter(data, Web == webs[i] & Individual_variation == "low" & Interaction_type == "trade_off")
      
      
      collapse.point.highvar.trade_off[i]<-max(t1.high.trade_off$Mutualistic_strength[which(round( (t1.high.trade_off$Total_network_abundance),2) < 1)])
      collapse.point.lowvar.trade_off[i]<-max(t2.low.trade_off$Mutualistic_strength[which(round((t2.low.trade_off$Total_network_abundance),2) < 1)])
      
      
      nestd.h.trade_off[i]<-t1.high.trade_off$Nestedness[1]
      nestd.l.trade_off[i]<-t2.low.trade_off$Nestedness[1]
      
      
      connc.h.trade_off[i]<-t1.high.trade_off$Connectance[1]
      connc.l.trade_off[i]<-t2.low.trade_off$Connectance[1]
      
      
      network.size.h.trade_off[i]<-t1.high.trade_off$network.size[1]
      network.size.l.trade_off[i]<-t2.low.trade_off$network.size[1]
      
      collapse.net.high.trade_off[i]<-rate.of.collapse(t1.high.trade_off)
      collapse.net.low.trade_off[i]<-rate.of.collapse(t2.low.trade_off)
      
      ## high, low and medium variation
    }
    freq.of.abrupt.high.tradeoff <-sum(collapse.net.high.trade_off)/length(collapse.net.high.trade_off)
    freq.of.abrupt.low.tradeoff <-sum(collapse.net.low.trade_off)/length(collapse.net.low.trade_off)
    
    output<- list(freq.of.abrupt.high.tradeoff=freq.of.abrupt.high.tradeoff,
                  freq.of.abrupt.low.tradeoff=freq.of.abrupt.low.tradeoff,
                  
                  network.size.h.trade_off=network.size.h.trade_off,
                  network.size.l.trade_off=network.size.l.trade_off,
                  
                  collapse.net.high.trade_off=collapse.net.high.trade_off,
                  collapse.net.low.trade_off=collapse.net.low.trade_off,
                  
                  collapse.point.highvar.trade_off=collapse.point.highvar.trade_off,
                  collapse.point.lowvar.trade_off=collapse.point.lowvar.trade_off,
                  
                  nestd.h.trade_off=nestd.h.trade_off,
                  nestd.l.trade_off=nestd.l.trade_off,
                  
                  connc.h.trade_off=connc.h.trade_off,
                  connc.l.trade_off=connc.l.trade_off)
                 
      
      
    }else {
      for(i in 1:length(webs)){
      t1.high.trade_off<- filter(data, Web == webs[i] & Individual_variation == "high" & Interaction_type == "trade_off")
      t2.low.trade_off<-filter(data, Web == webs[i] & Individual_variation == "low" & Interaction_type == "trade_off")
      
      t1.high.ntrade_off<- filter(data, Web == webs[i] & Individual_variation == "high" & Interaction_type == "no_trade_off")
      t2.low.ntrade_off<-filter(data, Web == webs[i] & Individual_variation == "low" & Interaction_type == "no_trade_off")
      
      collapse.point.highvar.trade_off[i]<-max(t1.high.trade_off$Mutualistic_strength[which(round( (t1.high.trade_off$Total_network_abundance),1) < 1)])
      collapse.point.lowvar.trade_off[i]<-max(t2.low.trade_off$Mutualistic_strength[which(round((t2.low.trade_off$Total_network_abundance),1) < 1)])
      
      collapse.point.highvar.ntrade_off[i]<-max(t1.high.ntrade_off$Mutualistic_strength[which(round( (t1.high.ntrade_off$Total_network_abundance),1) < 1)])
      collapse.point.lowvar.ntrade_off[i]<-max(t2.low.ntrade_off$Mutualistic_strength[which(round((t2.low.ntrade_off$Total_network_abundance),1) < 1)])
      
      nestd.h.trade_off[i]<-t1.high.trade_off$Nestedness[1]
      nestd.l.trade_off[i]<-t2.low.trade_off$Nestedness[1]
      
      nestd.h.ntrade_off[i]<-t1.high.ntrade_off$Nestedness[1]
      nestd.l.ntrade_off[i]<-t2.low.ntrade_off$Nestedness[1]
      
      connc.h.trade_off[i]<-t1.high.trade_off$Connectance[1]
      connc.l.trade_off[i]<-t2.low.trade_off$Connectance[1]
      
      connc.h.ntrade_off[i]<-t1.high.ntrade_off$Connectance[1]
      connc.l.ntrade_off[i]<-t2.low.ntrade_off$Connectance[1]
      
      network.size.h.trade_off[i]<-t1.high.trade_off$network.size[1]
      network.size.l.trade_off[i]<-t2.low.trade_off$network.size[1]
      
      network.size.h.ntrade_off[i]<-t1.high.ntrade_off$network.size[1]
      network.size.l.ntrade_off[i]<-t2.low.ntrade_off$network.size[1]
      
      
      collapse.net.high.trade_off[i]<-rate.of.collapse(t1.high.trade_off)
      collapse.net.low.trade_off[i]<-rate.of.collapse(t2.low.trade_off)
      
      
      collapse.net.high.ntrade_off[i]<-rate.of.collapse(t1.high.ntrade_off)
      collapse.net.low.ntrade_off[i]<-rate.of.collapse(t2.low.ntrade_off)
      
      ## high, low and medium variation
    }
    freq.of.abrupt.high.tradeoff <-sum(collapse.net.high.trade_off)/length(collapse.net.high.trade_off)
    freq.of.abrupt.low.tradeoff <-sum(collapse.net.low.trade_off)/length(collapse.net.low.trade_off)
    freq.of.abrupt.high.no.tradeoff<-sum(collapse.net.high.ntrade_off)/length(collapse.net.high.ntrade_off)
    freq.of.abrupt.low.no.tradeoff<-sum(collapse.net.low.ntrade_off)/length(collapse.net.low.ntrade_off)
    
    output<- list(freq.of.abrupt.high.tradeoff=freq.of.abrupt.high.tradeoff,
                freq.of.abrupt.low.tradeoff=freq.of.abrupt.low.tradeoff,
                freq.of.abrupt.high.no.tradeoff=freq.of.abrupt.high.no.tradeoff,
                freq.of.abrupt.low.no.tradeoff=freq.of.abrupt.low.no.tradeoff,
                
                network.size.h.trade_off=network.size.h.trade_off,
                network.size.l.trade_off=network.size.l.trade_off,
                network.size.h.ntrade_off=network.size.h.ntrade_off,
                network.size.l.ntrade_off=network.size.l.ntrade_off,
                
                collapse.net.high.trade_off=collapse.net.high.trade_off,
                collapse.net.low.trade_off=collapse.net.low.trade_off,
                collapse.net.high.ntrade_off=collapse.net.high.ntrade_off,
                collapse.net.low.ntrade_off=collapse.net.low.ntrade_off,
                
                collapse.point.highvar.trade_off=collapse.point.highvar.trade_off,
                collapse.point.lowvar.trade_off=collapse.point.lowvar.trade_off,
                collapse.point.highvar.ntrade_off=collapse.point.highvar.ntrade_off,
                collapse.point.lowvar.ntrade_off=collapse.point.lowvar.ntrade_off,
                
                nestd.h.trade_off=nestd.h.trade_off,
                nestd.l.trade_off=nestd.l.trade_off,
                nestd.h.ntrade_off=nestd.h.ntrade_off,
                nestd.l.ntrade_off=nestd.l.ntrade_off,
                
                connc.h.trade_off=connc.h.trade_off,
                connc.l.trade_off=connc.l.trade_off,
                connc.h.ntrade_off=connc.h.ntrade_off,
                connc.l.ntrade_off=connc.l.ntrade_off)
    
    }
    
      return(output)
}



#function for evaluating whether a collapse was abrupt or not
rate.of.collapse<-function(dat)
{
  
  dbiomass <- (dat$Total_network_abundance[2:length(dat$Total_network_abundance)]- 
                    dat$Total_network_abundance[1:(length(dat$Total_network_abundance)-1)])
  
  dbiomass[is.na(dbiomass)]<-0
  
  if(max(dbiomass) > 45){
    collapse=1
  }else {collapse = 0
  }
  return(collapse)
  
}


collapse.points<-collapse_point(webs=newfiles, data=new_gamma_df)



#creating the data frame for summarizing abrupt collapse for individual variation and different tradeoffs
abrupt_collapse_dat<-data.frame(Abrup_collapse=c(collapse.points$freq.of.abrupt.high.tradeoff,
                                                 collapse.points$freq.of.abrupt.low.tradeoff,
                                                 collapse.points$freq.of.abrupt.high.no.tradeoff,
                                                 collapse.points$freq.of.abrupt.low.no.tradeoff
),
indvariation=c(rep("high",1),rep("low",1)),
int_type=c("trade-off", "trade-off", "no trade-off", "no trade-off"))


abrupt_collapse_dat_trade_off<-abrupt_collapse_dat %>% filter(int_type =="trade-off")

m2<-ggplot(abrupt_collapse_dat_trade_off, aes(x=indvariation, y=Abrup_collapse, fill=indvariation))+
  geom_bar(position = "dodge",stat="identity", color="black")+
  scale_fill_manual(values=c( "#E69F00", "#56B4E9"))+ggtitle("B")+
  theme_cowplot()+ ylab("Abrupt collapse in \n proportion of networks") +
  xlab("")+ ylim(c(0,1))+
  theme(legend.position = "none")#+labs(fill='individual variation') 

m2


#logistic regression : creating a data frame of abrupt collapse vs nestedness, connectance, individual variation
abrupt_collapse_for_glm<-data.frame(Abrupt_collapse=c(collapse.points$collapse.net.high.trade_off,
                                                      collapse.points$collapse.net.low.trade_off,
                                                      collapse.points$collapse.net.high.ntrade_off,
                                                      collapse.points$collapse.net.low.ntrade_off
),
indvariation=c(rep("high",length(collapse.points$collapse.net.high.trade_off)),
               rep("low",length(collapse.points$collapse.net.low.trade_off)),
               rep("high",length(collapse.points$collapse.net.high.ntrade_off)),
               rep("low",length(collapse.points$collapse.net.low.ntrade_off))),
type=c(rep("trade-off",202), rep("no trade-off",202)),
network.size=c(collapse.points$network.size.h.trade_off,
               collapse.points$network.size.l.trade_off,
               collapse.points$network.size.h.ntrade_off,
               collapse.points$network.size.l.ntrade_off),
nestedness=c(collapse.points$nestd.h.trade_off,
             collapse.points$nestd.l.trade_off,
             collapse.points$nestd.h.ntrade_off,
             collapse.points$nestd.l.ntrade_off #collapse.points$nestd.l.evoyes,collapse.points$nestd.l.evono
),
connectance=c(collapse.points$connc.h.trade_off,
              collapse.points$connc.l.trade_off,
              collapse.points$connc.h.ntrade_off,
              collapse.points$connc.l.ntrade_off
              #collapse.points$connc.l.evoyes,collapse.points$connc.l.evono
))


abrupt_collapse_tradeoff<-abrupt_collapse_for_glm %>% filter(type == "trade-off")


#glm
glm.fit<-glm(Abrupt_collapse~ nestedness+indvariation+network.size ,
             family = 'binomial', data = abrupt_collapse_tradeoff)

summary(glm.fit)

glm.fit$deviance / glm.fit$df.residual

#trying to plot the same through ggplot2
m3_nestedness<-logihist(abrupt_collapse_tradeoff$nestedness, abrupt_collapse_tradeoff$Abrupt_collapse,
             scale.hist = 3, breaks = "Sturges", 
             counts = TRUE, intervalo = 0,
             ylab2 = "Frequency", fillb = 4, colob = 4, sizeb = 0.5, pglm = FALSE, se = FALSE,
             sizeglm = 1, colglm = 1)+ggtitle("B")+
  stat_smooth(method = "glm", formula = y ~x, method.args = list(family = "quasibinomial"))+ 
  theme_cowplot()+xlab("Nestedness (NODF)") + ylab(" Chances of abrupt collapse")
m3_nestedness

m5_networksize<-logihist(abrupt_collapse_tradeoff$network.size, abrupt_collapse_tradeoff$Abrupt_collapse,
             scale.hist = 3, breaks = "Sturges", 
             counts = TRUE, intervalo = 0,
             ylab2 = "Frequency", fillb = 4, colob = 4, sizeb = 0.5, pglm = FALSE, se = FALSE,
             sizeglm = 1, colglm = 1)+ggtitle("C")+
  stat_smooth(method = "glm", formula = y ~x, method.args = list(family = "quasibinomial"))+ 
  theme_cowplot()+xlab("Network size") + ylab(" Chances of abrupt collapse")

m5_networksize


m7<-logihist(abrupt_collapse_tradeoff$connectance, abrupt_collapse_tradeoff$Abrupt_collapse,
             scale.hist = 5, breaks = "Sturges", 
             counts = TRUE, intervalo = 0,
             ylab2 = "Frequency", fillb = 4, colob = 4, sizeb = 0.5, pglm = FALSE, se = FALSE,
             sizeglm = 1, colglm = 1)+ggtitle("B")+
  stat_smooth(method = "glm", formula = y ~x, method.args = list(family = "binomial"))+ 
  theme_cowplot()+xlab("Connectance") + ylab(" Chances of abrupt collapse")

m7



plot(allEffects(glm.fit))
plot(predictorEffects(glm.fit, ~ indvariation*network.size),
     axes=list(grid=F,
               x=list(rug=FALSE),
               y=list(type="response"),
               main = "", xlab ="Nestedness", ylab="Probability of abrupt collapse"))



#data frame for the time point of collapse
tempo<-data.frame(point_of_collapse=c(collapse.points$collapse.point.highvar.trade_off,
                                      collapse.points$collapse.point.lowvar.trade_off,
                                      collapse.points$collapse.point.highvar.ntrade_off,
                                      collapse.points$collapse.point.lowvar.ntrade_off),
                  indvariation=c(rep("high",length(collapse.points$collapse.point.highvar.trade_off)),
                                 rep("low",length(collapse.points$collapse.point.lowvar.trade_off)),
                                 rep("high",length(collapse.points$collapse.point.highvar.ntrade_off)),
                                 rep("low",length(collapse.points$collapse.point.lowvar.ntrade_off))),
                  type=c(rep("trade-off",202), rep("no trade-off",202)),
                  nestedness=c(collapse.points$nestd.h.trade_off,
                               collapse.points$nestd.l.trade_off,
                               collapse.points$nestd.h.ntrade_off,
                               collapse.points$nestd.l.ntrade_off #collapse.points$nestd.l.evoyes,collapse.points$nestd.l.evono
                  ),
                  connectance=c(collapse.points$connc.h.trade_off,
                                collapse.points$connc.l.trade_off,
                                collapse.points$connc.h.ntrade_off,
                                collapse.points$connc.l.ntrade_off
                  ))

temp_trade_off<-tempo %>% filter(type == "trade-off")

m4<-ggplot(data=temp_trade_off, 
           aes(x=indvariation,y=point_of_collapse))+
  geom_violin(aes(fill = indvariation), trim = FALSE) + 
  geom_boxplot(width = 0.2)+
  ggtitle("D")+
  geom_jitter(shape=1, size=3.5,alpha=0.75,position=position_jitter(0.2))+
  scale_fill_manual(values = c( "#E69F00", "#56B4E9"))+xlab("")+
  theme_cowplot() + coord_flip()+
  labs(y = expression(paste("Point of collapse, ",gamma[0])))+ theme(legend.position = "none")
m4


#summarising mean and std error for point of collapse
tempo %>% group_by(indvariation,type) %>%
  dplyr::summarize(mean_size = mean(point_of_collapse, na.rm = TRUE),
                   sd_size=1.96*sd(point_of_collapse, na.rm = TRUE)/sqrt(101))



#plotting everything together
lay_out(
      list(m1, 1, 1),
       list(m2, 1, 2),
       list(m2_species, 1, 3),
      list(m4, 1, 4),
      list(t1,2,1:2),
     list(m3_nestedness,2,1),
      list(m5_networksize,2,2)
      )



#figure 2
#grid.arrange(m2,m4, nrow=2)

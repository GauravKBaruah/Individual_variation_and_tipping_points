#code for Ecology Letters paper: 
#" The impact of individual variation on abrupt collapses in mutualistic networks" 2021. Gaurav Baruah
#email: gbaruahecoevo@gmail.com



rm(list=ls())
source('~/02_functions_tipping_point.R', echo=F)#converting mutualism matrix to ones and zeros
require(deSolve) ## for integrating ordinary differential equations
require(tidyverse) ## for efficient data manipulation & plotting
require(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
require(ggplot2)
require(reshape2)
require(flux)
require(akima)
library(directlabels)
library(gridExtra)
library(grid)
library(igraph)
library(bipartite)


#five mutualistic networks :60, 34, 8
#nestedness: 0.635, 0.246, 0

newfiles <-c("plant_pollinator/M_PL_046.csv",
             "plant_pollinator/M_PL_061_33.csv", 
             "plant_pollinator/M_PL_061_18.csv")

mydir = 'plant_pollinator'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
#myfiles<-myfiles[1:101]


fact<- expand.grid(`Strength_mutualism`=seq(0,7.5,0.4), 
                   `web` =newfiles,
                   `interaction_type`= "trade_off",
                   `noise`="none",
                   `individual.variation` = c("high","low"),
                   `range_competition` = c(50,70,100,200,500,
                                           700,1000,2000,5000,7000,10000),
                   `random_seed`=4327+(1:30)*100) %>%
  as_tibble %>%
  mutate(feasibility = 0,
         richness = 0)




model.t<-list()



for(r in 1:nrow(fact)){
  #print(r)
  
  # extracting the interaction matrix
  g<-adj.mat(myfiles[which(myfiles == fact$web[r])]) #network web names
  # g<-g[-1,-1] 
  
  
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] # no of plant species
  degree.animals<-degree.plants<-numeric()
  
  #degree of plants and anichmals
  for(i in 1:Plantspecies){
    degree.plants[i]<-sum(g[i,])} # degree of plants
  for(j in 1:Aspecies){
    degree.animals[j]<-sum(g[,j]) # degree of animals
  }
  
  ##control loop for selecting whether variation is high or low
  if(fact$individual.variation[r] == "low"){
    sig <-runif((Aspecies+Plantspecies),0.0001,0.001)}else if(fact$individual.variation[r] == "high"){
      sig <-runif((Aspecies+Plantspecies),0.05,0.5)}
  
  
  #heritability used in the model
  h2<-runif((Aspecies+Plantspecies),0.4,0.4)
  
  ## vector of species trait standard deviations
  N <- runif( (Aspecies+Plantspecies) , 1,1)  ## initial species densities
  nainit<- N[1:Aspecies]
  npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
  
  #vector of mean phenotypic values sampled from a uniform distribution
  muinit <-runif((Aspecies+Plantspecies), -1,1)
  mainit<-muinit[1:Aspecies]
  mpinit<-muinit[(Aspecies+1): (Aspecies+Plantspecies)]
  
  #matrix of coefficients for competition
  Amatrix<-mat.comp_feasibility(g,strength=fact$range_competition[r])$Amatrix
  Pmatrix<-mat.comp_feasibility(g,strength=fact$range_competition[r])$Pmatrix
  gamma=0.5
  mut.strength<-fact$Strength_mutualism[r]
  nestedness<-nestedness_NODF(g)
  C<-Connectance(g)
  web.name<-myfiles[r]
  
  #growth rates of species
  ba<-runif(Aspecies, -0.05,-0.05)
  bp<-runif(Plantspecies,-0.05,-0.05)
  dganimals<-degree.animals
  dgplants<-degree.plants
  
  
  ic <-c(nainit, npinit, mainit,mpinit)
  
  
  params <- list(time=time,matrix=g,sig=sig,Amatrix=Amatrix,
                 Pmatrix=Pmatrix,w=gamma,
                 mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                 web.name=web.name,h2=h2, ba=ba,bp=bp,dganimals=dganimals,
                 dgplants=dgplants,interaction_type=fact$interaction_type[r],
                 noise=fact$noise[r])
  
  
  
  start.time = 3000
  model.t<-lapply(1, Mcommunity,time=start.time,state=ic,
                  pars=params)
  
  abvector<-c(model.t[[1]]$Plants[3000,],model.t[[1]]$Animals[3000,])
  percent_survived<-length(which(abvector > 0.001))/length(abvector)
  
  
  fact$feasibility[r] <-  percent_survived
  fact$richness[r] <- length(which(abvector>1e-3))
  
  print(r)
  
  
}
#save(fact, file="Feasibility_data.RData")


# plot and aanalyis for the feasibility data

#load("Feasibility_data.RData")


for(i in 1:nrow(fact)){
  
  g<-adj.mat(myfiles[which(myfiles == fact$web[i])]) #network web names
  #g<-g[-1,-1] 
  fact$network.size[i]<- nrow(g)+ncol(g)
  fact$nesdtedness[i]<- nestedness_NODF(g)
  fact$connectance[i]<-Connectance(g)
}


#web1
a1<-fact %>% filter(individual.variation == "low", web == "plant_pollinator/M_PL_046.csv") %>%  
  group_by(range_competition,Strength_mutualism,interaction_type) %>% 
  summarise(mean_feasibility = mean(feasibility, na.rm = TRUE))


plot1<-feasibility_plot(a1)

plot1


a2<-fact %>% filter(individual.variation == "high", web == "plant_pollinator/M_PL_046.csv") %>%  
  group_by(range_competition,Strength_mutualism) %>% 
  summarise(mean_feasibility = mean(feasibility, na.rm = TRUE))


plot2<-feasibility_plot(a2)

plot2



# web2


a5<-fact %>% filter(individual.variation == "low") %>%  
  group_by(range_competition,Strength_mutualism) %>% 
  summarise(mean_feasibility = mean(feasibility, na.rm = TRUE))


plot5<-feasibility_plot(a5)

plot5


a6<-fact %>% filter(individual.variation == "high") %>%  
  group_by(range_competition,Strength_mutualism) %>% 
  summarise(mean_feasibility = mean(feasibility, na.rm = TRUE))


plot6<-feasibility_plot(a6)

plot6


#web3

a3<-fact %>% filter(individual.variation == "low", web == "plant_pollinator/M_PL_061_18.csv") %>%  
  group_by(range_competition,Strength_mutualism) %>% 
  summarise(mean_feasibility = mean(feasibility, na.rm = TRUE))


plot3<-feasibility_plot(a3)

plot3


a4<-fact %>% filter(individual.variation == "high", web == "plant_pollinator/M_PL_061_18.csv") %>%  
  group_by(range_competition,Strength_mutualism) %>% 
  summarise(mean_feasibility = mean(feasibility, na.rm = TRUE))


plot4<-feasibility_plot(a4)

plot4



net1<-adj.mat("plant_pollinator/M_PL_046.csv")
net3<-adj.mat("plant_pollinator/M_PL_061_33.csv")
net2<-adj.mat("plant_pollinator/M_PL_061_18.csv")

par(mfrow=(c(3,1)))
web1<-plotweb(net1,
        method="normal",ybig=0.1, y.width.low = 0.1,
        col.interaction="wheat4",
        bor.col.interaction="white", 
        arrow="no",  col.high="lightblue",
        col.low="tomato",labsize=0.1)
web2<-plotweb(net2,
              method="normal",ybig=0.1, y.width.low = 0.1,
              col.interaction="wheat4",
              bor.col.interaction="white", 
              arrow="no",  col.high="lightblue",
              col.low="tomato",labsize=0.1)
web3<-plotweb(net3,
              method="normal",ybig=0.1, y.width.low = 0.1,
              col.interaction="wheat4",
              bor.col.interaction="white",
              arrow="no",  col.high="lightblue",
              col.low="tomato",labsize=0.1)
  
grid_arrange_shared_legend(plot1,plot2,
                           plot3,plot4,
                           plot5,plot6,
                           nrow=3,ncol=2)




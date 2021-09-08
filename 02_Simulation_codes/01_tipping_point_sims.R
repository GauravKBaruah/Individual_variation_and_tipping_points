#code for Ecology Letters paper: 
#" The impact of individual variation on abrupt collapses in mutualistic networks" 2021. Gaurav Baruah
#email: gbaruahecoevo@gmail.com


rm(list=ls())
source('~/02_functions_tipping_point.R', echo=F)#converting mutualism matrix to ones and zeros
require(deSolve) ## for integrating ordinary differential equations
require(tidyverse) ## for efficient data manipulation & plotting
require(cowplot) ## for arranging plots in a grid
require(dplyr)
require(readr)
require(beepr)
require(viridis)


theme_set(theme_classic()) 

# reading all the datasets
# calculating nestedness and connectance
mydir = '~/Network data/' #path to the directory where all the network .csv files are located
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)

newfiles<-myfiles[1:101]

#empty dataframe creation
fact<- expand.grid( `web` =newfiles,
                    `noise`=c("no_noise"), # could be also c("additive", "multiplicative", "no_noise)
                    `interaction_type`= "trade_off", # could be also c("trade_off", "no_trade off")
                    `individual.variation` = c("high","low"),
                    `evolution`="yes",
                    `random_seed`=4327+1*100) %>%
  as_tibble %>%
  mutate(`tipping.point`=0,
         degree.plants=0,
         degree.animals=0,
         nestedness=0,
         connectance=0)
set.seed(1234)
model.t<-list()
new_ddf<-NULL
new_gamma_df<-NULL
mut.strength <- seq(0,7,0.15)
for(r in 1:nrow(fact)){
  #print(r)
  
  
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
    sig <-runif((Aspecies+Plantspecies),0.0005,0.001)}else if(fact$individual.variation[r] == "high"){
      sig <-runif((Aspecies+Plantspecies),0.05,0.5)}
  
  
  #heritabilities
  h2<-runif((Aspecies+Plantspecies),0.4,0.4)
  
  ## initial species densities
  N <- runif( (Aspecies+Plantspecies) , 1,1)  
  nainit<- N[1:Aspecies]
  npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
  
  #initial mean trait values
  muinit <-runif((Aspecies+Plantspecies), -1,1)
  mainit<-muinit[1:Aspecies]
  mpinit<-muinit[(Aspecies+1): (Aspecies+Plantspecies)]
  
  #competition coefficients  
  Amatrix<-mat.comp(g)$Amatrix
  Pmatrix<-mat.comp(g)$Pmatrix
  
  #w^2 value
  gamma=0.5
  
  #gamma_0 value
  mut.strength<- mut.strength
  
  nestedness<-nestedness_NODF(g)
  C<-Connectance(g)
  
  web.name<-fact$web[r]
  
  #growth rates
  ba<-runif(Aspecies, -0.05,-0.05)
  bp<-runif(Plantspecies,-0.05,-0.05)
  
  #degrees
  dganimals<-degree.animals
  dgplants<-degree.plants
  
  
  ic <-c(nainit, npinit, mainit,mpinit)
  
  #fact$noise
  params <- list(time=time,matrix=g,sig=sig,Amatrix=Amatrix,
                 Pmatrix=Pmatrix,w=gamma,noise=fact$noise[r],
                 mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                 web.name=web.name,h2=h2, ba=ba,bp=bp,dganimals=dganimals,
                 dgplants=dgplants,
                 interaction_type=fact$interaction_type[r])
  
  
  
   start.time =5000
   model.t<-lapply(1, Mcommunity_1,time=start.time,state=ic,
                               pars=params)
  
  

  # plot_snapshot(Na = model.t[[1]]$Animals[1000,],
  #               Np = model.t[[1]]$Plants[1000,],
  #               m = c(model.t[[1]]$Animal.trait[1000,], model.t[[1]]$Plant.trait[1000,]),
  #              sigma =sig, moment=0, limits=c(-1, 1), res=1001)
  
  
  ddf <-as.data.frame(cbind( rep(as.character(web.name),each=(Aspecies+Plantspecies)), 
                             rep(seq(1,(nrow(g)+ncol(g)),1)),
                             c(model.t[[1]]$Na_tipping.point,model.t[[1]]$Np_tipping.point),
                             c(model.t[[1]]$degree.animals,model.t[[1]]$degree.plants),
                             rep(nestedness_NODF(g), each=((Aspecies+Plantspecies)) ),
                             rep(Connectance(g), each=((Aspecies+Plantspecies)) ),
                             rep(model.t[[1]]$collapse, each=(Aspecies+Plantspecies) ),
                             rep(as.character(fact$individual.variation[r]),each=(Aspecies+Plantspecies)),
                             rep(as.character(fact$interaction_type[r]),each=(Aspecies+Plantspecies))))
  
  
  colnames(ddf)<-c("Web","Species","Tipping_points","Degree","Nestedness", "Connectance","Abrupt_collapse",
                   "Individual_variation", "Interaction_type")
  
  
  new_ddf<-rbind(new_ddf,ddf)
  ddf_gamma<-as.data.frame(cbind( rep(as.character(web.name),each=(length(mut.strength))),   
                                  rep(as.character(fact$interaction_type[r]),each =length(mut.strength)),
                                  rep(as.character(fact$individual.variation[r]),each=length(mut.strength)),
                                  mut.strength,
                                  model.t[[1]]$trait_matching,
                                  model.t[[1]]$total.community.abundance,
                                  rep(nestedness_NODF(g), each=length(mut.strength)),
                                  rep(Connectance(g), each=length(mut.strength))))
  colnames(ddf_gamma)<-c("Web","Interaction_type","Individual_variation","Mutualistic_strength","Trait_matching",
                         "Total_network_abundance", "Nestedness","Connectance")
  
  new_gamma_df<-rbind(new_gamma_df, ddf_gamma )
  
  
  
  print(r)
  

  
}
#save(new_gamma_df, file = "01_Appendix_network_level_tipping_point_trade_off_multiplicative_noise_1.RData")
#save(new_ddf, file = "01_species_level_tipping_point_trade_off_multiplicative_noise_1.RData")


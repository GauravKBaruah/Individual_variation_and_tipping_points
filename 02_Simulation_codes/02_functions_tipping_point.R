#code for Ecology Letters paper: 
#" The impact of individual variation on abrupt collapses in mutualistic networks" 2021. Gaurav Baruah
#email: gbaruahecoevo@gmail.com
require(statmod)

cutoff <- function(x) ifelse(x<1, (1*(x>0))*(x*x*x*(10+x*(-15+6*x))), 1)



#conversion to a matrix
adj.mat<-function(data){
  #dat <- paste('network.csv',sep='')
  d <- read.csv(file=data,header=FALSE )
  dat<-as.matrix(d)
  dat[dat > 0] = 1
  dat[dat < 0] = 
    dat<-apply(dat,2,as.numeric)
  return(dat)}



#function for estimating trait matching for a network
trait.matching<-function(mA,mP,adj.mat_1,gamma){
  tm<-numeric()
  for(i in 1:nrow(adj.mat_1)){
    tm[i] <- mean(adj.mat_1[i,]*exp(-(mA-mP[i])^2)/gamma)
    
  }
  return(tm=mean(tm))
}



#gaussian quadrature approximation used to numerically determine the integrals.
gausquad.animals<-function(m,sigma,w,h,np,na,mut.strength,points,mat,degree.animal,interaction_type){
  
  
  temp2<-dat2<-x2<-x3<-array(dim=c(points))
  if(mat == 0){
   return(list(G= 0, B = 0))
}
  else if(mat == 1){
  #nodes oir points in the abscissa where the integral will be evaluated numerically
  if(interaction_type == "trade_off"){
  z1<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sigma$sa)$nodes #z'
  z2<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sigma$sp)$nodes #z''
  
  #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
  w1<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$ma,sigma =sigma$sa)$weights #pi(z')
  w2<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$mp,sigma =sigma$sp)$weights #pj(z'')
  
  
  #for the pairwise model however there are only two species interacting and hence i and j
  #or in other words the integral goes over z and z'
  for (i in 1: points){
    
    
    temp2[i]<- sum(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i])
       x2[i]<- sum(h*w2*exp(-(z1[i]-z2)^2/w^2))
       
    dat2[i]<- sum((z1[i]-m$ma)*(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]))
    
     # x3[i] <-sum(h*w2*exp(-(z1[i]-z2)^2/w^2))
      
  }
  G = sum(temp2*np*(mut.strength/degree.animal)/(1+x2*(mut.strength/degree.animal)*np))
  B = sum(dat2*np*(mut.strength/degree.animal)/(1+x2*(mut.strength/degree.animal)*np)) 
  
    } else {
      
      z1<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sigma$sa)$nodes #z'
      z2<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sigma$sp)$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w1<-gauss.quad.prob(points, dist = "normal", 
                          mu=m$ma,sigma =sigma$sa)$weights #pi(z')
      w2<-gauss.quad.prob(points, dist = "normal", 
                          mu=m$mp,sigma =sigma$sp)$weights #pj(z'')
      
      
      #for the pairwise model however there are only two species interacting and hence i and j
      #or in other words the integral goes over z and z'
      for (i in 1: points){
        
        
        temp2[i]<- sum(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i])
        x2[i]<- sum(h*w2*exp(-(z1[i]-z2)^2/w^2))
        
        dat2[i]<- sum((z1[i]-m$ma)*(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]))
        
        # x3[i] <-sum(h*w2*exp(-(z1[i]-z2)^2/w^2))
        
      }
      G = sum(temp2*np*(mut.strength/1)/(1+x2*(mut.strength/1)*np))
      B = sum(dat2*np*(mut.strength/1)/(1+x2*(mut.strength/1)*np)) 
    }
    
  
    return(list(G= G, B = B))
  }
}

#gaussian quadrature approximation used to numerically determine the integrals.
gausquad.plants<-function(m,sigma,w,h,np,na,mut.strength,points,mat,degree.plant,interaction_type){
  
  temp2<-dat2<-x3<-x4<-array(dim=c(points))
  
  if(mat == 0){
    
    return(list(G= 0, 
                B = 0))
    
  }
  else if (mat==1){
    if(interaction_type == "trade_off"){
  #nodes oir points in the abscissa where the integral will be evaluated numerically
  z1<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sigma$sp)$nodes #z'
  z2<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sigma$sa)$nodes #z''
  
  #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
  w1<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$mp,sigma =sigma$sp)$weights #pi(z')
  w2<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$ma,sigma =sigma$sa)$weights #pj(z'')
  
  
  #for the pairwise model however there are only two species interacting and hence i and j
  #or in other words the integral goes over z and z'
  for (i in 1: points){
    temp2[i]<- sum(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i])
    x3[i]<-sum(h*exp(-(z1[i]-z2)^2/w^2)*w2)
    dat2[i]<- sum((z1[i]-m$mp)*(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]))
    
   # x4[i]<-sum(h*w2*exp(-(z1[i]-z2)^2/w^2))
  }
  
  G = sum(temp2*na*(mut.strength/degree.plant)/(1+x3*(mut.strength/degree.plant)*na))
  B = sum(dat2*na*(mut.strength/degree.plant)/(1+x3*(mut.strength/degree.plant)*na))
    } else {
      
      z1<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sigma$sp)$nodes #z'
      z2<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sigma$sa)$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w1<-gauss.quad.prob(points, dist = "normal", 
                          mu=m$mp,sigma =sigma$sp)$weights #pi(z')
      w2<-gauss.quad.prob(points, dist = "normal", 
                          mu=m$ma,sigma =sigma$sa)$weights #pj(z'')
      
      
      #for the pairwise model however there are only two species interacting and hence i and j
      #or in other words the integral goes over z and z'
      for (i in 1: points){
        temp2[i]<- sum(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i])
        x3[i]<-sum(h*exp(-(z1[i]-z2)^2/w^2)*w2)
        dat2[i]<- sum((z1[i]-m$mp)*(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]))
        
        # x4[i]<-sum(h*w2*exp(-(z1[i]-z2)^2/w^2))
      }
      
      G = sum(temp2*na*(mut.strength/1)/(1+x3*(mut.strength/1)*na))
      B = sum(dat2*na*(mut.strength/1)/(1+x3*(mut.strength/1)*na))
      
    }
    
    
  return(list(G= G, 
              B = B))
  }
  
}



#this function is used specifically to evaluate tipping points for each species in a network
#and evaluate whether a species in a network went through an abrupt collapse.
# Function dynamically calculates co-evolutionary dynamics over time and evaluates
# time: timepoints of simulation
# state: initial abundance vector for plants and animals which is fixed to 1
# pars: list of parameters
#output: Na_tipping point, Np_tipping_point:  returns tipping points for animals and plants for a mutualistic network
#degrees: degree of animals and plants
# abrupt collapse: whether the mutualistic network went through an abrupt collapse
# point.of.network.collapse: mutualistic strength at which the total equilibrium network abundance falls below 1.
#trait_matching: estimates average trait matching for a mutualistic network
#total.community.abundance: equlibrium community abundance

eqs_species <- function(time, state, pars) {
  A <- dim(pars$matrix)[2]  
  P <-dim(pars$matrix)[1]
  a <- state[1:A] ## species densities of animals
  p <- state[(A+1):(A+P)] ## species densities of plants
  s <- pars$sig ## species' trait standard deviations
  ma<-state[(A+P+1):(A+P+A)]
  mp<-state[(A+P+A+1):(A+P+A+P)]
  ## define g, where g[i] is the selection pressure on species i from growth
  alpha.a<-pars$Amatrix ## alpha matrix
  alpha.p<-pars$Pmatrix
  dt<-0.09
  aij<-bij<-matrix(0, nrow=A,ncol=P) 
  aji<-bji<-matrix(0, nrow=P,ncol=A) 
  Na<-muA<-matrix(0, nrow = time, ncol = A )
  Np<-muP<-matrix(0, nrow =time, ncol = P )
  Np[1,]<-state[(A+1):(A+P)]
  Na[1,]<-state[1:A]
  muA[1,]<-ma
  muP[1,]<-mp
  aj<-bj<-ai<-bi<-numeric()
  
  Np_tipping.point<-Na_tipping.point<-degree.plants<-degre.animals<-total.community.abundance<-trait_matching<-max_change<-
  tipping.point.index_P<-tipping.point.index_A<-numeric()
  mean.animal.abundance<-matrix(NA,nrow=length(pars$mut.strength), ncol=A)
  mean.plant.abundance<-matrix(NA,nrow=length(pars$mut.strength), ncol=P)
  rate_of_change_A<-matrix(NA,nrow=A,ncol=(length(pars$mut.strength)-1))
  rate_of_change_P<-matrix(NA,nrow=P,ncol=(length(pars$mut.strength)-1))
  
  for(x in 1:length(pars$mut.strength)){
  for (t in 1:(time-1)){
    for(r in 1:A){
      for(l in 1:P){
        #
        m.temp<-list(ma=muA[t,r],mp=muP[t,l])
        sigma1<-list(sa=s[r],sp=s[(A)+l])
        temp1<-gausquad.animals(m=m.temp,sigma=sigma1,w=pars$w,h=0.25,np=Np[t,l],na=Na[t,r],
                                mut.strength=pars$mut.strength[x], points=7
                                ,mat=pars$matrix[l,r],degree.animal = pars$dganimals[r],
                                interaction_type = pars$interaction_type)
        aij[r,l] <-temp1$G
        bij[r,l]<-temp1$B
        
      }
      ai[r]<-sum(aij[r,])
      bi[r]<-sum(bij[r,])
    }
    for(k in 1:P){
      for(m in 1:A){
        m2.temp<-list(ma=muA[t,m],mp=muP[t,k])
        sigma2<-list(sa=s[m],sp=s[(A+k)])
        temp2<-gausquad.plants(m=m2.temp,sigma=sigma2,w=pars$w,h=0.25,np=Np[t,k],na=Na[t,m],
                               mut.strength=pars$mut.strength[x],
                               points=7,mat=pars$matrix[k,m], degree.plant =pars$dgplants[k],
                               interaction_type = pars$interaction_type)
        aji[k,m] <-temp2$G
        bji[k,m]<-temp2$B
      }
      aj[k]<-sum(aji[k,])
      bj[k]<-sum(bji[k,])
    }
    #print(t)
    if(pars$noise == "multiplicative"){
      
      ba<-rnorm(A, pars$ba[1], sd = 0.15)
      bp<-rnorm(P, pars$bp[1], sd = 0.15)
      Na[t+1, ]<-Na[t,] + Na[t,]*(ba-alpha.a%*%Na[t,]+ai)*dt ## density eqs
      Np[t+1, ]<-Np[t,] + Np[t,]*(bp-alpha.p%*%Np[t,]+aj)*dt # density eqs
      muA[t+1, ]<-muA[t,] +pars$h2[1:A]*(bij%*%Np[t,])*dt #+ rnorm(A, 0,sd=0)*dt ## trait mean eqs
      muP[t+1, ]<- muP[t,]+pars$h2[(A+1):(A+P)]*(bji%*%Na[t,])*dt #+ rnorm(P, 0,sd=0)*dt ## trait mean eqs
      
      Na[t+1,which(Na[t+1,] < 0)]<-0
      Np[t+1,which(Np[t+1,] < 0)]<-0
      
    }else if (pars$noise == "additive") {
      
      Na[t+1, ]<-Na[t,] + Na[t,]*(pars$ba-alpha.a%*%Na[t,]+ai)*dt + rnorm(A, 0,sd=0.15)*Na[t,]*dt## density eqs
      Np[t+1, ]<-Np[t,] + Np[t,]*(pars$bp-alpha.p%*%Np[t,]+aj)*dt + rnorm(P, 0,sd=0.15)*Np[t,]*dt ## trait mean eqs
      muA[t+1, ]<-muA[t,] +pars$h2[1:A]*(bij%*%Np[t,])*dt #+ rnorm(A, 0,sd=0)*dt ## trait mean eqs
      muP[t+1, ]<- muP[t,]+pars$h2[(A+1):(A+P)]*(bji%*%Na[t,])*dt #+ rnorm(P, 0,sd=0)*dt ## trait mean eqs
      
      Na[t+1,which(Na[t+1,] < 0)]<-0
      Np[t+1,which(Np[t+1,] < 0)]<-0
      
      
    } else {
      
      Na[t+1, ]<-Na[t,] + Na[t,]*(pars$ba-alpha.a%*%Na[t,]+ai)*dt #+ rnorm(A, 0,sd=0.15)*Na[t,]*dt## density eqs
      Np[t+1, ]<-Np[t,] + Np[t,]*(pars$bp-alpha.p%*%Np[t,]+aj)*dt #+ rnorm(P, 0,sd=0.15)*Np[t,]*dt ## trait mean eqs
      muA[t+1, ]<-muA[t,] +pars$h2[1:A]*(bij%*%Np[t,])*dt #+ rnorm(A, 0,sd=0)*dt ## trait mean eqs
      muP[t+1, ]<- muP[t,]+pars$h2[(A+1):(A+P)]*(bji%*%Na[t,])*dt #+ rnorm(P, 0,sd=0)*dt ## trait mean eqs
      
      Na[t+1,which(Na[t+1,] < 0)]<-0
      Np[t+1,which(Np[t+1,] < 0)]<-0
      
    }
    
  } 
    
    for( i in 1:P){
      mean.plant.abundance[x,i]<-mean(Np[700:time,i],na.rm=T)
    }
    for(j in 1:A){
      mean.animal.abundance[x,j]<-mean(Na[700:time,j],na.rm=T)
    }
    
    total.community.abundance[x]<- sum( colMeans(Na[700:time,])) + sum(colMeans(Np[600:time,]))
    
    
    trait_matching[x]<- trait.matching(mA=muA[time,],
                   mP = muP[time,], 
                   adj.mat = pars$matrix, gamma=pars$w)
    #print(x)
  }
      
      rate.of.community.change<-abs(total.community.abundance[2:length(pars$mut.strength)]-
      total.community.abundance[1:(length(pars$mut.strength)-1)])
      
      #threshold at which a network collapses
      point.of.network.collapse<-max(pars$mut.strength[which( total.community.abundance < 1)])
     
      #inflection point of the network
      tipping.point.community.index<-max(which(rate.of.community.change > 45))
      
      #if there is no inflection point, or abprut change in consecutive network abundances, the network does not have a tipping point
      #and is termed as "NA"
      tipping.point.community.index[!is.finite(tipping.point.community.index)] <- NA
      
      abundance.at.tipping.point<-total.community.abundance[tipping.point.community.index-1]
     
      #whether collapse is abrupt or not determined by the rate of change in abundances
      if(max(rate.of.community.change) > 45){
        collapse = 1
      }else { collapse = 0
      }
      
      
      if (collapse == 0){
        
        for(l in 1:A){
          rate_of_change_A[l,]<- 0
          tipping.point.index_A[l]<- NA
          Na_tipping.point[l]<-NA
          degree.animals[l]<- sum(pars$matrix[,l])
        }
        for(j in 1:P){
          rate_of_change_P[j,]<- 0
          Np_tipping.point[j]<- NA
          degree.plants[j]<-sum(pars$matrix[j,])
        }
      } else if (collapse == 1){
          for(l in 1:A){
              rate_of_change_A[l,]<- abs(mean.animal.abundance[2:length(pars$mut.strength),l]-mean.animal.abundance[1:(length(pars$mut.strength)-1),l] )
              tipping.point.index_A[l]<- max(which(rate_of_change_A[l,] > 5))
              tipping.point.index_A[!is.finite(tipping.point.index_A[l])] <- NA
              Na_tipping.point[l]<- pars$mut.strength[tipping.point.index_A[l]]
              degree.animals[l]<- sum(pars$matrix[,l])
              #biomass_change<-abs(total.community.abundance[2:length(pars$mut.strength)] - total.community.abundance[1:(length(pars$mut.strength)-1)])
            }
        for(j in 1:P){
              rate_of_change_P[j,]<- abs(mean.plant.abundance[2:length(pars$mut.strength),j] -  mean.plant.abundance[1:(length(pars$mut.strength)-1),j])
              tipping.point.index_P[j]<- max(which(rate_of_change_P[j,] > 5))
              tipping.point.index_P[!is.finite(tipping.point.index_P)] <- NA
              Np_tipping.point[j]<- pars$mut.strength[tipping.point.index_P[j]]
              degree.plants[j]<-sum(pars$matrix[j,])
             }
      }
  
  
  
  
 output= list( Np_tipping.point=Np_tipping.point,
               Na_tipping.point=Na_tipping.point, 
               degree.animals=degree.animals,
               degree.plants=degree.plants,
               collapse=collapse,
               point.of.network.collapse=point.of.network.collapse,
               trait_matching=trait_matching,
               total.community.abundance=total.community.abundance
               )
  return(output)
}

Mcommunity_1 = function(iter, time, ...){
  set.seed(rnorm(1,as.numeric(Sys.time())-Sys.getpid(),10000)) 
  init = time
  replicate =try(eqs_species(time=init, ...))
  replicate$start = init
  replicate$iter = iter
  return(replicate)
}





## Organize simulation results into tidy table (code adapted from Barabas and D'Andrea 2016 Eco.Letts paper)
## Input:
## - sol: output produced by the function ode()
## - pars: list of parameters, with the following elements:
## Output:
## - a tibble with columns: time, species, n (density), m (trait mean),
##   sigma
organize_results <- function(sol, pars) {
  S <- length(pars$sigma) ## number of species
  A<-dim(pars$matrix)[2] # no. of animals
  P<-dim(pars$matrix)[1] # no. of plants
  temp<- sol %>% as.data.frame %>% as_tibble ## convert to tibble
  ## name the first column "time"
  temp<- temp %>% filter(time >= pars$cutoff.time)
  names(temp)[2:(A+1)] <- paste0("N_", 1:(A)) ## name abundance columns (n_k)
  names(temp)[1] <- "time"
  names(temp)[(A+2):(A+1+P)] <- paste0("P_", 1:P) ## name trait mean columns
  temp <- temp %>%
    gather("variable", "v", 2:ncol(temp)) %>% ## normalize the data
    separate(variable, c("type", "species"), sep="_") %>%
    #spread(type, v) %>% ## separate columns for animal densities n and plant densities m
    dplyr::select(time, type, species,v) %>% ## rearrange columns
    mutate(species=as.integer(species), Mut_strength=pars$mut.strength,
           Nestedness=pars$nestedness, Connectance=pars$C,
           theta=pars$theta,Web.name=pars$web.name) ## add params
  return(as_tibble(temp))
}


## Plot time series of densities, time series of trait values, and
## snapshot of the trait distributions at time = moment (code adapted from Barabas and D'Andrea 2016 Eco.Letts paper)
## Input:
## - dat: data generated by organize_results()
## - moment: time at which trait distribution should be plotted
## - limits: a vector of two entries (x_low, x_high) for the x-axis limits
## - res: number of evenly spaced sampling points along the trait axis
##               for the trait distribution plot
## Output:
## - a ggplot2 plot with three panels in one column: abundance time series,
##   trait value time seties, and snapshot of trait distribution
plot_all <- function(dat, moment=0, limits=c(-1, 1), res=1001) {
  plot_grid(plot_density(dat), ncol=1, align="hv") %>%
    return
}



## Plot species densities through time
## Input:
## - dat: data generated by organize_results()
## Output:
## - a ggplot2 plot
## used to produce figure 1.

## code adapted from Barabas and D'Andrea 2016 Eco.Letts paper
plot_density<- function(dat) {
  dat %>%
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species))) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    theme(legend.position="none") + facet_wrap(.~type) %>%
    return
}


#lay out function for multiple plots
lay_out = function(...) {    
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))    
  
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]], vp = grid::viewport(layout.pos.row = x[[i]][[2]], 
                                           layout.pos.col = x[[i]][[3]]))
  }
} 



#multiplot of ggplot2 figures with a common shared legend. Code taken from :https://rpubs.com/sjackman/grid_arrange_shared_legend
grid_arrange_shared_legend <- function(..., ncol, nrow, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol =ncol , nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}


# plots the density distribution  at a particular timepoint. This function was used to produce figure 1.
# Na: Abundance of animals at equilibrium
# Np: Abundance of plants at equilibrium
# m: mean traits at equilibrium
# sigma: variance of traits
# moment: mean
# limits: limits of the mean trait axis which in the study are -1,1

plot_snapshot <- function(Na, Np, m, sigma, moment=0, limits=c(-1, 1), res=1001) {
  Sa <- length(Na) ## number of species
  Sp <- length(Np)
  ma<- m[1:(Sa)]
  mp<- m[Sa+1:Sp]
  sigma_a <-sigma[1:(Sa)]
  sigma_p <- sigma[Sa+1:Sp]
  traitaxis <- seq(limits[1], limits[2], l=res) ## sampling the trait axis
  #snap <- dat %>% filter(time==moment) %>% select(-time) ## time = moment
  traits_a <- expand.grid(species=1:Sa, trait=traitaxis) %>% as_tibble ## trait table
  traits_p <- expand.grid(species=Sa+1:Sp, trait=traitaxis) %>% as_tibble ## trait table
  
  traits_a["density"] <- 0 ## add column for population densities
  traits_p["density"] <- 0
  
  for (i in 1:Sa) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_a$density[(traits_a$species==i)] <- Na[i]*
      dnorm(traits_a$trait[(traits_a$species==i)], ma[i], sigma_a[i]) ## times density
  }
  traits_a$density[traits_a$density<max(traits_a$density)/1e3] <- NA
  
  for (i in 1:Sp) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_p$density[(traits_p$species==(Sa+i))] <- Np[i]*dnorm(traits_p$trait[(traits_p$species==(Sa+i))], 
                                                                mp[i], sigma_p[i]) ## times density
  }
  traits_p$density[traits_p$density<max(traits_p$density)/1e3] <- NA
  

  traits<-data.frame(rbind(traits_a,traits_p), 
                     species_group=c(rep("Animals", nrow(traits_a)),
                                     rep("Plants", nrow(traits_p))))
  
  ggplot(traits) + ## generate plot
    geom_line(aes(x=trait, y=density, colour=factor(species)), na.rm=TRUE) +
    geom_ribbon(aes(x=trait, ymin=0, ymax=density, fill=factor(species)),
                alpha=0.15, colour=NA)+scale_fill_viridis_d()+
    facet_wrap(.~species_group, nrow = 2)+
    theme(legend.title = element_text(size = 14, face = "bold"), 
          legend.position = "right", panel.background = element_blank(), 
          axis.text = element_text(colour = "black", size = 14, face = "bold"), 
          axis.title = element_text(size = 14, face = "bold"), 
          legend.text = element_text(size = 14), legend.key = element_blank(),
          strip.text.x = element_text(size= 14, face ="bold"))+
    #geom_line(data=landscape, aes(x=trait, y=r), linetype="dashed",
    #         colour="darkred", alpha=0.5, na.rm=TRUE) +
    scale_x_continuous(name="trait value", limits=limits) +
    scale_y_continuous(name="density", limits=c(0, NA)) +
    theme(legend.position="none") %>%
    return }




#computes the raw NODF taken from Song et al 2017 J. Animal Ecology
#input: web = mutualistic network
#output: raw NODF of the given network
nestedness_NODF <- function(web){
  web[web > 0] = 1
  SA <- nrow(web)
  SP <- ncol(web)
  N <- t(web) %*% web
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SP,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SP)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele == 0] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)] = 1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n1 <- sum(nes)
  
  N <- web %*% t(web)
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SA,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SA)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele ==0 ] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)]=1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n2 <- sum(nes)
  out <- 2*(n1 + n2) / (SA*(SA-1)+SP*(SP-1))
  return(out)
}



# measures connectance of a web network
# web: interaction network
Connectance<-function(web)
{
  return(sum(web)/(ncol(web)*nrow(web)))}



# function for sampling competitive coefficients from random uniform distribution 
# competitive interactions  were scaled by the total number of species within a guild as Dakos & Bascompte 2014 PNAS.
# matrix: network of interactions which are 0 or 1. 
# strength: average competition strength

mat.comp<-function(matrix){
  Aspecies<- dim(matrix)[2]
  Plantspecies<- dim(matrix)[1]
  
  Amatrix<-matrix(runif(Aspecies^2, 0.0001, 0.001), nrow=Aspecies, ncol = Aspecies)/Aspecies #scaled by number of competitors within a guild
  diag(Amatrix)<-1 #intraspecific competition for animals
  Pmatrix<-matrix(runif(Plantspecies^2, 0.0001, 0.001), nrow=Plantspecies, ncol = Plantspecies)/Plantspecies ##scaled by number of competitors within a guild
  diag(Pmatrix)<-1 #intraspecific competion for plants
  
  out<-return(list(Amatrix=Amatrix,Pmatrix=Pmatrix))
  
}


# function for sampling competitive coefficients from exponential distribution 
# competitive interactions  were scaled by the total number of species within a guild as Dakos & Bascompte 2014 PNAS.
# matrix: network of interactions which are 0 or 1. 
# strength: average competition strength

mat.comp_exp<-function(matrix){
  Aspecies<- dim(matrix)[2]
  Plantspecies<- dim(matrix)[1]
  
  Amatrix<- matrix(rexp(Aspecies^2, 70),nrow=Aspecies,ncol=Aspecies)/Aspecies #scaling by no. of competitors within a guild
   
  diag(Amatrix)<-1
  Pmatrix<-matrix(rexp(Plantspecies^2, 70),nrow=Plantspecies,ncol=Plantspecies)/Plantspecies  #scaling by no. of competitors within a guild
  
  diag(Pmatrix)<-1
  
  out<-return(list(Amatrix=Amatrix,Pmatrix=Pmatrix))
  
}

# function for competition coefficients within a guild for feasibility analysis.
# competitive interactions  were scaled by the total number of species within a guild as Dakos & Bascompte 2014 PNAS.
# matrix: network of interactions which are 0 or 1. 
# strength: average competition strength
mat.comp_feasibility<-function(matrix,strength){
  Aspecies<- dim(matrix)[2]
  Plantspecies<- dim(matrix)[1]
  
  Amatrix<-strength*matrix(runif(Aspecies^2, 0.0001, 0.001), nrow=Aspecies, ncol = Aspecies)/Aspecies #scaled by number of competitors within a guild
  diag(Amatrix)<-1 #intraspecific competition for animals
  Pmatrix<-strength*matrix(runif(Plantspecies^2, 0.0001, 0.001), nrow=Plantspecies, ncol = Plantspecies)/Plantspecies ##scaled by number of competitors within a guild
  diag(Pmatrix)<-1 #intraspecific competion for plants
  
  out<-return(list(Amatrix=Amatrix,Pmatrix=Pmatrix))
  
}


#plot for feasibility analysis that takes data frame as input
feasibility_plot<-function(dat){
  akima.R<-with(dat,interp(dat$Strength_mutualism,log(dat$range_competition),dat$mean_feasibility,
                           xo=seq(min(dat$Strength_mutualism),max(dat$Strength_mutualism),length=20), 
                           yo=seq(min(log(dat$range_competition)),max(log(dat$range_competition)), length=20)))
  
  gdat1<-interp2xyz(akima.R, data.frame=TRUE)
  colnames(gdat1)<-c("mutualism_strength","range_competition","Feasibility")
  a<-ggplot(gdat1 , aes(x=mutualism_strength,y=range_competition,
                        z=Feasibility))+
    geom_raster(aes(fill=Feasibility),show.legend =T)+ 
    #geom_contour(aes(colour = ..level..),bins=1)+
    theme_cowplot()+
    theme(legend.title = element_text(size = 9, face = "bold"), 
          legend.text=element_text(size=rel(0.5)),
          legend.position = "bottom", panel.background = element_blank(), 
          axis.text = element_text(colour = "black", size = 9, face = "bold"), 
          axis.title = element_text(size = 9, face = "bold"), 
          legend.key = element_blank())+
    scale_fill_continuous(low = "#BFE1B0", high = "#137177") +
    labs(x = expression(paste("Mutualistic strength, ",gamma[0])),
         y = expression(paste(log(rho)))) +
    scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
  
  a
  
  
  return(a)
}


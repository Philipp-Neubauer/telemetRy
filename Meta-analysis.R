## Simulate fish movement on a linear coastline or river -----

mvmt_TS <- function(timesteps,N_fish,sex,mean_dist,gender_diff,sensitivity, domain,num_receivers){
  
  # pre-allocate vectors
  likes <- matrix(NA,length(domain),timesteps)
  pos <- matrix(NA,N_fish,timesteps)
  det_pos <- matrix(NA,N_fish,timesteps)

  # suming a linear coastline/river, the position at time t will depend on temperature (normal over domain)
  for (t in 1:timesteps)
  { cat(t,'\n')
    # simulate time-dependent temperatures suitedness along locations
    
    likes[,t]<- dlnorm(domain,log(t*100/timesteps),1)
    likes[,t]<- likes[,t]/max(likes[,t])
    maxlike = which.max(likes[,t])
    # if t=1 draw position at random along the same temperature like gradient, 
    #else the position will depend on time t plus a movement of random aplitude linked to temp dislike 
    # (move more and toward better if dislike strong, move more if you're a female)
    
    if (t==1){
      pos[,t] <- rlnorm(N_fish,log(t*100/timesteps),1)
    
      while (any(pos[,t]>max(domain)))
      {
        pos[pos[,t]>max(domain) ,t] <- rlnorm(length(pos[pos[,t]>max(domain) ,t]),log(t*100/timesteps),1)
      }
      
      } else {
        
        # weighted interpolation of temp
       likes_int<-(likes[floor(pos[,t-1])+1,t]*(pos[,t-1]-floor(pos[,t-1]))+likes[ceil(pos[,t-1])+1,t]*(ceil(pos[,t-1])-pos[,t-1]))/2
       signs <- sign(maxlike-pos[,t-1])
      
      pos[,t] <- pos[,t-1] + rnorm(N_fish,sensitivity*(likes[maxlike,t-1]-likes_int)*signs,mean_dist+sex*gender_diff) 
      
      # for simplicity, make sure fish stay in domain, only update fish that are outside doamin
      while (any(pos[,t]>max(domain) | pos[,t]<min(domain)))
      { ixxs <- which(pos[,t]>max(domain) | pos[,t]<min(domain))
        pos[ixxs,t] <- pos[ixxs,t-1] + rnorm(length(ixxs),sensitivity[ixxs]*(likes[maxlike,t-1]-likes_int[ixxs])*signs[ixxs],mean_dist+gender_diff[ixxs])
      }
     # detection only if fish swims by antenna detection radius define by a exponential distribution with rate r
       
      # can't properly use rbinom with a matrix, so use a workaround
      detections=matrix(rbinom(N_fish*length(seq(0,max(domain),num_receivers)),1,1-pexp(sapply(pos[,t],function(x){abs(seq(0,max(domain),num_receivers)-x)}),r)),length(seq(0,max(domain),num_receivers)),N_fish)
      det_pos[which(detections==1, arr.ind = T)[,2],t] = seq(0,max(domain),num_receivers)[which(detections==1, arr.ind = T)[,1]]
      
    }
  }
  list(det_pos=det_pos,pos=pos,likes=likes)
}

# Try the function ----

# N_fish = number of fish in 1 study
N_fish=25

# genders
sex = rbinom(N_fish,1,0.5)

# mean movement range
mean_dist = 2

# timesteps
timesteps =50

# sensitivity to temperature
sensitivity = rnorm(N_fish,0.7,0.2)

##-
## how much more (proportionally) do females move than males##
##-

gender_diff = mean_dist*rnorm(N_fish,0.5,0.3)

# number of receivers along a linear coastline/river
num_receivers = 10
num_receivers = num_receivers-1 # because I use it for itnervals there'll be as many as chosen...
# detection probability rate
r = 0.5 

# plot detection radius
plot(seq(0,10,length.out=1000),1-pexp(seq(0,10,length.out=1000),0.7),ylab='proba',xlab='distance',t='l',col=1,lwd=2)

domain = 0:100

mvmt <- mvmt_TS(timesteps,N_fish,sex,mean_dist,gender_diff,sensitivity, domain,num_receivers)
  
par(mfcol=c(3,1))
k=5
# detected positions
plot(mvmt$det_pos[k,],xlab='time',ylab='position')
# actual positions
lines(mvmt$pos[k,],lwd=2,col=(4-sex[k]))
k=3
plot(mvmt$det_pos[k,],xlab='time',ylab='position')
lines(mvmt$pos[k,],lwd=2,col=(4-sex[k]))
k=1
plot(mvmt$det_pos[k,],xlab='time',ylab='position')
lines(mvmt$pos[k,],lwd=2,col=(4-sex[k]))

# Prepare Winbugs data -------------------------------

Movement_model <- function(){
  for(i in 1:N_fish){
    
    # 'meta-analytic' beta and sex effects
    beta[i] ~ dnorm(mu_beta,tau_beta)
   
    
    for(t in 2:T){
      
      # beta is essentially a 'regression' for the effect of temperature
      means[i,t] <- pos[i,t-1] + beta[i]*temp_grad[i,t-1]
      
      # process model
      pos[i,t] ~ dnorm(means[i,t],tau[sex[i]]) %_% I(domain[1],domain[2])
      
      # temp 'preference' at pos[i] relative to best position
      rpos[i,t] <- round(pos[i,t])
      
      # now a gradient thanks to my little trick ! 
      temp_grad[i, t] <-  1 - temp[rpos[i, t], t]
      
      # detection model
      for (r in 1:n_rec){
        dist[i,t,r] <- abs(r-pos[i,t])
        # exponential detection model
        p[i,t,r] <- lam*exp(-lam*dist[i,t,r])/lam
        # measurement model
        pos_det[r, t,i] ~ dbern(p[i, t, r])
      }
    }
  }
   
  ##-- priors
  
  # sex specific variance
  meansig ~ dlnorm(0,0.0001)
  sig[1] <- meansig
  sig[2]<- sex_diff*meansig
  tau[1] <- 1/sig[1]
  tau[2] <- 1/sig[2]
  
  sex_diff ~ dlnorm(mu_sex, tau_sex)
  
  mu_sex~ dnorm(0.00000E+00, 1.00000E-04)
  
  mu_beta ~ dnorm(0,0.0001) 
  
  tau_beta ~ dgamma(0.01,0.01)
  tau_sex ~ dgamma(0.01,0.01)
  
  lam ~ dgamma(0.01,0.01)
  
}

write.model(Movement_model, 'Movement_model.bug')
#file.show('Movement_model.bug')

## prepare data ---

# we know where we tagged the fish....
pos1 <- mvmt$pos/10+1
pos1[,2:ncol(pos1)] <- NA

# need position of temp optimum
maxtemp <- apply(mvmt$likes,2,which.max)

# need receiver psoitions normalized to 1:num_receivers
recs=11

# detection positions

pos_det <- array(0,dim=c(recs,timesteps,N_fish))

for (i in 1:N_fish){
  for (j in 1:timesteps){
    pos_det[mvmt$det_pos[i,j]/10+1,j,i] <- 1}}

# initia temp_optimum gradient 
temp_grad <- matrix(NA,N_fish,timesteps)
signs <- sign(maxtemp[1]-pos1[,1])
tmp_pos <- likes[round(pos1[, 1]),1]
temp_grad[,1] <- signs * (1 - tmp_pos) 

# use trick to make MCMC faster
temp=mvmt$likes
for (t in 1:timesteps){
  temp[(maxtemp[t]+1):nrow(temp),t] <- 1+(1-temp[(maxtemp[t]+1):nrow(temp),t])}

mvmt_data <- list(pos_det=pos_det,sex=sex+1,pos=pos1,N_fish=N_fish,T=timesteps,n_rec=recs,domain=c(1,11),temp=temp,temp_grad=temp_grad)

## prepare inits ---------------

# initialize positions at actual positions....cheating but making sure that things work ok...
pos.init <- mvmt$pos/10+1
pos.init[,1] <- NA

sex_diff.init <- 1

mvmt_inits <- list(list(
  pos=pos.init,
  beta = rep(0,N_fish),
  mu_sex=1,
  sex_diff=sex_diff.init,
  mu_beta=0,
  tau_beta=1,
  tau_sex=1,
  lam=0.7,
  meansig=1))

# output parameters - 'pos_det' is not output since it is HUGE with jsut a moderate number of timesteps and iterations
params = c('pos','mu_beta','meansig','sex_diff','lam')

#main function to run winbugs-----------


mvmt.mcmc <- bugs(data=mvmt_data,
    inits=mvmt_inits,
    parameters.to.save=params,
    model.file='Movement_model.bug',
    n.chains=1,
    n.burnin=1000,
    n.iter=11000,
    bugs.directory = "C:\\Program Files (x86)\\WinBUGS14", # change this part - remove th x86 for non-64 bit systems
    working.directory=getwd())


# get chains for some parameters
post_pos <- mvmt.mcmc$sims.list$pos
mu_beta <- mvmt.mcmc$sims.list$mu_beta
sex_diffs <- mvmt.mcmc$sims.list$sex_diff
lam <- mvmt.mcmc$sims.list$lam

# plot posterior distributions
par(mfcol=c(1,1))
hist(mu_beta)
hist(sex_diffs)
hist(meansig)
hist(lam)

# p value of temp effect
mean(mu_beta>0)

# p value of sex diff
mean(sex_diffs>1)

# plot inferred, smoothed trajectories

plot_funct <- function (k) {
  
  # detected positions
  plot(mvmt$pos[k,],lwd=2,col=(4-sex[k]),xlab='time',ylab='position',t='l')
  
  for (i in seq(1,1000,10))
  {lines((post_pos[i,k,]-1)*10,lwd=1,col=2)}
  
  # detected positions
  points(mvmt$det_pos[k,])
  # actual positions
  lines(mvmt$pos[k,],lwd=2,col=(4-sex[k]))
}

plot_funct(1)


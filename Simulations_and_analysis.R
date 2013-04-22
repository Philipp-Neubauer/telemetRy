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
      
      # can't properly use rbinom with a matrix, so use a workaround
      detections=matrix(rbinom(N_fish*length(seq(0,max(domain),num_receivers)),1,1-pexp(sapply(pos[,t],function(x){abs(seq(0,max(domain),num_receivers)-x)}),r)),length(seq(0,max(domain),num_receivers)),N_fish)
      det_pos[which(detections==1, arr.ind = T)[,2],t] = seq(min(domain),max(domain),num_receivers)[which(detections==1, arr.ind = T)[,1]]
      
    } else {
      
      # weighted interpolation of temp
      likes_int<-(likes[floor(pos[,t-1])+1,t]*(pos[,t-1]-floor(pos[,t-1]))+likes[ceiling(pos[,t-1])+1,t]*(ceiling(pos[,t-1])-pos[,t-1]))/2
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
      det_pos[which(detections==1, arr.ind = T)[,2],t] = seq(min(domain),max(domain),num_receivers)[which(detections==1, arr.ind = T)[,1]]
      
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
num_receivers = num_receivers # because I use it for itnervals there'll be as many as chosen...
# detection probability
r = 0.7

# plot detection radius
plot(seq(0,10,length.out=1000),exp(-r * seq(0,10,length.out=1000)),ylab='proba',xlab='distance',t='l',col=1,lwd=2)

domain = 1:100

mvmt <- mvmt_TS(timesteps,N_fish,sex,mean_dist,gender_diff,sensitivity, domain,num_receivers)


# Prepare data -------------------------------

# we know where we tagged the fish....
pos1 <- mvmt$pos
pos1[,2:ncol(pos1)] <- NA

#need receiver psoitions normalized to 1:num_receivers
recs=10

# detection positions

pos_det <- array(0,dim=c(recs,timesteps,N_fish))

for (i in 1:N_fish){
  for (j in 1:timesteps){
    pos_det[(mvmt$det_pos[i,j]+9)/10,j,i] <- 1}}

# need position of temp optimum
maxtemp <- apply(mvmt$likes,2,which.max)

# initial temp_optimum gradient 
temp_grad <- matrix(NA,N_fish,timesteps)
signs <- sign(maxtemp[1]-pos1[,1])
tmp_pos <- mvmt$likes[ceiling(pos1[, 1]),1]
temp_grad[,1] <- signs * (1 - tmp_pos) 

# use trick to make MCMC faster
temp=mvmt$likes
for (t in 1:timesteps){
  temp[(maxtemp[t]+1):nrow(temp),t] <- 1+(1-temp[(maxtemp[t]+1):nrow(temp),t])}

# plot examples - open circles are no detections, filled cicles are detections, lines are simualted movements
# colors show evolving temp field

k=5 # fish number
par(mfcol=c(1,1))
# detected positions
image(x=1:50-1,y=1:100,t(temp),xlab='time',ylab='position (river kilometer)')
points(expand.grid(1:50,seq(1,100,10)))
points(mvmt$det_pos[k,],pch=16)
# actual positions
lines(mvmt$pos[k,],lwd=2,col=(4-sex[k]))

# plot movement for all simualted fish
for (k in 1:N_fish){
  #points(mvmt$det_pos[k,],xlab='time',ylab='position')
  lines(mvmt$pos[k,],lwd=2,col=(4-sex[k]))
}


mvmt_data <- list(pos_det=pos_det,sex=sex+1,pos=pos1,N_fish=N_fish,rec_pos = seq(min(domain),max(domain),num_receivers),T=timesteps,n_rec=recs,domain=c(min(domain),max(domain)),temp=temp,temp_grad=temp_grad)

## prepare inits ---------------

# initialize positions at actual positions....cheating but making sure that things work ok...
pos.init <- mvmt$pos
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


require(rjags)

# output parameters - 'pos_det' is not output since it is HUGE with just a moderate number of timesteps and iterations
mvmt_bugs = jags.model('Movement_model.bug',data=mvmt_data,inits=mvmt_inits,n.adapt=100)

params = c('pos','mu_beta','meansig','sex_diff','lam')

mvmt_mcmc <- coda.samples(mvmt_bugs,variable.names=params,n.iter=10000,thin=10)

# get chains for some parameters
# predicted positions
post_pos <- mvmt_mcmc[[1]][,grep('pos',colnames((mvmt_mcmc[[1]])))]
# population temp effect
mu_beta <- mvmt_mcmc[[1]][,grep('mu_beta',colnames((mvmt_mcmc[[1]])))]
# sex differences in movement
sex_diffs <- mvmt_mcmc[[1]][,grep('sex_diff',colnames((mvmt_mcmc[[1]])))]
# detection function parameter (exponential)
lam <- mvmt_mcmc[[1]][,grep('lam',colnames((mvmt_mcmc[[1]])))]

# plot posterior distributions
par(mfcol=c(1,1))
hist(mu_beta,xlab='temperature effect',main='')
hist(sex_diffs,xlab='variance ratio female/male',main='')
hist(lam,xlab=expression(lambda),main='')

# p value of temp effect
mean(mu_beta>0)

# p value of sex diff
mean(sex_diffs>1)

# plot inferred, smoothed trajectories
pos_ar <- array(post_pos,c(1000,N_fish,timesteps))
plot_funct <- function (k) {
  
  # detected positions
  plot(mvmt$pos[k,],lwd=2,col=(4-sex[k]),xlab='time',ylab='position',t='l')
  
  for (i in seq(1,1000,10))
  {lines((pos_ar[i,k,]),lwd=1,col=2)}
  
  # detected positions
  points(mvmt$det_pos[k,])
  # actual positions
  lines(mvmt$pos[k,],lwd=2,col=(4-sex[k]))
}

plot_funct(8)




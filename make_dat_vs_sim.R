#make_dat_vs_sim.R
#construct figure of histograms and simulated data plots from model

set.seed(12345)

#helper function
logistic <- function(x){1/(1+exp(-x))}

#load up data
load("ppv_results") #load data output
source("all_ppv_data.R") #raw data

#some useful reference variables
vnames=colnames(xx)
pnames=colnames(pp)
numcat=4
nummonk=8

#the plan is to do posterior predictive checking:
#draw a whole bunch of samples from Bayesian posteriors of the *actual data*
#and compare these to *fictitious* session data drawn de novo from the day-to-day 
#distributions inferred from the data
#in other words, show that the session choice curves we see are typical for our model

#ready the svg!
svg(file='dat_vs_sim.svg')
par(mfcol=c(3,2))

#first, we want to sample the output of the model:
sel.inds=grepl("^v",vnames) #which columns contain value data?
sel.rows=sample(dim(xx)[1],sum(sel.inds)) #one row for each variable -- sampling of real data
pse=1000*xx[cbind(sel.rows,which(sel.inds))] #bind into PSEs
hp=hist(pse,breaks=100,xlab="Image value (ms juice)",main="Histogram of image value (data)")

#now sample scale parameters for choice curves
sel.inds=grepl("choice.scale",vnames)
sel.rows=sample(dim(xx)[1],sum(sel.inds)) #one row for each variable -- sampling of real data
scale=1000*xx[cbind(sel.rows,which(sel.inds))] #one scale for each session
scale=rep(scale,times=numcat) #each session has multiple categories, v[session,category]
hs=hist(scale[scale<200],breaks=100,xlab="Curve width (ms juice)",main="Histogram of curve width (data)")

#finally, use these data to plot choice curves
for (ind in 1:length(pse)){
  doadd=ifelse(ind==1,FALSE,TRUE)
  curve(logistic((x-pse[ind])/scale[ind]),from=-120,to=120,add=doadd,xlab="",ylab="",ylim=c(0,1),
        col=rgb(0,0,0,alpha=0.1))
}
title(main="Choice curves for all sessions (data)",xlab="Difference in juice value (ms)",
      ylab="Percent choose image")


#now select random pses and scales
#numgrab=ceiling(numsess/nummonk) #each line contains both monks, one pse per category
#for a better estimate, sample monks in proportion to actual number of sessions performed!
pse=c()
scale=c()
for (ind in 1:nummonk){
  #first, grab pses
  numgrab=length(which(sessvec==ind)) #number of sessions for this monk
  sel.rows=sample(dim(pp)[1],numgrab) #select random rows
  sel.inds=grepl(paste("ppred.v\\[",ind,sep=""),pnames) #get the right monk
  this.pse=c(as.matrix(1000*pp[sel.rows,sel.inds])) #sample pses
  pse=c(pse,this.pse) #concatenate pses
  
  #now, grab scales for sessions
  numgrab=numcat*length(which(sessvec==ind)) #one per category per session for this monk
  sel.rows=sample(dim(pp)[1],numgrab)
  sel.inds=grepl(paste("ppred.choice.scale\\[",ind,sep=""),pnames) #get the right monk
  this.scale=c(as.matrix(1000*pp[sel.rows,sel.inds]))
  scale=c(scale,this.scale)
}
#now plot stuff
hist(pse[(pse<max(hp$breaks))&(pse>min(hp$breaks))],
     breaks=hp$breaks,xlab="pse",main="Histogram of pse (simulated)")
hist(scale[scale<max(hs$breaks)],breaks=union(0,hs$breaks),
     xlab="width",main="Histogram of curve width (simulated)")
for (ind in 1:length(pse)){
  doadd=ifelse(ind==1,FALSE,TRUE)
  curve(logistic((x-pse[ind])/scale[ind]),from=-120,to=120,add=doadd,xlab="",ylab="",ylim=c(0,1),
        col=rgb(0,0,0,alpha=0.1))
}
title(main="Choice curves for all sessions (simulated)",xlab="Difference in juice value (ms)",
      ylab="Percent choose image")
par(mfrow=c(1,1))
dev.off()
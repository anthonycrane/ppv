#dat_vs_sim
#histograms and simulated data plots from model
logistic <- function(x){1/(1+exp(-x))}
#first, plot some fitted choice curves from data
pdf(file='dat_vs_sim.pdf')
par(mfcol=c(3,2))
sel.inds=grepl("^v",vnames)
sel.rows=sample(dim(xx)[1],sum(sel.inds)) #one row for each variable -- sampling of real data
pse=1000*xx[cbind(sel.rows,which(sel.inds))]
#pse=1000*qq[sel.inds,3] #get all median values for pses -- but this underestimates variance
hp=hist(pse,breaks=100,xlab="pse",main="Histogram of pse (data)")
sel.inds=grepl("choice.scale",vnames)
sel.rows=sample(dim(xx)[1],sum(sel.inds)) #one row for each variable -- sampling of real data
scale=1000*xx[cbind(sel.rows,which(sel.inds))] #one scale for each session
scale=rep(scale,times=numcat) #each session has multiple categories, v[session,category]
#scale=1000*qq[sel.inds,3] #median values for pse scale -- but this underestimates variance
hs=hist(scale[scale<200],breaks=100,xlab="width",main="Histogram of curve width (data)")
for (ind in 1:length(pse)){
  doadd=ifelse(ind==1,FALSE,TRUE)
  curve(logistic((x-pse[ind])/scale[ind]),from=-120,to=120,add=doadd,xlab="",ylab="",ylim=c(0,1))
}
title(main="Choice curves for all sessions (data)",xlab="difference in juice value (ms)",
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
  curve(logistic((x-pse[ind])/scale[ind]),from=-120,to=120,add=doadd,xlab="",ylab="",ylim=c(0,1))
}
title(main="Choice curves for all sessions(simulated)",xlab="difference in juice value (ms)",
      ylab="Percent choose image")
par(mfrow=c(1,1))
dev.off()
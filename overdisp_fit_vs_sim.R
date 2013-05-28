#overdisp_fit_vs_sim.R
#compare overdispersion as fit to data and overdispersion as simulated in model
par(mfrow=c(2,1))
#first, do data
sel.inds=grepl("resid.lp",vnames)
sel.rows=sample(dim(xx)[1],sum(sel.inds)) #one row for each variable -- sampling of real data
omega=1000*xx[cbind(sel.rows,which(sel.inds))]
#omega=1000*qq[sel.inds,3] #by taking only medians, this underestimates variance
hist(omega,breaks=100,main="Overdispersion noise (data)",xlab=expression(lp-hat(lp)))
orng=range(omega)
#now get simulated data
numgrab=ceiling(length(omega)/(10*nummonk)) #each row in data contains 10 draws from omega for each monk
#again, would be better to sample each monk in proportion to actual session numbers
sel.rows=sample(dim(pp)[1],numgrab)
sel.inds=grepl("ppred.resid.lp",pnames)
omega=c(as.matrix(1000*pp[sel.rows,sel.inds]))
hist(omega,breaks=100,main="Overdispersion noise (simulated)",xlab=expression(lp-hat(lp)),xlim=orng)
par(mfrow=c(1,1))
#plots_ppv.R
#exploratory plots of jags output of karli's ppv data

# osname=as.character(Sys.info()['sysname'])
# dirpath=switch(osname,
#                Windows="C:/Documents and Settings/pearson.CCN/My Documents/My Dropbox",
#                Darwin="~/Dropbox",
#                Linux="~/Dropbox")
# setwd(paste(dirpath,"/frbayes/ppv",sep=""))
library(rjags)
library(vioplot)
library(MASS)

load("ppv_trend_results") #load data output
source("all_ppv_data.R") #raw data
source("ppv_fitting.R")  #useful plotting functions
source("conventional_fits.R") #more conventional day-to-day fits of data

#set up some useful variables
vnames=colnames(xx)
pnames=colnames(pp)
cnames=c("gray","peri","dom","sub")
mnames=c("Ernie","Oskar","Otto","Dart","Sherry","Cajal","Broome","Niko")

###### fit diagnostics ############
#counts level:
#R^2
sel.inds=grepl("resid.N",vnames) #draws from distributions for residuals
xsel=as.matrix(xx[,sel.inds])
SS.resid=mean(apply(xsel,1,var)) #expected value of the variance across draws
SS.tot=var(Nimg)
R2.N=1-SS.resid/SS.tot
#pooling
pooling.N=1-var(apply(xsel,2,mean))/SS.resid

#R^2: probabilities level
sel.inds=grepl("resid.lp",vnames)
xsel=as.matrix(xx[,sel.inds])
SS.resid=mean(apply(xsel,1,var)) #expected value of the variance across draws
SS.pool=var(apply(xsel,2,mean))
sel.inds=grepl("linpred",vnames)
xsel=as.matrix(xx[,sel.inds])
SS.tot=mean(apply(xsel,1,var)) #expected value of the variance across draws
R2.lp=1-SS.resid/SS.tot
#pooling
pooling.lp=1-SS.pool/SS.resid

#R^2: session level
#note: only one dispersion parameter per monk for all categories, so don't break them out separately
sel.inds=grepl("resid.v",vnames)
xsel=1000*xx[,sel.inds] #scale in ms
SS.resid=mean(apply(xsel,1,var))
SS.pool=var(apply(xsel,2,mean))
sel.inds=grepl("v",vnames)
xsel=1000*xx[,sel.inds] #scale in ms
SS.tot=mean(apply(xsel,1,var))
R2.v=1-SS.resid/SS.tot
#pooling
pooling.v=1-SS.pool/SS.resid
##############################################

#plot session-by-session values for a given category
catnum=2
#plot pooled multi-level fits
plotsbys(catnum,xx,sessvec,datevec)
#mean.eff=mean(quants[,2])
#mean.spread=mean(quants[,3]-quants[,1])
#hist(1000*qq[v.inds,3],breaks=20) #histogram of day-to-day values (check for unimodality)

#plot standard fits alongside
#requires having fit.sets from explore_ppv_data.R already run
sc.mat=cbind(rep(1:numsess,each=numcat),rep(1:numcat,times=numsess)) #unique combinations of session and category
plotcat(catnum,fit.sets,which.fit=1,add.plot=T,jitt=0.25,colset=rep("gray",length(fit.sets)))

#plot select categories and sessions
mnum=6
this.sessions=which(sessvec==mnum)
this.cats=3:4
plot.selections(this.cats,1000*qq,this.sessions)

#plot curve widths
sel.inds=grepl("choice.scale",vnames)
xsel=1000*xx[,sel.inds]
quants=array(dim=c(dim(xsel)[2],3))
for (ind in 1:dim(xsel)[2]){
  quants[ind,]=quantile(xsel[,ind],c(0.025,0.5,0.975))
}
plot(quants[,2],pch=16,xlab="Session",ylab="Choice curve width (ms juice)")
segments(x0=1:dim(xsel)[2],y0=quants[,1],y1=quants[,3])

#plot curve distribution parameters
sel.inds=grepl("^t\\.",vnames)
xsel=qq[sel.inds,]
nn=rownames(xsel)
xsel=diag(c(1,1000,1000))%*%xsel #rescale scale and location parameters to ms juice
rownames(xsel)=nn
plotquants(xsel)

#ANOVA table for superpopulation variances: image value
sel.inds=grepl("sess.std",vnames)
super.sd=1000*qq[sel.inds,]
plotquants(super.sd,varlabs=mnames) 
title(main="ANOVA: Session-by session image value standard deviations",xlab='ms juice')

#ANOVA table for superpopulation variances: choice variability
sel.inds=grepl("omega.std",vnames)
super.sd=1000*qq[sel.inds,]
plotquants(super.sd,varlabs=mnames) 
title(main="ANOVA: Session-by session performance standard deviations",xlab='logits')

#ANOVA table for finite population variances
fin.sd=c()
qvec=c(0.025,0.25,0.5,0.75,0.975)
for (ind in 1:nummonk){
  sel.inds=grepl("resid.v",vnames)
  xsel=xx[,sel.inds]
  xsel=xsel[,which(sessvec==ind)] #get only sessions for this monk
  post.var=apply(xsel,1,sd) #posterior for finite population variance
  fin.sd=rbind(fin.sd,quantile(post.var,qvec)) #add quantiles of this
  rownames(fin.sd)[dim(fin.sd)[1]]=paste("v.std.fin[",ind,"]",sep="")
}
for (ind in 1:nummonk){
  sel.inds=grepl("resid.lp",vnames)
  xsel=xx[,sel.inds]
  xsel=xsel[,which(sessvec[session]==ind)] #get only sessions for this monk
  post.var=apply(xsel,1,sd) #posterior for finite population variance
  fin.sd=rbind(fin.sd,quantile(post.var,qvec)) #add quantiles of this
  rownames(fin.sd)[dim(fin.sd)[1]]=paste("omega.std.fin[",ind,"]",sep="")
}
fin.sd=1000*fin.sd #convert to ms juice scale
plotquants(fin.sd)

#compare
plotquants(rbind(super.sd,fin.sd))

#plot posterior values by monk
mnum=1
V.inds=grepl(paste("V\\[",mnum,sep=""),vnames)
xsel=1000*xx[,V.inds]
vioplot(xsel[,1],xsel[,2],xsel[,3],xsel[,4],names=cnames,col="gray")
abline(h=0)
title(paste("Posterior means for image value\n Monkey ",mnames[mnum],sep=""))

#all monks together
V.inds=grepl("V",vnames)
xsel=1000*xx[,V.inds]
colmap=rainbow(nummonk)
plot.new()
plot.window(xlim=c(0,sum(V.inds)+1),ylim=range(xsel))
par(lab=c(9,5,7))
for (ind in 1:nummonk){
  skip=nummonk
  offset=(ind-1)*numcat
  vioplot(xsel[,ind],xsel[,ind+skip],xsel[,ind+2*skip],xsel[,ind+3*skip],
          at=offset+(1:numcat),col=colmap[ind],add=T)
  text(x=2.5+offset,y=-100,mnames[ind],col=colmap[ind])
}
axis(side=1,at=1:sum(V.inds),labels=rep(cnames,nummonk))
axis(side=2)
abline(h=0)
title(main="Posteriors for mean image value", ylab="Juice value (ms)")

#all monks together -- time trend
V.inds=grepl("vslope",vnames)
xsel=1000*xx[,V.inds]
colmap=rainbow(nummonk)
plot.new()
plot.window(xlim=c(0,sum(V.inds)+1),ylim=range(xsel))
par(lab=c(9,5,7))
for (ind in 1:nummonk){
  skip=nummonk
  offset=(ind-1)*numcat
  vioplot(xsel[,ind],xsel[,ind+skip],xsel[,ind+2*skip],xsel[,ind+3*skip],
          at=offset+(1:numcat),col=colmap[ind],add=T)
  text(x=2.5+offset,y=-100,mnames[ind],col=colmap[ind])
}
axis(side=1,at=1:sum(V.inds),labels=rep(cnames,nummonk))
axis(side=2)
abline(h=0)
title(main="Posteriors for mean image value slope", ylab="Juice value (ms/session)")

#correlations among PSEs
mnum=6
V.inds=grepl(paste("V\\[",mnum,sep=""),vnames)
xsel=1000*xx[,V.inds]
colnames(xsel)=c("gray","peri","dom","sub")
pairs(~gray+peri+dom+sub,data=xsel)
cor(xsel)

#better: correlations among day-to-day PSEs
vtemp=list(c(),c(),c(),c())
for (catnum in 1:4){
  v.inds=grepl(paste("^v\\[\\d*,",catnum,sep=""),vnames)
  xsel=1000*xx[,v.inds]
  vtemp[[catnum]]=as.vector(as.matrix(xsel))
}
xsel=as.data.frame(vtemp)
colnames(xsel)=c("gray","peri","dom","sub")
pairs(~gray+peri+dom+sub,data=xsel,
      subset=sample(dim(xsel)[1],1e4,replace=F)) #take only a small subsample
cor(xsel,method="spearman")

############# predictive data sets #################
#compare various stats for fitted vs simulated data from model
source('dat_vs_sim.R')

#compare overdispersions in data vs simulated
source('overdisp_fit_vs_sim.R')

#compare variance around choice curves
dvaxis=c(-40,-20,0,20,40)
numresid=length(dvaxis)
numcurves=3
par(mfrow=c(2,numcurves))
#sample sessions and categories
sel.sess=sample(numsess,numcurves)
sel.cat=sample(numcat,numcurves)
for (ind in 1:numcurves){sessfit(sel.sess[ind],sel.cat[ind])}
#now, get simulations
sel.inds=grepl("^ppred.v",pnames)
xsel=1000*pp[,sel.inds] #big list of all pses 
#works because monks look pretty much the same, otherwise need to match pses with monk and draw
#scales and omegas from relevant distributions
sel.rows=sample(dim(xsel)[1],numcurves)
sel.cols=sample(dim(xsel)[2],length(sel.rows))
pse=xsel[cbind(sel.rows,sel.cols)]
#now get scales
sel.inds=grepl("ppred.choice.scale",pnames)
xsel=c(as.matrix(1000*pp[,sel.inds])) #big list of all scales 
scale=xsel[sel.rows]
#finally, get omegas
sel.inds=grepl("ppred.resid.lp",pnames)
xsel=pp[,sel.inds] #matrix of residuals (=omega)
omega=xsel[sel.rows,1:numresid]
#make linear predictor
pse=rep(pse,numresid)
dim(pse)=c(numcurves,numresid)
dvaxis=rep(dvaxis,each=numcurves)
dim(dvaxis)=c(numcurves,numresid)
scale=rep(scale,numresid)
dim(scale)=c(numcurves,numresid)
lp=(dvaxis+pse)/scale
pmat=logistic(lp+omega)
for(ind in 1:dim(pmat)[1]){
  plot(dvaxis[ind,],pmat[ind,],main="simulated session",xlab="Percent choose image",ylab="dv")
  curve(logistic((x-pse[ind,1])/scale[ind,1]),from=min(dvaxis),to=max(dvaxis),add=T)
}
par(mfrow=c(1,1))
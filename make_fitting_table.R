#make_fitting_table.R
#make table of fit statistics given JAGS output 

load("ppv_results") #load data output
source("all_ppv_data.R") #raw data

#set up some useful variables
vnames=colnames(xx)
pnames=colnames(pp)

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

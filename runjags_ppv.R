#runjags_frbayes.R
#run jags for a Bayesian hierarchical firing rate regression model
#on ppv data

ptm=proc.time() #get start time

# osname=as.character(Sys.info()['sysname'])
# dirpath=switch(osname,
#                Windows="C:/Documents and Settings/pearson.CCN/My Documents/My Dropbox",
#                Darwin="~/Dropbox")
# setwd(paste(dirpath,"/frbayes/ppv",sep=""))
library(rjags)

modstr="ppv6"

bugstr=paste(modstr,".bug",sep="")

d <- read.jagsdata("all_ppv_data.R")
load.module("glm")
m <- jags.model(bugstr, d, n.chains=5,n.adapt=1000) #,inits=initfun)

update(m,10000)

#x <- coda.samples(m, c("group.mean","group.cov","V","v","choice.std","choice.dof",
#                       "Q","monk.sig"), n.iter=5000,thin=2)
#x <- coda.samples(m, c("V","v","choice.std","choice.dof"), n.iter=20000,thin=10)
#x <- coda.samples(m, c("V","v","choice.scale","omega.std","sess.std", #ppv5
#                       "resid.v","resid.lp","resid.N","linpred"), n.iter=20000,thin=10)
x <- coda.samples(m, c("V","v","choice.scale","omega.std","sess.std", #ppv5
                       "resid.v","resid.lp","resid.N","linpred",
                       "t.loc","t.scale","t.df"), n.iter=20000,thin=20)

#now, do diagnostics
jagssum=summary(x)
ss<-jagssum$statistics
qq<-jagssum$quantiles
rejectionRate(x)
ff=effectiveSize(x)

# calculate Rhat convergence stat
getstat <- function(x,fun){apply(as.data.frame(x),2,fun)} #calculate standard deviation of sample chains
nn=dim(x[[1]])[1] #number of samples
chain_var = lapply(x,function(x)getstat(x,var))
within_chain_mat = do.call(rbind,chain_var)
within_chain_var = apply(within_chain_mat,2,mean)
chain_mean = lapply(x,function(x){getstat(x,mean)})
between_chain_mat=do.call(rbind,chain_mean)
between_chain_var=nn*apply(between_chain_mat,2,var)
totvar=(1/nn)*((nn-1)*within_chain_var+between_chain_var)
Rhat=sqrt(totvar/within_chain_var)

# bind all samples into one data frame
xx=as.data.frame(do.call(rbind,x))

#now draw some fake (posterior predictive) data
p <- coda.samples(m, c("ppred.v","ppred.choice.scale","ppred.resid.lp"), n.iter=5000,thin=10)
pp=as.data.frame(do.call(rbind,p))

dumplist=c("xx","ss","qq","ff","pp","Rhat")
fname="ppv_results"
save(list=dumplist,file=fname)


#how long did this take?
endtime=proc.time()-ptm
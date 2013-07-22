#explore_ppv_data.R
#look at pse fitting and diagnostics for ppv behavioral data
library("MASS")
library("arm")

source("all_ppv_data.R")
source("ppv_fitting.R")

cnames=c("gray","peri","dom","sub")

set.seed(12345)

#first, combine all data, plot choice curve of everything
Ntot.tab=table(Ntot,dv)
Ntot.vals=as.integer(rownames(Ntot.tab))
Ntot.comb=Ntot.vals%*%Ntot.tab
Nimg.tab=table(Nimg,dv)
Nimg.vals=as.integer(rownames(Nimg.tab))
Nimg.comb=Nimg.vals%*%Nimg.tab
plot(sort(unique(dv)),Nimg.comb/Ntot.comb)

#same thing for selected picture category
catnum=4
Ntot.tab=table(Ntot[piccat==catnum],dv[piccat==catnum])
Ntot.vals=as.integer(rownames(Ntot.tab))
Ntot.comb=Ntot.vals%*%Ntot.tab
Nimg.tab=table(Nimg[piccat==catnum],dv[piccat==catnum])
Nimg.vals=as.integer(rownames(Nimg.tab))
Nimg.comb=Nimg.vals%*%Nimg.tab
udv=as.numeric(colnames(Nimg.comb))
plot(udv,Nimg.comb/Ntot.comb)
#try fitting a curve
beta0=c(0,1,10)
lb=c(-20,0.1,0.01)
ub=c(20,100,100)
optimfn <- function(beta){-LLrobit(beta,udv,Nimg.comb,Ntot.comb)}
fit=optim(beta0,optimfn,method="L-BFGS-B",lower=lb,upper=ub)
curve(pt(fit$par[2]*(x+fit$par[1]),fit$par[3]),from=min(udv),to=max(udv),add=TRUE)

#plot fits for individual session and category
sessplot(sess=44,catnum=1)

#fit an overdispersed binomial model to a given session and category
sessplot(4,4,model='odlogit')

#fit all sessions
fit.sets = list()
fit.objs = list()
sc.mat=cbind(rep(1:numsess,each=numcat),rep(1:numcat,times=numsess)) #unique combinations of session and category
modlist = c('logit','odlogit')
for (ind in 1:numsess){
  for (ind2 in 1:numcat){
    for (model in modlist){
      fit = sessfit(ind,ind2,model)
      if (!is.na(fit)){
        pf = process.fit(fit)
      }
      fit.objs[[model]][[(ind-1)*numcat+ind2]] = fit
      fit.sets[[model]][[(ind-1)*numcat+ind2]] = pf
    }
  }
}
dump(c("fit.sets","fit.objs"),file="conventional_fits.R")

#plot estimates of image value
catnum=2
aa=plotcat(catnum,fit.sets,1)
bb=plotcat(catnum,fit.sets,2,1)

#histogram within category
catnum=2
sel=(sc.mat[,2]==catnum)
imval=as.vector(sapply(fit.sets[[1]],function(x){x$coef[1]},simplify=T))
hist(imval[sel],breaks=c(-1e6,seq(from=-40,to=40,by=5),1e6),xlim=c(-40,40))
sd(imval[sel])
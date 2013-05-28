#make_figure_sessfit.R
#construct a figure illustrating:
#1) well-behaved behavioral data with a standard logistic fit
#2) ill-behaved data with an ill-behaved logistic fit
#3) partial pooling fits to both of these
#based on explore_ppv_data.R and ppv_fitting.R

#things we need
library(arm)

#load up data
source("all_ppv_data.R")
load('ppv_results')

logistic <- function(x){1/(1+exp(-x))}

findmap <- function(x,N){
  #find maximum a posteriori value of a given variable given x, a vector of samples
  #from its distribution
  #N is the number of points to use, and should be a power of 2 (as per help from ?density)
  
  dd = density(x,n=N)
  
  max_ind = which.max(dd$y)
  
  mapest = dd$x[max_ind]
  
  return(mapest)
  
}

#which monk is which
mvec=c('E','Os','Ot','D','S','C','B','N')

#which category
cnames=c("Neutral","Female","Dominant Male","Subordinate Male")

#this function fits a logit model to a given session and category
sessfit <- function(sess,catnum){
  sel=(piccat==catnum)&(session==sess)
  if (any(sel)){
    Ntot.tab=table(Ntot[sel],dv[sel])
    Ntot.vals=as.integer(rownames(Ntot.tab))
    Ntot.comb=Ntot.vals%*%Ntot.tab
    Nimg.tab=table(Nimg[sel],dv[sel])
    Nimg.vals=as.integer(rownames(Nimg.tab))
    Nimg.comb=Nimg.vals%*%Nimg.tab
    pp=Nimg.comb/Ntot.comb
    qq=1-pp
    sd=sqrt(pp*qq/Ntot.comb)
    udv=1000*sort(unique(dv[sel]))
    plot(udv,pp,xlab="Value Difference (ms juice)",ylab="Percent choose image",ylim=c(0,1),pch=19)
    segments(udv,pp-sd,udv,pp+sd,lwd=2)
    y=t(rbind(Nimg.comb,Ntot.comb-Nimg.comb))
    #try to fit; data may be bad
    fit=try(glm(y~udv,family=binomial(link="logit"))) 
    if (!(class(fit)[1]=="try-error")){
      beta=coef(fit)
      curve(logistic(beta[1]+beta[2]*x),from=min(udv),to=max(udv),add=T,lty=2)
      #title(paste("Session=",sess,"  Category=",cnames[catnum],"\n Width=",1/beta[2],
      #            "  Image=",beta[1]/beta[2]))
      cf=c(beta[1]/beta[2],1/beta[2])
      names(cf)=c("Image","Width")
      simdat=sim(fit,1000)
      bb=simdat@coef #mc coefficient draws
      ci=array(dim=c(2,2))
      ci[1,]=quantile(bb[,1]/bb[,2],c(.025,.975))
      ci[2,]=quantile(1/bb[,2],c(.025,.975))
      rownames(ci)=c("Image","Width")
      out=list(coef=cf,ci=ci,session=sess,category=catnum)
    }
    else {out=list(coef=rep(NA,2),ci=matrix(NA,2,2),session=sess,category=catnum)}
  }
  else {out=list(coef=rep(NA,2),ci=matrix(NA,2,2),session=sess,category=catnum)}
}

jags_sess_plot <- function(sess,catnum){
  #plot jags fitting result for a given session and category
  #of course, jags is bayesian, so we will take MAP estimates of parameters
  #THIS FUNCTION WILL ADD TO PLOT, NOT START A NEW PLOT WINDOW!
  
  sel=(piccat==catnum)&(session==sess)
  if (any(sel)){
    udv=1000*sort(unique(dv[sel]))
  }
  
  omegaind = grep(paste('choice.scale\\[',sess,'\\]',sep=''),names(xx)) #width parameter
  vvind = grep(paste('^v\\[',sess,',',catnum,sep=''),names(xx)) #pse/value
  
  #now grab median values of these (same as MAP estimate, since these are)
  omegasamp = xx[,omegaind]
  vvsamp = xx[,vvind]
  
  #use the density function to calculate map estimate
  omega = 1000*findmap(omegasamp,1024)
  vv = 1000*findmap(vvsamp,1024)
  
  #plot a curve
  curve(logistic((x+vv)/omega),from=min(udv),to=max(udv),add=T)
}



########################## let's make a plot #############################

#ready some params

svg(file='figure_sessfit.svg',width=18,height=10)

op = par()
par(mfrow = c(1,2),lwd=3,pch=19,font=1,ps=20,mar=c(5.1,7.1,4.1,2.1),cex.lab=1.25)

#plot good-looking data
sess = 39 
cat = 4
sessfit(sess,cat)
jags_sess_plot(sess,cat)

#plot crappy data
sess = 5  
cat = 4
sessfit(sess,cat)
jags_sess_plot(sess,cat)

#set graphics back right
par(op)
dev.off()

########################## plot session-by-session fits and original fits #############################

#ready some params
source("ppv_fitting.R")  #useful plotting functions
source("conventional_fits.R") #more conventional day-to-day fits of data
mnames=mvec

svg(file='figure_allsess_corrected.svg')
#plot session-by-session values for a given category
catnum=2
#plot pooled multi-level fits
plotsbys(catnum,xx,sessvec,datevec)

#plot standard fits alongside
#requires having fit.sets from explore_ppv_data.R already run
sc.mat=cbind(rep(1:numsess,each=numcat),rep(1:numcat,times=numsess)) #unique combinations of session and category
plotcat(catnum,fit.sets,which.fit=1,add.plot=T,jitt=0.25,colset=rep("gray",length(fit.sets)))

dev.off()
#ppv_fitting.R
#miscellaneous functions used to fit ppv choice curves
#called by explore_ppv_data.R

library("arm", quietly=T)
library("robust", quietly=T)

logistic <- function(x){1/(1+exp(-x))}

###### model log likelihood functions ######
LLrobit <- function(beta,x,n,N){
  #calculate log-likelihood of data for robit regression
  #x is the independent variable
  #beta=[mu, beta, dof] are the parameters of the t-distribution (t=beta*(x-mu))
  #n is the number of successes in N tries
  tt=beta[1]+beta[2]*x
  logp = pt(tt,beta[3],log.p=TRUE) #log(p)
  logq = pt(tt,beta[3],log.p=TRUE,lower.tail=FALSE) #log(1-p)
  LL = n*logp+(N-n)*logq + log(choose(N,n))
  out=sum(LL)
}

LLblogit <- function(beta,x,n,N){
  #calculate log likelihood for expanded logit model with biases
  #x is the independent variable
  #beta=[mu, beta, gamma, alpha] are the params, such that
  #lim{x->Inf} p(x) = alpha and lim{x->-Inf} p(x) = gamma
  #p(x) = {\gamma + \alpha e^{\beta x}} \over {1 + e^{\beta x}}
  # = gamma * (1 + exp(beta * x - log(alpha/gamma))) / (1 + exp(beta *x))
  #n is the number of successes in N tries
  zz = beta[1] + beta[2]*x
  logp = log(beta[3]) + log1p(exp(zz + log(beta[4]) - log(beta[3]))) - log1p(exp(zz))
  logq = log(1-beta[3]) + log1p(exp(zz + log(1-beta[4]) - log(1-beta[3]))) - log1p(exp(zz))
  LL = n*logp+(N-n)*logq + log(choose(N,n))
  out=sum(LL)
}

###### code to fit models by max log likelihood ######
fitLL <- function(data, LLfun, inits, modstr,...){
  # fit a model given its data (as retrieved by assemble.data),
  # log-likelihood function (LLfun), and initial guesses for params (inits)
  n = data[,1]
  N = rowSums(data)
  x = as.numeric(rownames(data)) #convert to ms for stability
  fitobj = optim(inits, function(b){-LLfun(b,x,n,N)}, control=list(maxit=10000),...)
  fit = list(coefficients = fitobj$par, LL = -fitobj$value, nobs = length(n), call=modstr,
             exit=fitobj$convergence)
  return(fit)
}

###### pull data corresponding to session and category ######
assemble.data <- function(sess, catnum){
  # marshal all data for a given session and image category into a form suitable for regression
  
  # which sessions to use
  sel=(piccat==catnum)&(session==sess)
  
  if (any(sel)){
    # accumulate total choices at each value of dv
    Ntot.tab=table(Ntot[sel],dv[sel])
    Ntot.vals=as.integer(rownames(Ntot.tab))
    Ntot.comb=Ntot.vals%*%Ntot.tab
    
    # accumulate total image choices at each value of dv
    Nimg.tab=table(Nimg[sel],dv[sel])
    Nimg.vals=as.integer(rownames(Nimg.tab))
    Nimg.comb=Nimg.vals%*%Nimg.tab
    
    y=t(rbind(Nimg.comb,Ntot.comb-Nimg.comb))
    return(y)
  }  
  else {return(NULL)}
}

###### fit session data ######
sessfit <- function(sess,catnum,model='logit'){
  # given session and category, fit model
  y = assemble.data(sess,catnum)
  
  if (!is.null(y)){
    udv = as.numeric(rownames(y))
    fit = switch(model,
                 logit = try(glm(y ~ udv, family=binomial(link='logit'))),
                 odlogit = try(glm(y ~ udv, family=quasibinomial(link='logit'))),
                 robust = try(glmRob(y ~ udv, family=binomial(link='logit'))),
                 robit = fitLL(y,LLrobit,c(0,100,3),'robit',method='L-BFGS-B',lower=c(-Inf,-Inf,1e-5)),
                 blogit = fitLL(y,LLblogit,c(0,100,0.1,0.9),'blogit'))
    
    if ((class(fit)[1]=="try-error")){
      return(NULL)
    }
    else {
      return(fit)
    }
  }
  else {
    return(NULL)
  }
}

###### extract params, confidence intervals, and model fit info ######
process.fit <- function(fit){
  # extract some useful quantities for plotting from a fit object
  modstr = grep('^\\w+',toString(fit$call[1]),value=T)
  beta=coef(fit)
  cf=c(beta[1]/beta[2],1/beta[2])
  names(cf)=c("Image","Width")
  if (('exit' %in% names(fit)) || any(is.na(beta)) || (fit$df.residual == 0) ){
    out = NULL
  }
  else {
    ci=array(dim=c(2,2))
    simdat = try(sim(fit,1000), silent=T)
    if ((class(simdat)[1]=="try-error")){ #give up
      ci[1,] = c(NA,NA)
      ci[2,] = c(NA,NA)
    }
    else { #calculate quantiles
      bb=simdat@coef #mc coefficient draws
      ci[1,]=quantile(bb[,1]/bb[,2],c(.025,.975))
      ci[2,]=quantile(1/bb[,2],c(.025,.975))
    }
    rownames(ci)=c("Image","Width")
    out=list(coef=cf,ci=ci)
  }
  
}

bic.fit <- function(fit){
  LLobj = try(logLik(fit), silent=T)
  if ((class(LLobj)[1]=="try-error")){
    # for fits I had to code myself, or for which 
    out = list(LL=fit$LL, N=fit$nobs, k=length(fit$coefficients))
  }
  else {
    out = list(LL=LLobj[1], N=nobs(LLobj), k=attributes(LLobj)$df)
  }
}

###### plotting functions ######
sessplot <- function(sess,catnum,model='logit'){
  # given a session, category number, and model fit, plot data
  y = assemble.data(sess,catnum)
  if (!is.null(y)){
    udv = as.numeric(rownames(y))
    Ntot = rowSums(y)
    pp = y[,1]/Ntot
    qq = 1-pp
    sd = sqrt(pp*qq/Ntot)
    plot(udv,pp,xlab="dv",ylab="Probability choose image")
    segments(udv,pp-sd,udv,pp+sd)
  }
  
  # now try to fit data
  fit = sessfit(sess,catnum,model)
  
  if (!is.null(fit)){
    beta=coef(fit)
    pf <- process.fit(fit)
    cf = pf$coef
    ci = pf$ci
    curve(logistic(beta[1]+beta[2]*x),from=min(udv),to=max(udv),add=T)
    title(paste("Session=",sess,"  Category=",cnames[catnum],"\n Width=",1/beta[2],
                "  Image=",beta[1]/beta[2]))
  }
}

plotcat <- function(catnum,fit.sets,which.fit,add.plot=FALSE,jitt=0,colset=rainbow(50)){
  imval=sapply(fit.sets[[which.fit]],
                         function(x){if (is.null(x)) NA else x$coef[1]},simplify=TRUE)
  imserr=sapply(fit.sets[[which.fit]],
                function(x){if (is.null(x)) c(NA,NA) else x$ci[1,]},simplify=TRUE)
  imserr=t(imserr)
  # convert to ms
  imval = 1000 * imval
  imserr = 1000 * imserr
  
  sel=(sc.mat[,2]==catnum)
  this.sess=sc.mat[sel,1]
  xrng=this.sess+jitt
  if(!add.plot){
    plot(xrng,imval[sel],xlab="Session",ylab="Image value (ms juice)",pch=16,
         ylim=c(-40,40),col=colset[which.fit])
  }
  else {
    points(xrng,imval[sel],xlab="Session",ylab="Image value (ms juice)",pch=16,
           ylim=c(-40,40),col=colset[which.fit])
  }
  segments(xrng,imserr[sel,1],this.sess,imserr[sel,2],col=colset[which.fit])
  abline(h=0)
  #title(paste("Image value for category=",cnames[catnum],sep=""))
  out=cbind(imval[sel],imserr[sel,])
}

plotquants <- function(qq,varlabs=c()){
  #plots quantiles of distributions in a vertical fashion a la Gelman and Hill
  #qq is a matrix with (possibly) named rows and columns corresponding to quantiles
  #if qq is the output of jags, these are the 2.5%, 25%, 50%, 75%, and 97.5% values
  namevec=rownames(qq) #get names of variables
  yvec=dim(qq)[1]:1 
  plot(qq[,3],yvec,pch=16,ylab="",xlab="",yaxt="n",xlim=range(qq),frame=F)
  if (is.null(varlabs)){
    axis(at=yvec,labels=namevec,2,las=2)
  } else {
    axis(at=yvec,labels=varlabs,2,las=2)
  }
  segments(x0=qq[,2],y0=yvec,x1=qq[,4],lwd=2)
  segments(x0=qq[,1],y0=yvec,x1=qq[,5],lwd=1)
}

plotsbys <- function(catnum,data,monkvec,datevec=rep(NA,length(monkvec)),add.plot=FALSE){
  #plot session-by-session values for each monk 
  #not necessarily sorted by date
  #catnum is the category number
  #data is a matrix of posterior samples, each column a variable
  #monkvec is the monk associated with each session
  #datevec is a date vector used or order sessions
  vnames=colnames(data)
  v.inds=grepl(paste("^v\\[\\d*,",catnum,sep=""),vnames)
  xsel=1000*data[,v.inds] #scale in ms
  quants=array(dim=c(dim(xsel)[2],3))
  for (ind in 1:dim(xsel)[2]){
    quants[ind,]=quantile(xsel[,ind],c(0.025,0.5,0.975))
  }
  umonk=unique(monkvec)
  colmap=rainbow(length(umonk))
  sstart=0
  if (!add.plot){
    plot.new()
    plot.window(xlim=c(0,length(monkvec)+1),ylim=range(quants))
  }
  for (ind in 1:length(umonk)){
    mm=umonk[ind]
    sel=(monkvec==mm)
    thisdates = datevec[sel]
    dperm = sort.int(thisdates,index.return=T) #permutation that puts date in date order
    sel = which(sel)[dperm$ix] #reorder to grab in correct date order
    xrng=sstart+(1:length(sel))
    points(xrng,quants[sel,2],pch=16,
           ylim=range(quants),col=colmap[ind])
    segments(x0=xrng,y0=quants[sel,1],y1=quants[sel,3],col=colmap[ind])
    text(x=mean(xrng),y=range(quants)[1],mnames[ind],col=colmap[ind])
    sstart=sstart+length(sel) #increase x offset by number of sessions this monk
  }
  axis(1)
  axis(2)
  title(main=paste('Session-by-session image values\n Category: ',cnames[catnum],sep=""),
        ylab='Image value (ms juice)',xlab='Session')
  abline(h=0)
}

plot.selections <- function(cats,data,sess){
  #plot select category session-by-session values for a given monk
  #cats is a vector of categories
  #data is a matrix of quantiles (quartiles plus 95% ci endpoints)
  #sess is a vector of sessions to use
  
  vnames=rownames(data)
  sel.inds=grepl("^v\\[",vnames) #get all session estimates
  xx=data[sel.inds,]
  xnames=rownames(xx)
  
  numsess=length(sess)
  sess.matches=regexpr("(?<=v\\[)\\d+",xnames,perl=T) #indices for session
  sess.inds=as.integer(regmatches(xnames,sess.matches))
  sess.sel=sess.inds %in% sess
  
  numcats=length(cats)
  cat.matches=regexpr("(?<=,)\\d+",xnames,perl=T) #indices for category
  cat.inds=as.integer(regmatches(xnames,cat.matches))
  ucats=unique(cats)
  
  jitter=0.1
  plot.new()
  plot.window(xlim=c(min(sess)-1,max(sess)+1),ylim=range(xx[sess.sel,]))
  colmap=c("gray","red","blue","green")
  
  for (ind in 1:numcats){
    this.sel=sess.sel&(cat.inds==ucats[ind])
    xrng=(ind-1)*jitter+sess.inds[this.sel]
    points(xrng,xx[this.sel,3],pch=16,col=colmap[ucats[ind]])
    segments(x0=xrng,y0=xx[this.sel,1],y1=xx[this.sel,5],col=colmap[ucats[ind]])
  }
  abline(h=0)
  axis(1)
  axis(2)
  title(main="Session by Session category values",xlab="Session",ylab="Image Value (ms juice)")
}
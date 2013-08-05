#load up data
source("all_ppv_data.R")
source("conventional_fits.R") #more conventional day-to-day fits of data

#which monk is which
mvec=c('E','Os','Ot','D','S','C','B','N')

#which category
cnames=c("Neutral","Female","Dominant Male","Subordinate Male")

#load up functions for plotting
source('ppv_fitting.R')

#plot curve widths
catnum = 2
plotcat(catnum,fit.sets,which.fit=1,sessvec,datevec,
        add.plot=F,which.coeff=2,colset=rainbow(length(mvec)),add.title=F)
title(main=paste('Session-by-session juice precisions\n Category: ',cnames[catnum],sep=""),
      ylab='Juice precision (ms juice)',xlab='Session')


######## Bayesian version #################
load('ppv_results')

doplot <- function(data,monkvec,datevec=1:length(monkvec),
                     add.plot=FALSE){
  #plot session-by-session values for each monk 
  #catnum is the category number
  #data is a matrix of posterior samples, each column a variable
  #monkvec is the monk associated with each session
  #datevec is a date vector used to order sessions
  vnames=colnames(data)
  v.inds=grepl("^choice.scale\\[\\d*",vnames)
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
    plot.window(xlim=c(0,length(monkvec)+1),ylim=c(0,100))
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
    text(x=mean(xrng),y=-110, mvec[ind],col=colmap[ind])
    sstart=sstart+length(sel) #increase x offset by number of sessions this monk
  }
  axis(1)
  axis(2)
  title(main=paste('Session-by-session image values\n Category: ',cnames[catnum],sep=""),
        ylab='Image value (ms juice)',xlab='Session')
  abline(h=0)
}

doplot(xx,sessvec,datevec,add.plot=FALSE)
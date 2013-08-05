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
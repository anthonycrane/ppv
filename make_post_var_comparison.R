#make_post_var_comparison.R
#compare day-to-day variances across monks

#based on plots_ppv.R

#IMPORTANT: modeling note: day-to-day variance is modeled as the same for all categories with
#each monk, so only one posterior per monk in this case

#load up some data
load("ppv_results") #load data output
source("all_ppv_data.R") #raw data

#some useful reference variables
vnames=colnames(xx)
pnames=colnames(pp)
cnames=c("Neutral","Female","Dominant Male","Subordinate Male")
mnames=c("E","Os","Ot","D","S","C","B","N")

###############################################
#get ready to plot!

#get quantiles for each monk
whichquants = c(0.025,0.5,0.975)
v.inds=grepl("^sess.std\\[\\d*",vnames)
xsel=1000*xx[,v.inds] #scale in ms
quants=array(dim=c(dim(xsel)[2],3))
for (ind in 1:dim(xsel)[2]){
  quants[ind,]=quantile(xsel[,ind],whichquants)
}
rownames(quants) = vnames[v.inds]
colnames(quants) = whichquants

svg(file='figure_post_var_comparison.svg')
par(lwd=3)
plot.new()
mjit = 0  #amount to space monks by horizontally
plot.window(xlim=c(0,nummonk+1),ylim = range(quants))

points(quants[,2],pch=16,col='black')
segments(x0=1:nummonk,y0=quants[,1],y1=quants[,3],col='black')

#add axes, title, etc.
xlabtick = 1:nummonk
axis(1,at=xlabtick,labels=mnames)
axis(2)
title(main='Image Value: Posterior Standard Deviation and Credible Intervals',
      ylab='Image value standard deviation (ms juice)',xlab='Subject')

dev.off()
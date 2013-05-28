#make_trend_allmonks.R
#make a plot of session-by-session image values for a given category
#based on plots_ppv_trend.R

load("ppv_trend_results") #load data output
source("all_ppv_data.R") #raw data
source('all_ppv_data_frame.R')
source("ppv_fitting.R")  #useful plotting functions

#set up some useful variables
vnames=colnames(xx)
pnames=colnames(pp)
cnames=c("Gray Square","Female Perinea","Dominant Male","Subordinate Male")
mnames=c("E","Os","Ot","D","S","C","B","N")

#set up device
svg(file='figure_trend_allmonks.svg')

#plot session-by-session values for a given category
catnum=2
#plot pooled multi-level fits
plotsbys(catnum,xx,sessvec,datevec)

#close device
dev.off()

#now plot posteriors on trendlines

svg(file='figure_post_trend_comparison.svg')
par(lwd=3)

#get quantiles for each monk
whichquants = c(0.025,0.5,0.975)
v.inds=grepl("vslope",vnames)
xsel=1000*xx[,v.inds] #scale in ms
quants=array(dim=c(dim(xsel)[2],3))
for (ind in 1:dim(xsel)[2]){
  quants[ind,]=quantile(xsel[,ind],whichquants)
}
rownames(quants) = vnames[v.inds]
colnames(quants) = whichquants

#now rescale quantiles
quants = quants/sdvec #remember: we z-scored session number before fitting

colmap = c('gray','red','blue','green')

plot.new()
mjit = 2 #amount to space monks by horizontally
plot.window(xlim=c(0,(numcat+mjit)*nummonk+1),ylim = c(-50,50))

for (ind in 1:numcat){
  whichrows = ((ind-1)*nummonk + 1):(ind*nummonk)
  thisx = ind + (numcat+mjit)*(0:(nummonk-1))
  points(thisx,quants[whichrows,2],pch=16,col=colmap[ind])
  segments(x0=thisx,y0=quants[whichrows,1],y1=quants[whichrows,3],col=colmap[ind])
}

#add axes, title, etc.
xlabtick = (numcat+mjit)*(0:(nummonk-1)) + 2.5
axis(1,at=xlabtick,labels=mnames)
axis(2)
title(main='Image Value: Posterior Time Trends and Credible Intervals',
      ylab='Rate of change (ms juice/session)',xlab='Subject')
abline(h=0)
dev.off()
#make_post_mean_comparison
#dot plot of posterior mean distributions for each monk and each category
#based on code in plots_ppv.R

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
v.inds=grepl("^V\\[\\d*,",vnames)
xsel=1000*xx[,v.inds] #scale in ms
quants=array(dim=c(dim(xsel)[2],3))
for (ind in 1:dim(xsel)[2]){
  quants[ind,]=quantile(xsel[,ind],whichquants)
}
rownames(quants) = vnames[v.inds]
colnames(quants) = whichquants

colmap = c('gray','red','blue','green')

svg(file='figure_post_mean_comparison.svg')
par(lwd=3)
plot.new()
mjit = 2 #amount to space monks by horizontally
plot.window(xlim=c(0,(numcat+mjit)*nummonk+1),ylim = range(quants),lwd=3)

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
title(main='Image Value: Posterior Means and Credible Intervals',
      ylab='Image value (ms juice)',xlab='Subject')
abline(h=0)
dev.off()
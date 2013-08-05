#make_figure_sessfit.R
#construct a figure illustrating:
#1) well-behaved behavioral data with a standard logistic fit
#2) ill-behaved data with an ill-behaved logistic fit
#3) partial pooling fits to both of these
#based on explore_ppv_data.R and ppv_fitting.R

#load up data
source("all_ppv_data.R")
load('ppv_results')

#which monk is which
mvec=c('E','Os','Ot','D','S','C','B','N')

#which category
cnames=c("Neutral","Female","Dominant Male","Subordinate Male")

#load up functions for plotting
source('ppv_fitting.R')


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

#plot messy data
sess = 5  
cat = 4
sessfit(sess,cat)
jags_sess_plot(sess,cat)

#set graphics back right
par(op)
dev.off()

########################## plot session-by-session fits and original fits #############################

#ready some params
source("conventional_fits.R") #more conventional day-to-day fits of data
mnames=mvec

svg(file='figure_allsess_corrected.svg')
#plot session-by-session values for a given category
catnum=2

#plot standard fits
#requires having fit.sets from explore_ppv_data.R already run
plotcat(catnum,fit.sets,which.fit=1,sessvec,datevec,
        add.plot=F,jitt=0.25,colset=rep("gray",length(fit.sets)))

#plot pooled multi-level fits on top
plotsbys(catnum,xx,sessvec,datevec,add.plot=T)

dev.off()

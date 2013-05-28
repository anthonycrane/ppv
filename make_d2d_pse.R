#make_d2d_pse.R

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
mnum = 2;

svg(file='figure_d2d_pse.svg')
vtemp=list(c(),c(),c(),c())
for (catnum in 1:4){
  #v.inds=grepl(paste("^v\\[\\d*,",catnum,sep=""),vnames) #use all monks
  usable = which(mnum==sessvec) #sessions belonging to this monk
  #the following regular expression takes all names of the form
  #v[<any session for this monk>,<category>
  rexpr = paste("^v\\[(",paste(usable,sep="",collapse="|"),"),",catnum,sep="")
  v.inds=grepl(rexpr,vnames)
  xsel=1000*xx[,v.inds]
  vtemp[[catnum]]=as.vector(as.matrix(xsel))
}
xsel=as.data.frame(vtemp)
xsel = xsel[,]
colnames(xsel)=c("gray","peri","dom","sub")
numpts = 1000
pairs(~gray+peri+dom+sub,data=xsel,
      subset=sample(dim(xsel)[1],numpts,replace=F),
      labels = cnames) #take only a small subsample
cor(xsel,method="spearman")
dev.off()
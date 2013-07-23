source("all_ppv_data.R")
source("ppv_fitting.R")

#fit all sessions
fit.sets = list()
fit.objs = list()
sc.mat=cbind(rep(1:numsess,each=numcat),rep(1:numcat,times=numsess)) #unique combinations of session and category
modlist = c('logit','odlogit')
for (ind in 1:numsess){
  for (ind2 in 1:numcat){
    for (model in modlist){
      fit = sessfit(ind,ind2,model)
      if (!is.null(fit)){
        pf = process.fit(fit)
      }
      else {
        pf = NULL
      }
      fit.objs[[model]][[(ind-1)*numcat+ind2]] = fit
      fit.sets[[model]][[(ind-1)*numcat+ind2]] = pf
    }
  }
}
dump(c("fit.sets","fit.objs","modlist"),file="conventional_fits.R")

for (mod in modlist){
  LLobjs = lapply(fit.objs[[mod]], function(x){if(is.null(x)) NA else logLik(x)})
  LL = sum(sapply(LLobjs, function(x){x[1]}, simplify=T), na.rm=T)
  N = sum(sapply(LLobjs, function(x){if(is.na(x)) NA else nobs(x)}, simplify=T), na.rm=T)
  k = sum(sapply(LLobjs, function(x){if(is.na(x)) NA else attributes(x)$df}, simplify=T), na.rm=T)
}
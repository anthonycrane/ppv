source("all_ppv_data.R")
source("ppv_fitting.R")

#fit all sessions
fit.sets = list()
fit.objs = list()
sc.mat=cbind(rep(1:numsess,each=numcat),rep(1:numcat,times=numsess)) #unique combinations of session and category
modlist = c('logit','odlogit','robit','blogit')
for (ind in 1:numsess){
  for (ind2 in 1:numcat){
    print(c(ind,ind2))
    for (model in modlist){
      fit = sessfit(ind,ind2,model)
      if (!is.null(fit)){
        pf = process.fit(fit)
        bicvars = bic.fit(fit)
      }
      else {
        pf = NULL
        bicvars = NULL
      }
      fit.sets[[model]][[(ind-1)*numcat+ind2]] = c(pf,bicvars)
    }
  }
}
dump(c("fit.sets","modlist"),file="conventional_fits.R")
rm(list=ls())
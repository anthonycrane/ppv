source("all_ppv_data.R")
source("ppv_fitting.R")

#allocate variables
fit.sets = list()
modlist = c('logit','odlogit','robit','blogit')

#fit all sessions
for (model in modlist){
  this_model = list()
  for (ind in 1:numsess){
    this_sess = list()
    for (ind2 in 1:numcat){
      print(c(model,ind,ind2))
      fit = sessfit(ind,ind2,model)
      if (!is.null(fit)){
        pf = process.fit(fit)
        bicvars = bic.fit(fit)
      }
      else {
        pf = NULL
        bicvars = NULL
      }
      this_sess[[ind2]] = c(pf,bicvars)
    }
    this_model[[ind]] = this_sess
  }
  fit.sets[[model]] = this_model
}
dump(c("fit.sets","modlist"),file="conventional_fits.R")
rm(list=ls())
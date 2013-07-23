# compare conventional fits with hierarchical models based on bic

source('conventional_fits.R')

LL = c()
N = c()
k = c()
for (ind in modlist){
  LL[ind] = sum(sapply(fit.sets[[ind]],function(x){if(is.null(x)) NA else x$LL}, simplify=T), na.rm=T)
  N[ind] = sum(sapply(fit.sets[[ind]],function(x){if(is.null(x)) NA else x$N}, simplify=T), na.rm=T)
  k[ind] = sum(sapply(fit.sets[[ind]],function(x){if(is.null(x)) NA else x$k}, simplify=T), na.rm=T)
}

bic = -2*LL + k * log(N)
aic = -2*LL + 2*k
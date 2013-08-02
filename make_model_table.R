#get model fitting metrics for the model comparison table

#need to load this library so dic objects will display correctly
library(rjags)

load('dic_results')
#dic is a list (in order of the table in the ms); element names correspond to .bug files
#each entry is a dic object listing Dbar, pD, and DIC


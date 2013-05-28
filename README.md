# Background:
The following code can be used to replicate the results of [this paper](http://insert.link.here). It fits a Bayesian hierarchical model to choice data from a series of experiments in which monkeys chose between juice rewards and juice rewards paired with social images.

# Data:
The data are contained in a csv file (`all_ppv_data.csv`) and as individual variables in R's "dump" data format (`all_ppv_data.R`). The latter is read by the modeling code to run inference.

# Dependencies:
The code is written in R and makes use of [JAGS](http://mcmc-jags.sourceforge.net) and the `rjags` package to run Markov Chain Monte Carlo. Very limited use is made of the `arm` and `MASS` packages for basic model fitting in `explore_ppv_data.R`.

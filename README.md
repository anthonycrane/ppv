# Background:
The following code can be used to replicate the results of [this paper](http://insert.link.here). It fits a Bayesian hierarchical model to choice data from a series of experiments in which monkeys chose between juice rewards and juice rewards paired with social images.

# Data:
The data are contained in a csv file (`all_ppv_data.csv`) and as individual variables in R's "dump" data format (`all_ppv_data.R`). The latter is read by the modeling code to run inference.

# Dependencies:
The code is written in R and makes use of [JAGS](http://mcmc-jags.sourceforge.net) and the `rjags` package to run Markov Chain Monte Carlo. Very limited use is made of the `arm` and `MASS` packages for basic model fitting in `explore_ppv_data.R`.

# Workflow:
## Running Analyses:
* `explore_ppv_data.R` uses conventional logit fits to model image values separately for each session. This is (more or less) the approach used in previous papers and is included for comparison. This file compiles parameters for all the fits and saves them in dump format in `conventional_fits.R` for later use.
* `runjags_ppv.R` performs MCMC inference on the model(s) and saves the samples drawn from the Markov chains to files for later analysis. __more detail later about how to select model and output file__
* __which bugs files to use__

## Reproducing Tables and Figures:
* Data in Table 1 in the paper, summarizing explained variance and pooling, can be retrieved by running `make_fitting_table.R` on the output of the above analyses.
* Figure 3, with comparisons between hierarchical and unpooled estimates, can be reproduced with the code in `make_figure_sessfit.R`.
* Figure 4, comparing fitted and simulated data, can be generated with `dat_vs_sim.R`.
* A similar comparison of overdispersion in real and fitted data can be produced via `overdisp_fit_vs_sim.R`.
* Figure 5, which displays comparisons of the posterior mean and variance for image values across subjects, can be reproduced using `make_post_mean_comparison.R` and `make_post_var_comparison.R`.
* Figure 6, the model including time trend, was made using `make_trend_allmonks.R` and `make_d2d_pse.R`.

## Other files:
* `ppv_fitting.R` contains helper functions for fitting and plotting data. 
 

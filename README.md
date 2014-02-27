# Background:
The following code can be used to replicate the results of [this paper](http://journal.frontiersin.org/Journal/10.3389/fnins.2013.00165/abstract). It fits a Bayesian hierarchical model to choice data from a series of experiments in which monkeys chose between juice rewards and juice rewards paired with social images.

# Data:
The data are contained in a csv file (`all_ppv_data.csv`) and as individual variables in R's "dump" data format (`all_ppv_data.R`). The latter is read by the modeling code to run inference.

# Dependencies:
The code is written in R and makes use of [JAGS](http://mcmc-jags.sourceforge.net) and the `rjags` package to run Markov Chain Monte Carlo. Very limited use is made of the `arm`, `robust`, and `MASS` packages for basic model fitting in `explore_ppv_data.R`.

# Workflow:
## Running Analyses:
* `make_conventional_fits.R` fits a variety of conventional models (logit, overdispersed logit, robust logit, logit with biases) and stores the results in `conventional_fits.R` in R's `dump` format.
* `explore_ppv_data.R` uses conventional logit fits to model image values separately for each session. Runs some basic plots. Relies on `conventional_fits.R` as produced by `make_conventional_fits.R`. 
* `runjags_ppv.R` performs MCMC inference on the model(s) and saves the samples drawn from the Markov chains to files for later analysis. Setting the `modstr` variable selects which model to fit for analysis. Correspondences between model strings and model numbers in the text are listed in the file and restated below. 
* Models are specified in BUGS/JAGS language in the models listed below. Models are (so far) only tested in JAGS. The correct model will be used in `runjags_ppv.R` so long as the `modstr` variable is specified correctly. Results will be saved in separate files depending on `modstr`. 
    - Model 0: `nopooling.bug`
    - Model 1: `monkcatpooling.bug`
    - Model 2: `monkpooling.bug`
    - Model 3: `nomonknocat.bug`
    - Model 4: `nocat.bug`
    - Model 5: `nomonk.bug`
    - Model 6: `model.bug`
    - Model 7: `trend.bug`

## Reproducing Tables and Figures:
* Data in Table 1 in the paper, summarizing explained variance and pooling, can be retrieved by running `make_fitting_table.R` on the output of the above analyses.
* Data in Table 2 are based on running `runjags_ppv.R` for each model and combining outputs. However, they can be directly reproduced by running `make_model_table.R`, which loads the DIC objects contained in the binary file `dic_results`.
* Figure 3, comparing fitted and simulated data, can be generated with `dat_vs_sim.R`.
* A similar comparison of overdispersion in real and synthetic data can be produced via `overdisp_fit_vs_sim.R`.
* Figure 4, with comparisons between hierarchical and unpooled estimates, can be reproduced with the code in `make_figure_sessfit.R`.
* The file `explore_curve_widths.R` plots both conventional and Bayesian hierarchical fits of juice precisions/curve widths for each session as a function of time.
* Figure 5, which displays comparisons of the posterior mean and variance for image values across subjects, can be reproduced using `make_post_mean_comparison.R` and `make_post_var_comparison.R`.
* Figure 6, the model including time trend, was made using `make_trend_allmonks.R` and `make_d2d_pse.R`.

## Other files:
* `ppv_fitting.R` contains helper functions for fitting and plotting data. 
 

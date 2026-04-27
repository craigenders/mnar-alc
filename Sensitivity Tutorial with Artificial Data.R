#------------------------------------------------------------------------------#
# LOAD R PACKAGES ----
#------------------------------------------------------------------------------#

options (scipen = 999)

# load packages
library(ggplot2)
library(patchwork)
library(rblimp)
set_blimp('/applications/blimp/blimp-nightly')

#------------------------------------------------------------------------------#
# READ DATA ----
#------------------------------------------------------------------------------#

# github url for raw data
filepath <- 'https://raw.githubusercontent.com/craigenders/mnar-alc/main/mnar.csv'

# create data frame from github data
mnar <- read.csv(filepath, stringsAsFactors = T)

# plotting functions
source('https://raw.githubusercontent.com/blimp-stats/blimp-book/main/misc/functions.R')
source('https://raw.githubusercontent.com/craigenders/mnar-mlm/main/mnar-plotting.R')

#------------------------------------------------------------------------------#
# CMAR MODELS ----
#------------------------------------------------------------------------------#

# MODEL 1: CMAR ----
model1_cmar <- rblimp(
  data = mnar,
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable
  latent = 'id = b0i b1i', # define random effects as latent variables
  fixed = 'month med', # complete predictors
  model = '
    # level-2 equations
    b0i ~ intercept@b0 med@b2;
    b1i ~ intercept@b1 med@b3;
    b0i b1i ~~ b0i b1i;
    # level-1 equation
    dpdd ~ intercept@b0i month@b1i;',
  parameters = ' 
    # compute time-specific effect sizes
    d_month0 = ( (b0+b2) - (b0) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month1 = ( (b0+b2 + 1*(b1+b3)) - (b0 + 1*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month2 = ( (b0+b2 + 2*(b1+b3)) - (b0 + 2*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month3 = ( (b0+b2 + 3*(b1+b3)) - (b0 + 3*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month4 = ( (b0+b2 + 4*(b1+b3)) - (b0 + 4*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);',
  seed = 90291, # random number seed
  chains = 4, # number of mcmc processes
  burn = 20000, # burn-in iterations
  iter = 20000, # iterations for analysis summaries
  nimps = 20 # save imputations for graphing trajectories
)

# print output
output(model1_cmar)

# plot model-predicted means
bivariate_plot(dpdd.predicted ~ month | med, 
               model = model1_cmar, 
               discrete_x = 'month', 
               points = F, ci = F) + ylim(0,15)

# MODEL 2: CMAR With Auxiliary Variable ----
model2_aux <- rblimp(
  data = mnar,
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable
  latent = 'id = b0i b1i', # define random effects as latent variables
  fixed = 'month med', # complete predictors
  model = '
    # level-2 equations
    b0i ~ intercept@b0 med@b2;
    b1i ~ intercept@b1 med@b3;
    b0i b1i ~~ b0i b1i;
    # level-1 equation
    dpdd ~ intercept@b0i month@b1i;
    # auxiliary variable equation
    ftnd ~ b0i b1i;',
  parameters = ' 
    # compute time-specific effect sizes
    d_month0 = ( (b0+b2) - (b0) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month1 = ( (b0+b2 + 1*(b1+b3)) - (b0 + 1*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month2 = ( (b0+b2 + 2*(b1+b3)) - (b0 + 2*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month3 = ( (b0+b2 + 3*(b1+b3)) - (b0 + 3*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month4 = ( (b0+b2 + 4*(b1+b3)) - (b0 + 4*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);',
  seed = 90291, # random number seed
  chains = 4, # number of mcmc processes
  burn = 20000, # burn-in iterations
  iter = 20000, # iterations for analysis summaries
  nimps = 20 # save imputations for graphing trajectories
)

# print output
output(model2_aux)

# plot model-predicted means
bivariate_plot(dpdd.predicted ~ month | med, 
               model = model2_aux, 
               discrete_x = 'month', 
               points = F, ci = F) + ylim(0,15)

#------------------------------------------------------------------------------#
# TIME-RELATED DROPOUT ----
#------------------------------------------------------------------------------#

# DROPOUT PROBS: Linear Time ----
time_linear <- rblimp(
  data = mnar,
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable 
  timeid = 'month', # time identifier for dropout and lagged effects 
  dropout = 'dropout = dpdd (monotone)', # create dropout indicator
  model = 'dropout ~ intercept@-3 month month*med | intercept@0',
  seed = 90291, # random number seed
  chains = 4, # number of mcmc processes
  burn = 20000, # burn-in iterations
  iter = 20000, # iterations for analysis summaries
  nimps = 20 # save imputations for graphing trajectories
)

# print output
output(time_linear)

# DROPOUT PROBS: Quadratic Time ----
time_quadratic <- rblimp(
  data = mnar,
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable 
  timeid = 'month', # time identifier for dropout and lagged effects 
  dropout = 'dropout = dpdd (monotone)', # create dropout indicator
  model = 'dropout ~ intercept@-3 month month^2 month*med month^2*med | intercept@0',
  seed = 90291, # random number seed
  chains = 4, # number of mcmc processes
  burn = 20000, # burn-in iterations
  iter = 20000, # iterations for analysis summaries
  nimps = 20 # save imputations for graphing trajectories
)

# print output
output(time_quadratic)


# DROPOUT PROBS: Dummy-Coded Time ----
time_dummy <- rblimp(
  data = mnar,
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable 
  timeid = 'month', # time identifier for dropout and lagged effects 
  dropout = 'dropout = dpdd (monotone)', # create dropout indicator
  # (month == 1) (month == 2) etc. are dummy codes for time
  model = 'dropout ~ intercept@-3 (month==1) (month==2) (month==3) (month==4)
      (month==1)*med (month==2)*med (month==3)*med (month==4)*med | intercept@0',
  seed = 90291, # random number seed
  chains = 4, # number of mcmc processes
  burn = 20000, # burn-in iterations
  iter = 20000, # iterations for analysis summaries
  nimps = 20 # save imputations for graphing trajectories
)

# print output
output(time_dummy)

#------------------------------------------------------------------------------#
# PLOT MISSINGNESS PROBABILITIES ----
#------------------------------------------------------------------------------#

ymin <- 0
ymax <- .5

# plot observed probabilities
plot_pmiss_obs <- plot_means(dropout ~ month | med, 
    model = time_dummy,
    ylab = "Probability",
    title = "Observed Probabilities",
    group_labels = c("0" = "0", "1" = "1")) + ylim(ymin,ymax) +
    theme(legend.position = "top",legend.justification = "center") +
    scale_linetype_manual(values = c("dashed", "solid")) +
    geom_line(linewidth = .25)

# plot predicted probabilities from linear model
plot_pmiss_linear <- plot_means(dropout.1.probability ~ month | med, 
    model = time_linear,
    ylab = "Probability",
    title = "Linear Time",
    group_labels = c("0" = "0", "1" = "1")) + ylim(ymin,ymax) +
    theme(legend.position = "top",legend.justification = "center") +
    scale_linetype_manual(values = c("dashed", "solid")) +
    geom_line(linewidth = .25)

# plot predicted probabilities from quadratic model
plot_pmiss_quadratic <- plot_means(dropout.1.probability ~ month | med, 
    model = time_quadratic,
    ylab = "Probability",
    title = "Quadratic Time",
    group_labels = c("0" = "0", "1" = "1")) + ylim(ymin,ymax) +
    theme(legend.position = "top",legend.justification = "center") +
    scale_linetype_manual(values = c("dashed", "solid")) +
    geom_line(linewidth = .25)

# plot predicted probabilities from dummy-coded model
plot_pmiss_dummy <- plot_means(dropout.1.probability ~ month | med, 
    model = time_dummy,
    ylab = "Probability",
    title = "Dummy-Coded Time",
    group_labels = c("0" = "0", "1" = "1")) + ylim(ymin,ymax) +
    theme(legend.position = "top",legend.justification = "center") +
    scale_linetype_manual(values = c("dashed", "solid")) +
    geom_line(linewidth = .25)

# combine plots with patchwork
plot_pmiss <- (plot_pmiss_obs | plot_pmiss_linear) / (plot_pmiss_quadratic | plot_pmiss_dummy)
plot_pmiss

#------------------------------------------------------------------------------#
# SHARED PARAMETER MODELS ----
#------------------------------------------------------------------------------#

# MODEL 3: Shared Parameter Model With Random Effects Predicting Dropout ----
model3_sharedparam <- rblimp(
  data = mnar,
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable
  timeid = 'month', # time identifier for dropout and lagged effects 
  dropout = 'dropout = dpdd (monotone)', # create dropout indicator
  latent = 'id = b0i b1i', # define random effects as latent variables
  fixed = 'month med', # complete predictors
  model = '
    # level-2 equations
    b0i ~ intercept@b0 med@b2;
    b1i ~ intercept@b1 med@b3;
    b0i b1i ~~ b0i b1i;
    # level-1 equation
    dpdd ~ intercept@b0i month@b1i;
    # dropout model
    preknot = ifelse(month == 0, 0, 1); # dummy code to activate dropout predictors post-baseline
    dropout ~ intercept@-3 (month==1) (month==2) (month==3) (month==4)
      (month==1)*med (month==2)*med (month==3)*med (month==4)*med
      b0i*preknot b1i*preknot | intercept@0;',
  parameters = ' 
    # compute time-specific effect sizes
    d_month0 = ( (b0+b2) - (b0) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month1 = ( (b0+b2 + 1*(b1+b3)) - (b0 + 1*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month2 = ( (b0+b2 + 2*(b1+b3)) - (b0 + 2*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month3 = ( (b0+b2 + 3*(b1+b3)) - (b0 + 3*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month4 = ( (b0+b2 + 4*(b1+b3)) - (b0 + 4*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);',
  seed = 90291, # random number seed
  chains = 4, # number of mcmc processes
  burn = 20000, # burn-in iterations
  iter = 20000, # iterations for analysis summaries
  nimps = 20 # save imputations for graphing trajectories
)

# print output
output(model3_sharedparam)

# plot model-predicted means
bivariate_plot(dpdd.predicted ~ month | med, 
               model = model3_sharedparam, 
               discrete_x = 'month', 
               points = F, ci = F) + ylim(0,15)

# MODEL 5: Shared Parameter Model With Random Effects and Their Squares Predicting Dropout ----
model5_sharedparam <- rblimp(
  data = mnar,
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable
  timeid = 'month', # time identifier for dropout and lagged effects 
  dropout = 'dropout = dpdd (monotone)', # create dropout indicator
  latent = 'id = b0i b1i', # define random effects as latent variables
  fixed = 'month med', # complete predictors
  model = '
    # level-2 equations
    b0i ~ intercept@b0 med@b2;
    b1i ~ intercept@b1 med@b3;
    b0i b1i ~~ b0i b1i;
    # level-1 equation
    dpdd ~ intercept@b0i month@b1i;
    # dropout model
    preknot = ifelse(month == 0, 0, 1); # dummy code to activate dropout predictors post-baseline
    dropout ~ intercept@-3 (month==1) (month==2) (month==3) (month==4)
      (month==1)*med (month==2)*med (month==3)*med (month==4)*med
      b0i*preknot b1i*preknot (b0i^2)*preknot (b1i^2)*preknot | intercept@0;',
  parameters = ' 
    # compute time-specific effect sizes
    d_month0 = ( (b0+b2) - (b0) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month1 = ( (b0+b2 + 1*(b1+b3)) - (b0 + 1*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month2 = ( (b0+b2 + 2*(b1+b3)) - (b0 + 2*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month3 = ( (b0+b2 + 3*(b1+b3)) - (b0 + 3*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month4 = ( (b0+b2 + 4*(b1+b3)) - (b0 + 4*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);',
  seed = 90291, # random number seed
  chains = 4, # number of mcmc processes
  chains = 4, # number of mcmc processes
  burn = 60000, # burn-in iterations
  iter = 60000, # iterations for analysis summaries
  nimps = 20 # save imputations for graphing trajectories
)

# print output
output(model5_sharedparam)

# plot model-predicted means
bivariate_plot(dpdd.predicted ~ month | med, 
               model = model5_sharedparam, 
               discrete_x = 'month', 
               points = F, ci = F) + ylim(0,15)

#------------------------------------------------------------------------------#
# SELECTION MODELS MODELS ----
#------------------------------------------------------------------------------#

# MODEL 6: Selection Model With Outcome Predicting Dropout ----
model6_selection <- rblimp(
  data = mnar,
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable
  timeid = 'month', # time identifier for dropout and lagged effects 
  dropout = 'dropout = dpdd (monotone)', # create dropout indicator
  latent = 'id = b0i b1i', # define random effects as latent variables
  fixed = 'month med', # complete predictors
  model = '
    # level-2 equations
    b0i ~ intercept@b0 med@b2;
    b1i ~ intercept@b1 med@b3;
    b0i b1i ~~ b0i b1i;
    # level-1 equation
    dpdd ~ intercept@b0i month@b1i;
    # dropout model
    preknot = ifelse(month == 0, 0, 1); # dummy code to activate dropout predictors post-baseline
    dropout ~ intercept@-3 (month==1) (month==2) (month==3) (month==4)
      (month==1)*med (month==2)*med (month==3)*med (month==4)*med
      dpdd*preknot | intercept@0;',
  parameters = ' 
    # compute time-specific effect sizes
    d_month0 = ( (b0+b2) - (b0) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month1 = ( (b0+b2 + 1*(b1+b3)) - (b0 + 1*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month2 = ( (b0+b2 + 2*(b1+b3)) - (b0 + 2*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month3 = ( (b0+b2 + 3*(b1+b3)) - (b0 + 3*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month4 = ( (b0+b2 + 4*(b1+b3)) - (b0 + 4*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);',
  seed = 90291, # random number seed
  chains = 4, # number of mcmc processes
  burn = 20000, # burn-in iterations
  iter = 20000, # iterations for analysis summaries
  nimps = 20 # save imputations for graphing trajectories
)

# print output
output(model6_selection)

# plot model-predicted means
bivariate_plot(dpdd.predicted ~ month | med, 
               model = model6_selection, 
               discrete_x = 'month', 
               points = F, ci = F) + ylim(0,15)

# MODEL 7: Selection Model With Outcome and its Square Predicting Dropout ----
model7_selection <- rblimp(
  data = mnar,
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable
  timeid = 'month', # time identifier for dropout and lagged effects 
  dropout = 'dropout = dpdd (monotone)', # create dropout indicator
  latent = 'id = b0i b1i', # define random effects as latent variables
  fixed = 'month med', # complete predictors
  model = '
    # level-2 equations
    b0i ~ intercept@b0 med@b2;
    b1i ~ intercept@b1 med@b3;
    b0i b1i ~~ b0i b1i;
    # level-1 equation
    dpdd ~ intercept@b0i month@b1i;
    # dropout model
    preknot = ifelse(month == 0, 0, 1); # dummy code to activate dropout predictors post-baseline
    dropout ~ intercept@-3 (month==1) (month==2) (month==3) (month==4)
      (month==1)*med (month==2)*med (month==3)*med (month==4)*med
      dpdd*preknot (dpdd^2)*preknot | intercept@0;',
  parameters = ' 
    # compute time-specific effect sizes
    d_month0 = ( (b0+b2) - (b0) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month1 = ( (b0+b2 + 1*(b1+b3)) - (b0 + 1*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month2 = ( (b0+b2 + 2*(b1+b3)) - (b0 + 2*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month3 = ( (b0+b2 + 3*(b1+b3)) - (b0 + 3*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month4 = ( (b0+b2 + 4*(b1+b3)) - (b0 + 4*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);',
  seed = 90291, # random number seed
  chains = 4, # number of mcmc processes
  burn = 20000, # burn-in iterations
  iter = 20000, # iterations for analysis summaries
  nimps = 20 # save imputations for graphing trajectories
)

# print output
output(model7_selection)

# plot model-predicted means
bivariate_plot(dpdd.predicted ~ month | med, 
               model = model7_selection, 
               discrete_x = 'month', 
               points = F, ci = F) + ylim(0,15)

# MODEL 8: Selection Model With Outcome and its Lag Predicting Dropout ----
model8_selection <- rblimp(
  data = mnar,
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable
  timeid = 'month', # time identifier for dropout and lagged effects 
  dropout = 'dropout = dpdd (monotone)', # create dropout indicator
  latent = 'id = b0i b1i', # define random effects as latent variables
  fixed = 'month med', # complete predictors
  model = '
    # level-2 equations
    b0i ~ intercept@b0 med@b2;
    b1i ~ intercept@b1 med@b3;
    b0i b1i ~~ b0i b1i;
    # level-1 equation
    dpdd ~ intercept@b0i month@b1i;
    # dropout model
    preknot = ifelse(month == 0, 0, 1); # dummy code to activate dropout predictors post-baseline
    dropout ~ intercept@-3 (month==1) (month==2) (month==3) (month==4)
      (month==1)*med (month==2)*med (month==3)*med (month==4)*med
      dpdd*preknot dpdd.lag*preknot | intercept@0;',
  parameters = ' 
    # compute time-specific effect sizes
    d_month0 = ( (b0+b2) - (b0) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month1 = ( (b0+b2 + 1*(b1+b3)) - (b0 + 1*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month2 = ( (b0+b2 + 2*(b1+b3)) - (b0 + 2*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month3 = ( (b0+b2 + 3*(b1+b3)) - (b0 + 3*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);
    d_month4 = ( (b0+b2 + 4*(b1+b3)) - (b0 + 4*b1) ) / sqrt(dpdd.totalvar + b0i.totalvar);',
  seed = 90291, # random number seed
  chains = 4, # number of mcmc processes
  burn = 20000, # burn-in iterations
  iter = 20000, # iterations for analysis summaries
  nimps = 20 # save imputations for graphing trajectories
)

# print output
output(model8_selection)

# plot model-predicted means
bivariate_plot(dpdd.predicted ~ month | med, 
               model = model8_selection, 
               discrete_x = 'month', 
               points = F, ci = F) + ylim(0,15)

#------------------------------------------------------------------------------#
# SUMMARIZE RESULTS ----
#------------------------------------------------------------------------------#

# function to extract key estimates
extract_params <- function(object, method) {
  
  tab <- object@estimates
  
  rows <- c(
    "b0i ~ Intercept",
    "b1i ~ Intercept",
    "b0i ~ med",
    "b1i ~ med",
    "b0i residual variance",
    "b1i residual variance",
    "Cor( b0i, b1i )",
    "dpdd residual variance",
    "Parameter: d_month0",
    "Parameter: d_month1",
    "Parameter: d_month2",
    "Parameter: d_month3",
    "Parameter: d_month4",
    "dropout R2: Coefficients"
  )
  
  res <- do.call(rbind, lapply(rows, function(r) {
    if (r %in% rownames(tab)) {
      round(tab[r, c("Estimate", "StdDev"), drop = FALSE], 2)
    } else {
      data.frame(Estimate = NA_real_, StdDev = NA_real_, row.names = r)
    }
  }))
  
  rownames(res) <- c(
    "Icept (M = 0)",
    "Slope (M = 0)",
    "Icept Diff.",
    "Slope Diff.",
    "Var(Icept)",
    "Var(Slope)",
    "Cor(Icept, Slope)",
    "Var(Residual)",
    "Std. Diff. T0",
    "Std. Diff. T1",
    "Std. Diff. T2",
    "Std. Diff. T3",
    "Std. Diff. T4",
    "Pseudo-Rsq"
  )
  
  colnames(res) <- c(
    paste0("Est_", method),
    paste0("SD_", method)
  )
  
  res
}

# main summary table ----
table_summary <- cbind(
  extract_params(model1_cmar, "Mod1"),
  extract_params(model2_aux, "Mod2"),
  extract_params(model3_sharedparam, "Mod3"),
  extract_params(model5_sharedparam, "Mod5"),
  extract_params(model6_selection, "Mod6"),
  extract_params(model7_selection, "Mod7"),
  extract_params(model8_selection, "Mod8")
)
table_summary

# changes in SE units ----

est_cmar <- table_summary[, "Est_Mod1"]
se_cmar  <- table_summary[, "SE_Mod1"]

compare_methods <- c("Mod2", "Mod3", "Mod5", "Mod6", "Mod7", "Mod8")

table_change <- sapply(compare_methods, function(m) {
  round((table_summary[, paste0("Est_", m)] - est_cmar) / se_cmar, 2)
})
table_change <- as.data.frame(table_change)
rownames(table_change) <- rownames(table_summary)
table_change

# diagnostics table ----

# function to extract convergence diagnostics
extract_convergence <- function(object, method) {
  
  neff <- object@estimates[, "N_Eff"]
  neff <- neff[is.finite(neff)]
  
  psr_row <- as.numeric(object@psr[20, ])
  psr_row <- psr_row[is.finite(psr_row)]
  
  data.frame(
    Min_Neff = round(min(neff),    3),
    Max_Neff = round(max(neff),    3),
    Min_PSR  = round(min(psr_row), 3),
    Max_PSR  = round(max(psr_row), 3),
    row.names = method
  )
}

# build table
table_diag <- rbind(
  extract_convergence(model1_cmar, "Model 1"),
  extract_convergence(model2_aux, "Model 2"),
  extract_convergence(model3_sharedparam, "Model 3"),
  extract_convergence(model5_sharedparam, "Model 5"),
  extract_convergence(model6_selection, "Model 6"),
  extract_convergence(model7_selection, "Model 7"),
  extract_convergence(model8_selection, "Model 8")
)

# add number of iterations
table_diag$Iterations <- c(
  nrow(model1_cmar@iterations),
  nrow(model2_aux@iterations),
  nrow(model3_sharedparam@iterations),
  nrow(model5_sharedparam@iterations),
  nrow(model6_selection@iterations),
  nrow(model7_selection@iterations),
  nrow(model8_selection@iterations)
)
table_diag







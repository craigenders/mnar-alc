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
# CMAR MODEL ----
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
ymax <- .4

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
  burn = 100000, # burn-in iterations
  iter = 100000, # iterations for analysis summaries
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
  burn = 50000, # burn-in iterations
  iter = 50000, # iterations for analysis summaries
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
# EXTRACT ESTIMATES ----
#------------------------------------------------------------------------------#

# function to extract key estimates
extract_growth_params <- function(object, method) {
  
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
  extract_growth_params(model1, "M1"),
  extract_growth_params(model2, "M2"),
  extract_growth_params(model3, "M3"),
  extract_growth_params(model4, "M4"),
  extract_growth_params(model5, "M5"),
  extract_growth_params(model6, "M6"),
  extract_growth_params(model7, "M7"),
  extract_growth_params(model8, "M8")
)
table_summary

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
  extract_convergence(model1, "Model 1"),
  extract_convergence(model2, "Model 2"),
  extract_convergence(model3, "Model 3"),
  extract_convergence(model4, "Model 4"),
  extract_convergence(model5, "Model 5"),
  extract_convergence(model6, "Model 6"),
  extract_convergence(model7, "Model 7"),
  extract_convergence(model8, "Model 8")
)

# add number of iterations
table_diag$Iterations <- c(
  nrow(model1@iterations),
  nrow(model2@iterations),
  nrow(model3@iterations),
  nrow(model4@iterations),
  nrow(model5@iterations),
  nrow(model6@iterations),
  nrow(model7@iterations),
  nrow(model8@iterations)
)
table_diag







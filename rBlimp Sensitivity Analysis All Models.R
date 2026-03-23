#------------------------------------------------------------------------------#
# LOAD R PACKAGES ----
#------------------------------------------------------------------------------#

options (scipen = 999)

# load packages
library(ggplot2)
library(patchwork)
library(rblimp)
set_blimp('/applications/blimp/blimp-nightly')
# remotes::update_packages('rblimp')

#------------------------------------------------------------------------------#
# READ DATA ----
#------------------------------------------------------------------------------#

# github url for raw data
filepath1 <- 'https://raw.githubusercontent.com/craigenders/mnar-mlm/main/growth-dropout.csv'

# create data frame from github data
growth <- read.csv(filepath1, stringsAsFactors = T)
growth <- growth[!is.na(growth$m),]


# plotting functions
source('https://raw.githubusercontent.com/blimp-stats/blimp-book/main/misc/functions.R')
source('https://raw.githubusercontent.com/craigenders/mnar-mlm/main/mnar-plotting.R')

# compute vars
growth$med <- growth$group
growth$month <- growth$time
growth$dpdd <- growth$y - min(growth$ycom, na.rm = TRUE)
mnar <- growth[growth$id < 166,c('id','med','month','dpdd')]
write.csv(mnar, file = '~/desktop/mnar.csv', row.names = F)

#------------------------------------------------------------------------------#
# CMAR MODELS ----
#------------------------------------------------------------------------------#

# cmar analysis model 1 ----
cmar_mod1 <- rblimp(
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
output(cmar_mod1)

# plot model-predicted means
bivariate_plot(dpdd.predicted ~ month | med, model = cmar_mod1, discrete_x = 'month', points = F, ci = F) + ylim(0,15)

#------------------------------------------------------------------------------#
# TIME-RELATED DROPOUT ----
#------------------------------------------------------------------------------#

# dummy codes with group by time interaction
dropout_time <- rblimp(
  data = mnar,
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable 
  timeid = 'month', # time identifier for dropout and lagged effects 
  dropout = 'dropout = dpdd', # create dropout indicator
  # (month == 1) (month == 2) etc. are dummy codes for time
  model = 'dropout ~ intercept@-3 (month==1) (month==2) (month==3) (month==4)
      (month==1)*med (month==2)*med (month==3)*med (month==4)*med | intercept@0',
  seed = 90291, # random number seed
  burn = 20000, # burn-in iterations
  iter = 20000, # iterations for analysis summaries
  nimps = 20 # save imputations for graphing trajectories
)

# print output
output(dropout_time)

# observed dropout rates
aggregate(dropout ~ month + med, data = dropout_time@average_imp, FUN = mean)

# model-predicted dropout rates
bivariate_plot(dropout.1.probability ~ month | med, model = dropout_time, discrete_x = 'month', points = F, ci = F) + ylim(0,.40)

#------------------------------------------------------------------------------#
# SHARED PARAMETER MODELS ----
#------------------------------------------------------------------------------#

# shared parameter model 3: intercepts and slopes predicting dropout ----
sp_mod3 <- rblimp(
  data = mnar,
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable
  timeid = 'month', # time identifier for dropout and lagged effects 
  dropout = 'dropout = dpdd', # create dropout indicator
  # transform = 'preknot = ifelse(month == 0, 0, 1)',
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
output(sp_mod3)

# plot model-predicted means
bivariate_plot(dpdd.predicted ~ month | med, model = sp_mod3, discrete_x = 'month', points = F, ci = F) + ylim(0,15)

# shared parameter model 3: intercepts and slopes and their squares predicting dropout ----
sp_mod5 <- rblimp(
  data = mnar,
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable
  timeid = 'month', # time identifier for dropout and lagged effects 
  dropout = 'dropout = dpdd', # create dropout indicator
  # transform = 'preknot = ifelse(month == 0, 0, 1)',
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
output(sp_mod5)

# plot model-predicted means
bivariate_plot(dpdd.predicted ~ month | med, model = sp_mod5, discrete_x = 'month', points = F, ci = F) + ylim(0,15)

#------------------------------------------------------------------------------#
# SELECTION MODELS MODELS ----
#------------------------------------------------------------------------------#

# diggle-kenward selection model 6: current (missing) drinking predicting dropout ----
dk_mod6 <- rblimp(
  data = mnar,
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable
  timeid = 'month', # time identifier for dropout and lagged effects 
  dropout = 'dropout = dpdd', # create dropout indicator
  # transform = 'preknot = ifelse(month == 0, 0, 1)',
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
output(dk_mod6)

# plot model-predicted means
bivariate_plot(dpdd.predicted ~ month | med, model = dk_mod6, discrete_x = 'month', points = F, ci = F) + ylim(0,15)

# diggle-kenward selection model 7: current (missing) drinking and its square predicting dropout ----
dk_mod7 <- rblimp(
  data = mnar,
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable
  timeid = 'month', # time identifier for dropout and lagged effects 
  dropout = 'dropout = dpdd', # create dropout indicator
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
output(dk_mod7)

# plot model-predicted means
bivariate_plot(dpdd.predicted ~ month | med, model = dk_mod7, discrete_x = 'month', points = F, ci = F) + ylim(0,15)

# diggle-kenward selection model 8: current (missing) and prior (observed) drinking predicting dropout ----
dk_mod8 <- rblimp(
  data = mnar,
  # transform = 'preknot = ifelse(month == 0, 0, 1)',
  ordinal = 'med', # define med as categorical
  clusterid = 'id', # person id variable
  timeid = 'month', # time identifier for dropout and lagged effects 
  dropout = 'dropout = dpdd', # create dropout indicator
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
output(dk_mod8)

# plot model-predicted means
bivariate_plot(dpdd.predicted ~ month | med, model = dk_mod8, discrete_x = 'month', points = F, ci = F) + ylim(0,15)







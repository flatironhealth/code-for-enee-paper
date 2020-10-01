# Simulations for ``A note on the amount of information borrowed from external data in hybrid controlled trials with time-to-event outcomes''
# by Brian D. Segal and W. Katherine Tan

# Simulation parameters

# Used in both exponential and Weibull simulations ----------------------------

B <- 1000  # number of iterations per beta_E
BETA_E_VEC <- log(seq(0.4, 1, 0.2))  # log HR of experimental effect
LAMBDA <- 1 / 12  # hazard in control arm
EPSILON <- 0.0001  # step size for finite difference
MIN_EXTERNAL_CONTROLS <- 10  # minimum number of patients in external control arm
PROP_INFO_AT_INTERIM <- 0.33  # proportion of patients who experience event before interim is triggered

# minimum number of events per arm before triggering interim analysis
MIN_EVENTS_PER_ARM <- 10

PATIENTS_PER_MONTH <- 2  # number of patients enrolled per month
RANDOMIZATION_RATIO <- 1  # randomization ratio (#experimental / #control)

PALETTE <- c("black", "grey")  # color palette for plots

D_EFF_INTERVAL_UPPER <- 1000  # upper bound of search interval
D_EFF_INTERVAL_LOWER <- 0.001 # added to -d_C to get lower bound of search interval

# Min and max trial size
MIN_TRIAL <- 60
MAX_TRIAL <- 100

# Min and max censoring rate
MIN_CENSOR <- 0.05
MAX_CENSOR <- 0.1

PAPER_PATH <- "../paper/"  # path to paper for saving plots


# Used only in Weibull simulations --------------------------------------------

# Min and max shape parameter
SHAPE <- 1.15

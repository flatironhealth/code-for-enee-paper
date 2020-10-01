# Simulations for ``A note on the amount of information borrowed from external data in hybrid controlled trials with time-to-event outcomes''
# by Brian D. Segal and W. Katherine Tan

# functions

GetDiffInPrecision <- function(d_eff,
                               prec_hyb,
                               trial_dat,
                               regression_type = "cox",
                               time = "time",
                               event = "event",
                               experimental = "experimental",
                               covariates = NULL) {
  # This function computes the difference in precision between a model
  # with d_eff effectively added events and the precision from the hybrid model.
  # Typically called by eneeRoot().
  #
  # Arguments:
  #   d_eff (num): A proposed effective number of events added
  #   prec_hyb (num): Precision from a hybrid model that borrows from external data
  #   trial_dat (data.frame): trial data containing columns:
  #                             time (num): time of event or censoring
  #                             event (logi or num): event indicator
  #                             experimental (logi or num): experimental indicator
  #   time (chr): name of time column in trial_dat
  #   event (chr): name of event column in trial_dat
  #   experimental (chr): name of experimental column in trial_dat
  #   regression_type (chr): Type of regression model for trial data. Must be either "cox" or an acceptable argument to `dist` in survival::survreg.
  #
  # Returns:
  #   dif (num): Difference between precision from the up-weighted data and 
  #              precision from the hybrid model

  d_control <- sum(trial_dat[, event] == 1 & trial_dat[, experimental] == 0)
  control_weight <- d_eff / d_control + 1
  trial_dat$weights <- ifelse(trial_dat[, experimental] == 0, control_weight, 1)
  
  if(is.null(covariates)) {

    formula <- as.formula(paste0("Surv(", time, ", ", event, ") ~ ", experimental))

  } else {

    formula <- as.formula(paste0("Surv(", time, ", ", event, ") ~ ", experimental, " + ", paste(covariates, collapse = " + ")))

  }

  if (regression_type == "cox") {

    fit <- coxph(formula, weights = weights, data = trial_dat)

  } else {

    fit <- survreg(formula, weights = weights, dist = regression_type, data = trial_dat)
  
  }

  prec_d_eff <- 1 / vcov(fit)[experimental, experimental]

  dif <- prec_d_eff - prec_hyb
  return(dif)
}

GetDiffInPrecision <- Vectorize(GetDiffInPrecision, vectorize.args = "d_eff")

eneeRoot <- function(trial_data,
                     prec_hyb,
                     d_eff_interval = NULL,
                     d_eff = NULL,
                     regression_type = "cox",
                     find_root = TRUE,
                     time = "time",
                     event = "event",
                     experimental = "experimental",
                     covariates = NULL) {
  # This function computes the effective number of external events (ENEE) when find_root == TRUE,
  # and the weighted precision of the experimental effect when find_root == FALSE. This is the
  # general solution proposed by Segal and Tan.
  #
  # Arguments:
  #   trial_dat (data.frame): trial data containing columns:
  #                             time (num): time of event or censoring
  #                             event (logi or num): event indicator
  #                             experimental (logi or num): experimental indicator
  #   prec_hyb (num): Precision from a hybrid model that borrows from external data
  #   d_eff_interval (num): interval over which to search for a solution for d_eff
  #   d_eff (num): number of additional control events at which to evaluate the precision
  #   regression_type (chr): Type of regression model for trial data. Must be either "cox" or an acceptable argument to `dist` in survival::survreg.
  #   find_root (log): If TRUE, solves for d_eff. If FALSE, returns value of weighted precision
  #   time (chr): name of time column in trial_dat
  #   event (chr): name of event column in trial_data
  #   experimental (chr): name of experimental column in trial_dat
  #   covariates (chr): vector of covariates to adjust for (can include interaction terms). All covariates must be in trial_dat.
  #
  # Returns:
  #   out (num): numeric vector length 1 if find_root == TRUE, and length of d_eff if find_root = FALSE.
  #              If estimation fails, returns NA.

  if (!is.null(covariates)) {
    unique_cols <- unique(unlist(strsplit(covariates, "\\*|:")))
    if(!all(unique_cols %in% colnames(trial_data))) {
      stop("all covariates must be columns in trial_data")
    }
  }

  if (!regression_type %in% c("cox", names(survreg.distributions))) {
    stop(paste0("regression_type must be either cox or a distribution defined in survreg.distributions (",
                paste0(names(survreg.distributions), collapse = ", "), ")"))
  }

  if (find_root & (length(d_eff_interval) != 2 | !is.numeric(d_eff_interval))) {
    stop("d_eff_interval must be a numeric vector of length 2")
  }

  if (!find_root & (is.null(d_eff) | !is.numeric(d_eff))) {
    stop("d_eff must be specified as a numeric vector if find_root == FALSE")
  }

  if (find_root & !is.null(d_eff)) {
    warning("find_root == TRUE so d_eff will be ignored")
  }

  if(find_root) {

    out <- NA

    try(out <- uniroot(GetDiffInPrecision,
                       interval = d_eff_interval,
                       prec_hyb = prec_hyb,
                       trial_dat = trial_data,
                       regression_type = regression_type,
                       covariates = covariates)$root)

  } else {

    out <- rep(NA, length(d_eff))

    if(!any(is.na(d_eff))) {
      try(out <- GetDiffInPrecision(d_eff = d_eff,
                                    prec_hyb = prec_hyb,
                                    trial_dat = trial_data,
                                    regression_type = regression_type,
                                    covariates = covariates) + prec_hyb)
    }
  }

  return(out)
}

eneeExact <- function(d_C, d_E, prec_hyb) {
  # This function calculates the exact effective number of external events (ENEE)
  # for unadjusted exponential models.
  #
  # Arguments:
  #   d_C (num): number of events in the randomized control arm
  #   d_E (num): number of events in the randomized experimental arm
  #   prec_hyb (num): Precision from a hybrid model that borrows from external data
  # 
  # Returns:
  #   out (num): ENEE
  out <- ((prec_hyb * (d_C + d_E)) - d_C * d_E) / (d_E - prec_hyb)
  return(out)

}

ehss <- function(n_C, n_E, prec_hyb, prec_ref) {
  # This function calculates the linear approximation for effective historical sample size (EHSS).
  # proposed by Hobbs et al. (2013).
  # If the number of events d_C and d_E are passed to n_C and n_E, respectively, then this function returns EHSS_d as described by Segal and Tan
  #
  # Arguments:
  #   n_C (num): number of events in the randomized control arm
  #   n_E (num): number of events in the randomized experimental arm
  #   prec_hyb (num): Precision from a hybrid model that borrows from external data
  #   prec_ref (num): Precision from the reference model that does not borrow from external data
  # 
  # Returns:
  #   out (num): EHSS

  out <- (n_C + n_E) * (prec_hyb / prec_ref - 1)
  return(out)

}

SimulateTrialData <- function(n_trial,
                              prop_info_at_interim,
                              h_C,
                              beta_E,
                              p_censor,
                              dist = "exponential",
                              ratio = 1,
                              shape = 1,
                              patients_per_month = 2,
                              min_events_per_arm = 10) {
  # This function simulates data for a 2-arm trial taking into account 
  # enrollment patterns. Assumes that h_C is the baseline hazard in terms 
  # of months.
  #
  # Arguments:
  #   n_trial (num): number of patients in trial
  #   prop_info_at_interim (num): Proportion of patients in trial with an event that will trigger interim analysis
  #   h_C (num): hazard in control arm
  #   beta_E (num): log HR of experimental effect
  #   p_censor (num): probability of censoring for exponential data
  #   dist (chr): specifies generating distribution (either "exponential" or "weibull")
  #   ratio (num): specifies randomization ratio (#experimental / #control). If not already an integer, will be rounded down to the nearest integer.
  #   shape (num): shape parameter of Weibull distribution
   #   patients_per_month (num): Number of patients that enroll per month (assumed linear)
   #   min_events_per_arm (num): Minimum number of events required in each arm prior to interim
  #
  # Returns:
  #   sim_data (data.frame, n_trial rows 4 columns): simulated data with columns:
  #     calendar_time (num): Calendar time since start of trial that event or censoring occurs
  #     time (num): Time from enrollment to event
  #     event (logi):  TRUE if event observed, FALSE if censored
  #     experimental (num): Equal to 1 if patient is in experimental group, 0 if in control

  enrollment_pattern <- rep(c(0, 1), times = c(1, ratio))

  sim_data <- data.frame(experimental = rep(enrollment_pattern, length.out = n_trial),
                         month_enter = (1:n_trial) / patients_per_month,
                         external = 0)

  # generate event times
  if (dist == "exponential") {

    t_trial <- with(sim_data, rexp(n = n_trial, rate = h_C * exp(beta_E)^experimental))
    c_trial <- with(sim_data, rexp(n = n_trial, rate = h_C * exp(beta_E)^experimental * p_censor / (1 - p_censor)))

  } else if (dist == "weibull") {

    t_trial <- with(sim_data, rweibull(n = n_trial,
                      scale = 1 / (h_C * exp(beta_E)^experimental), 
                      shape = shape))

    # Note: This won't give the expected probability of censoring for Weibull data,
    # but that's ok -- it will give some censoring rate,
    # and the censoring rate is not used to present or interpret results
    c_trial <- with(sim_data, rweibull(n = n_trial,
                      scale = 1 / (h_C * exp(beta_E)^experimental * p_censor / (1 - p_censor)),
                      shape = shape))

  } else {

    stop("dist must be equal to either 'exponential' or 'weibull'")

  }

  sim_data$time <- pmin(t_trial, c_trial)
  sim_data$event <- (t_trial <= c_trial)

  # get calendar time after start of study at which event occurs
  sim_data$calendar_time <- with(sim_data, time + month_enter)

  # order patients by calendar time of event
  sim_data <- sim_data[order(sim_data$calendar_time), ]

  # censor events that occur after date of data freeze for interim
  # In addition to target amount of information, 
  # Ensure that there is a minimum number of events in each arm to prevent
  # pathological simulations
  rows_to_censor <- cumsum(sim_data$event) >= floor(prop_info_at_interim * n_trial) &
                    with(sim_data, cumsum(event * experimental) > min_events_per_arm) &
                    with(sim_data, cumsum(event * !experimental) > min_events_per_arm)

  min_row_to_censor <- ifelse(sum(rows_to_censor) > 0, min(which(rows_to_censor)), NA)

  if (!is.na(min_row_to_censor)) {
    calendar_time_at_interim <- sim_data$calendar_time[min_row_to_censor]
    sim_data$event[rows_to_censor] <- FALSE
    sim_data$time[rows_to_censor] <- with(sim_data[rows_to_censor, ],
                                          pmin(time, calendar_time_at_interim - month_enter))

  # Remove patients who wouldn't have enrolled until after interim
  sim_data <- sim_data[sim_data$month_enter <= calendar_time_at_interim, ]
  }

  return(sim_data[, c("calendar_time", "time", "event", "experimental")])
}

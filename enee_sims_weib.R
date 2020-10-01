# Simulations for ``A note on the amount of information borrowed from external data in hybrid controlled trials with time-to-event outcomes''
# by Brian D. Segal and W. Katherine Tan

# Weibull simulations

library(survival)
library(ggplot2)
library(dplyr)
library(reshape2)
library(scales)

source("enee_functions.R")
source("enee_parameters.R")

# Labels plot legends (will be rendered with scales::parse_format())
ENEE_LABEL <- "d[eff] (5)"
EHSS_LABEL <- "EHSS[d] (6)"

# Heuristic cutoff for defining ENEE stability
DERIV_CUT_OFF <- exp(-3)

# data frame to hold results
results <- data.frame(beta_E = BETA_E_VEC[rep(1:length(BETA_E_VEC), each = B)]) %>%
  mutate(shape = NA,
         prec_ref = NA,
         prec_hyb = NA,
         beta_ref = NA,
         beta_hyb = NA,
         n_trial = NA,
         n_control = NA,
         n_experimental = NA,
         n_external = NA,
         p_censor = NA,
         d_C = NA,
         d_E = NA,
         d_hyb = NA,
         enee = NA, 
         enee_deriv = NA)

set.seed(123)

for (b in 1:nrow(results)) {

  print(paste(b, "of", nrow(results)))

  results$n_trial[b] <- sample(MIN_TRIAL:MAX_TRIAL, size = 1)
  results$p_censor[b] <- runif(n = 1, min = MIN_CENSOR, max = MAX_CENSOR)
  results$shape[b] <- SHAPE
  # runif(n = 1, min = MIN_SHAPE, max = MAX_SHAPE)

  # generate trial data
  trial_data <- SimulateTrialData(n_trial = results$n_trial[b],
                                  prop_info_at_interim = PROP_INFO_AT_INTERIM,
                                  shape = results$shape[b],
                                  h_C = LAMBDA,
                                  dist = "weibull",
                                  beta_E = results$beta_E[b],
                                  p_censor = results$p_censor[b],
                                  ratio = RANDOMIZATION_RATIO,
                                  patients_per_month = PATIENTS_PER_MONTH,
                                  min_events_per_arm = MIN_EVENTS_PER_ARM)

  results$d_C[b] <- with(trial_data, sum(event[experimental == 0]))
  results$d_E[b] <- with(trial_data, sum(event[experimental == 1]))
  results$n_control[b] <- sum(trial_data$experimental == 0)
  results$n_experimental[b] <- sum(trial_data$experimental == 1)

  # At most, select no more external controls than trial controls.
  results$n_external[b] <- sample(MIN_EXTERNAL_CONTROLS:results$n_control[b], size = 1)

  max_trial_followup <- max(trial_data$time)

  # Obtain external data by simulating new data
  external_data <- SimulateTrialData(n_trial = results$n_external[b],
                                     prop_info_at_interim = 1,
                                     shape = results$shape[b],
                                     h_C = LAMBDA,
                                     dist = "weibull",
                                     beta_E = 0,  # both control and experimental from same dist'n
                                     p_censor = results$p_censor[b],
                                     ratio = RANDOMIZATION_RATIO,
                                     patients_per_month = PATIENTS_PER_MONTH) %>%
    # set all patients to be controls and censor at max trial follow-up
    mutate(experimental = 0, 
           event = ifelse(time <= max_trial_followup, event, FALSE),
           time = ifelse(time <= max_trial_followup, time, max_trial_followup))

  results$d_hyb[b] <- sum(external_data$event)

  hyb_data <- rbind(trial_data, external_data)

  # fit models ----------------------------------------------------------------
  ref <- coxph(Surv(time, event) ~ experimental, data = trial_data)

  hyb <- coxph(Surv(time, event) ~ experimental, data = hyb_data)

  results$prec_ref[b] <- 1 / ref$var
  results$prec_hyb[b] <- 1 / hyb$var

  results$beta_ref[b] <- coef(ref)
  results$beta_hyb[b] <- coef(hyb)

  d_eff_interval <- c(-results$d_C[b] + D_EFF_INTERVAL_LOWER, D_EFF_INTERVAL_UPPER)
  
  # ENEE
  results$enee[b] <- eneeRoot(trial_data,
                              prec_hyb = results$prec_hyb[b],
                              d_eff_interval = d_eff_interval,
                              regression_type = "cox")

# Central finite difference
  curve <- eneeRoot(trial_data,
                    prec_hyb = results$prec_hyb[b],
                    d_eff = results$enee[b] + c(-EPSILON, EPSILON),
                    regression_type = "cox",
                    find_root = FALSE)

  results$enee_deriv[b] <- diff(curve) / (2 * EPSILON)

}

# calculate EHSS
results <- results %>% 
  mutate(ehss_d = ehss(d_C, d_E, prec_hyb, prec_ref))

# plots and tables ------------------------------------------------------------

# view distribution of derivatives
hist(log(results$enee_deriv))

# Number of unstable results
n_na <- sum(is.na(results$enee))
n_deriv_too_small <- sum(results$enee_deriv < DERIV_CUT_OFF, na.rm = TRUE)

n_unstable <- c(na = n_na,
                deriv_too_small = n_deriv_too_small,
                total = n_na + n_deriv_too_small)

n_unstable
n_unstable["total"] / nrow(results)

# prepare data for plots and tables
results_melted <- results %>%
  mutate(enee_status = ifelse(enee_deriv < DERIV_CUT_OFF | is.na(enee), "unstable", "stable"),
         enee_status = factor(enee_status, levels = c("stable", "unstable"))) %>%
  melt(id.vars = c("beta_E", "d_hyb", "enee_status"),
                   measure.vars = c("ehss_d", "enee")) %>%
  mutate(variable = ifelse(variable == "ehss_d", EHSS_LABEL, ENEE_LABEL),
         variable = factor(variable, levels = c(ENEE_LABEL, EHSS_LABEL)))

# Table of bias and variance
results_melted %>%
  mutate(value = ifelse(enee_status == "unstable" & variable == ENEE_LABEL, NA, value)) %>%
  group_by(enee_status, HR = exp(beta_E), variable) %>%
  summarise(n = sum(!is.na(value)),
            mse = sum(mean((value - d_hyb)^2)),
            bias = sum(mean(value - d_hyb)),
            sd = sqrt(mse - bias^2)) %>%
  select(-mse) %>%
  filter(!is.na(bias))

# removing unstable ENEE results, keeping all EHSS results
results_melted %>%
  mutate(value = ifelse(enee_status == "unstable" & variable == ENEE_LABEL, NA, value)) %>%
  ggplot(aes(x = d_hyb, y = value, color = variable, shape = variable),
       data = .) +
    geom_point(size = 2) +
    theme_bw(20) +
    facet_grid(enee_status ~ beta_E,
               switch = "y",
               labeller = label_bquote(cols = paste("HR = ", .(signif(exp(beta_E), 2))),
                                       rows = paste(.(levels(enee_status)[enee_status]), " ", d[eff]))) +
    scale_color_manual(values = PALETTE, labels = parse_format()) +
    scale_shape_discrete(labels = parse_format()) +
    scale_x_continuous(lim = c(0, max(results_melted$d_hyb))) +
    scale_y_continuous(lim = c(0, max(results_melted$value))) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(x = "True number of events borrowed",
         y = "Estimate") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.title = element_blank(),
          legend.position = "bottom")
ggsave(file.path(PAPER_PATH, "EHSS_ENEE_weibull_cox.png"), width = 8, height = 5.5)

# All results
results_melted %>%
  ggplot(aes(x = d_hyb, y = value, color = variable, shape = variable),
       data = .) +
    geom_point(size = 2) +
    theme_bw(20) +
    facet_grid(. ~ beta_E,
               scales = "free_y",
               switch = "y",
               labeller = label_bquote(cols = paste("HR = ", .(signif(exp(beta_E), 2))))) +
    scale_color_manual(values = PALETTE, labels = parse_format()) +
    scale_shape_discrete(labels = parse_format()) +
    scale_x_continuous(lim = c(0, max(results_melted$d_hyb))) +
    scale_y_continuous(lim = c(0, max(results_melted$value))) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(x = "True number of events borrowed",
         y = "Estimate") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.title = element_blank(),
          legend.position = "bottom")
ggsave(file.path(PAPER_PATH, "EHSS_ENEE_weibull_cox_all.png"), width = 8, height = 4)

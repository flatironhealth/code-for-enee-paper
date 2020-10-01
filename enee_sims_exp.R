# Simulations for ``A note on the amount of information borrowed from external data in hybrid controlled trials with time-to-event outcomes''
# by Brian D. Segal and W. Katherine Tan

# exponential simulations

library(survival)
library(ggplot2)
library(dplyr)
library(reshape2)
library(scales)

source("enee_functions.R")
source("enee_parameters.R")

# Labels plot legends (will be rendered with scales::parse_format())
ENEE_LABEL <- "d[eff] (4)"
EHSS_LABEL <- "EHSS[d] (6)"

# data frame to hold results
results <- data.frame(beta_E = BETA_E_VEC[rep(1:length(BETA_E_VEC), each = B)]) %>%
  mutate(prec_ref = NA,
         prec_hyb = NA,
         beta_ref = NA,
         beta_hyb = NA,
         n_trial = NA,
         n_external = NA,
         n_control = NA,
         n_experimental = NA,
         p_censor = NA,
         d_C = NA,
         d_E = NA,
         d_hyb = NA,
         enee = NA)

set.seed(123)

for (b in 1:nrow(results)) {

  print(paste(b, "of", nrow(results)))

  results$n_trial[b] <- sample(MIN_TRIAL:MAX_TRIAL, size = 1)
  results$p_censor[b] <- runif(n = 1, min = MIN_CENSOR, max = MAX_CENSOR)

  # generate trial data
  trial_data <- SimulateTrialData(n_trial = results$n_trial[b],
                                  prop_info_at_interim = PROP_INFO_AT_INTERIM,
                                  h_C = LAMBDA,
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
                                     h_C = LAMBDA,
                                     dist = "exponential",
                                     beta_E = 0,  # both control and experimental from same distribution
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
  ref <- survreg(Surv(time, event) ~ experimental, dist = "exponential", data = trial_data)

  hyb <- survreg(Surv(time, event) ~ experimental, dist = "exponential", data = hyb_data)

  results$prec_ref[b] <- 1 / ref$var["experimental", "experimental"]
  results$prec_hyb[b] <- 1 / hyb$var["experimental", "experimental"]

  results$beta_ref[b] <- -coef(ref)["experimental"]
  results$beta_hyb[b] <- -coef(hyb)["experimental"]

}

# calculate ENEE and EHSS
results <- results %>% 
  mutate(enee = eneeExact(d_C, d_E, prec_hyb),
         ehss_d = ehss(d_C, d_E, prec_hyb, prec_ref))

# plots -----------------------------------------------------------------------

# prepare data
results_melted <- melt(results,
                       id.vars = c("beta_E", "d_hyb"),
                       measure.vars = c("ehss_d", "enee")) %>%
  mutate(variable = ifelse(variable == "ehss_d", EHSS_LABEL, ENEE_LABEL),
         variable = factor(variable, levels = c(ENEE_LABEL, EHSS_LABEL)))

ggplot(aes(x = d_hyb, y = value, color = variable, shape = variable),
       data = results_melted) +
  geom_point(size = 2) +
  theme_bw(20) +
  facet_grid( ~ beta_E,
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
ggsave(file.path(PAPER_PATH, "EHSS_ENEE_exponential_exp.png"), width = 8, height = 4)

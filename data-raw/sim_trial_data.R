set.seed(1)
n <- 500

vaccination <- rbinom(n, 1, 0.5)
age <- rnorm(n, mean = 55, sd = 15)
age <- pmax(18, pmin(age, 90))
baseline_risk <- sample(c("low", "medium", "high"), n,
                        replace = TRUE,
                        prob = c(0.3, 0.5, 0.2))

# Simulate time-to-event using Weibull distribution

beta_age <- 0.02
beta_risk_medium <- 0.5
beta_risk_high <- 1.0
beta_vaccination <- -0.7

linear_pred <- beta_age * (age - 55) +
  beta_risk_medium * (baseline_risk == "medium") +
  beta_risk_high * (baseline_risk == "high") +
  beta_vaccination * vaccination

scale <- exp(-linear_pred) * 50
event_time <- rweibull(n, shape = 1.5, scale = scale)

admin_censor <- 100
dropout_time <- rexp(n, rate = 0.01)
censor_time <- pmin(admin_censor, dropout_time)
time <- pmin(event_time, censor_time)

status <- as.numeric(event_time <= censor_time)

sim_trial_data <- data.frame(
  id = 1:n,
  vaccination = vaccination,
  age = round(age, 1),
  baseline_risk = factor(baseline_risk, levels = c("low", "medium", "high")),
  time = round(time, 2),
  status = status
)

usethis::use_data(sim_trial_data, overwrite = TRUE)

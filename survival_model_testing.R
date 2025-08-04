source("survival_object.R")
## Cox-regression
# NUll
cox1 <- coxph(Surv(time, status) ~ strata(percentile), data = survival_data)
cox2 <- coxph(Surv(time, status) ~ strata(percentile) + site, data = survival_data)
cox3 <- coxph(Surv(time, status) ~ strata(percentile) + gender, data = survival_data)
cox4 <- coxph(Surv(time, status) ~ strata(percentile) + gender + site, data = survival_data)

model_comparison <-
  rbind(
    broom::glance(cox1), broom::glance(cox2),
    broom::glance(cox3), broom::glance(cox4)) %>%
  mutate(Predictor = c("Risk", "Site", "Gender", "Gender + Site"))
  # select(Predictor,)
  # relocate()

rempsyc::nice_table(model_comparison)

summary(cox4)
fit <- survfit(cox4)

## Test Cox model assumption
# https://www.sthda.com/english/wiki/cox-model-assumptions
# Schoenfeld residuals: proportional hazards assumption
ph_test <- cox.zph(cox4)
ggcoxzph(ph_test, ggtheme = theme_publish())
# notes: if the variable is categorical, we expect to see
# two clusters of betas

# Martingale residual: nonlinearity
ggcoxdiagnostics(
  cox,
  type = "deviance",
  linear.predictions = FALSE,
  ggtheme = theme_publish())

# Deviance residual: influential observations
ggcoxfunctional(fit)

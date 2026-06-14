# Survival probabilities
source("C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

hist <-
  readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0))

wd <-
  inner_join(survival_data, hist, by = "IID") %>%
  select(-IID)

cox <- coxph(Surv(time, status) ~ strata(percentile) + any_hist + gender + site, data = wd)
fit1 <- survfit(cox)

cox <- coxph(Surv(time, status) ~ strata(percentile) + any_hist + site, data = filter(wd, gender == "Female"))
fit2 <- survfit(cox)

cox <- coxph(Surv(time, status) ~ strata(percentile) + any_hist + site, data = filter(wd, gender == "Male"))
fit3 <- survfit(cox)

sm <- summary(fit1, times = c(5, 10, 12, 15, 20, 25))
tab <-
  data.frame(
    strata   = sm$strata,
    time     = sm$time,
    n_risk   = sm$n.risk,
    n_event  = sm$n.event,
    surv     = sm$surv,
    lower_ci = sm$lower,
    upper_ci = sm$upper) %>%
  mutate(cuminc = 1 - surv) %>%
  mutate(
    estimate = sprintf(
      "%.1f%% (%.1f–%.1f)",
      cuminc * 100)) %>%
  select(strata, time, estimate) %>%
  pivot_wider(
    names_from = strata,
    values_from = estimate)

final <-
  flextable(tab) %>%
  add_header_row(
    values = c("Age (yr)", "Bottom 10%", "Else", "Top 10%")) %>%
	align(align = "center", part = "all") %>%
	autofit()

save_as_docx(
  "Supplementary Table S24 - Cumulative survival probabilities per age (all)" = final,
  path = "sup_tab24_pt1.docx")


sm <- summary(fit2, times = c(5, 10, 12, 15, 20, 25))
tab <-
  data.frame(
    strata   = sm$strata,
    time     = sm$time,
    n_risk   = sm$n.risk,
    n_event  = sm$n.event,
    surv     = sm$surv,
    lower_ci = sm$lower,
    upper_ci = sm$upper) %>%
  mutate(cuminc = 1 - surv) %>%
  mutate(
    estimate = sprintf(
      "%.1f%% (%.1f–%.1f)",
      cuminc * 100)) %>%
  select(strata, time, estimate) %>%
  pivot_wider(
    names_from = strata,
    values_from = estimate)

final <-
  flextable(tab) %>%
  add_header_row(
    values = c("Age (yr)", "Bottom 10%", "Else", "Top 10%")) %>%
	align(align = "center", part = "all") %>%
	autofit()

save_as_docx(
  "Supplementary Table S24 - Cumulative survival probabilities per age (fem)" = final,
  path = "sup_tab24_pt2.docx")

sm <- summary(fit3, times = c(5, 10, 12, 15, 20, 25))
tab <-
  data.frame(
    strata   = sm$strata,
    time     = sm$time,
    n_risk   = sm$n.risk,
    n_event  = sm$n.event,
    surv     = sm$surv,
    lower_ci = sm$lower,
    upper_ci = sm$upper) %>%
  mutate(cuminc = 1 - surv) %>%
  mutate(
    estimate = sprintf(
      "%.1f%% (%.1f–%.1f)",
      cuminc * 100)) %>%
  select(strata, time, estimate) %>%
  pivot_wider(
    names_from = strata,
    values_from = estimate)

final <-
  flextable(tab) %>%
  add_header_row(
    values = c("Age (yr)", "Bottom 10%", "Else", "Top 10%")) %>%
	align(align = "center", part = "all") %>%
	autofit()

save_as_docx(
  "Supplementary Table S24 - Cumulative survival probabilities per age (man)" = final,
  path = "sup_tab24_pt3.docx")
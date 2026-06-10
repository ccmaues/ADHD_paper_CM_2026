# Main Cox proportional hazards model (90th only stratified fig6F)
pacman::p_load(flextable, gtsummary)
source("C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

hist <-
  readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0))

wd <-
  inner_join(survival_data, hist, by = "IID") %>%
  select(-IID)

all <-
  tbl_regression(
    coxph(Surv(time, status) ~ percentile + gender + site + any_hist, data = wd),
    exponentiate = TRUE,
    conf.int = TRUE)

females <-
  tidy(
    coxph(Surv(time, status) ~ percentile + any_hist, data = filter(wd, gender == "Female")),
    exponentiate = TRUE,
    conf.int = TRUE) %>%
  mutate(group = "Females")

males <-
  tidy(
    coxph(Surv(time, status) ~ percentile + site + any_hist, data = filter(wd, gender == "Male")),
    exponentiate = TRUE,
    conf.int = TRUE) %>%
  mutate(group = "Males")

for_plot_HR <-
  rbind(all, females, males) %>%

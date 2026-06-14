# Main Cox proportional hazards model (90th only stratified fig6F)
pacman::p_load(flextable, gtsummary, survival, broom)
source("C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

hist <-
  readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0))

wd <-
  inner_join(survival_data, hist, by = "IID") %>%
  select(-IID) %>%
  filter(percentile == "90th")

all_90 <-
  coxph(
    Surv(time, status) ~ site + gender + any_hist,
    data = wd)

females_90 <-
  coxph(
    Surv(time, status) ~ site + any_hist,
    data = wd)

males_90 <-
  coxph(
    Surv(time, status) ~ site + any_hist,
    data = wd)

final <- tbl_merge(
    list(
      tbl_regression(all_90, exponentiate = TRUE, conf.int = TRUE), 
      tbl_regression(females_90, exponentiate = TRUE, conf.int = TRUE), 
      tbl_regression(males_90, exponentiate = TRUE, conf.int = TRUE)), 
    tab_spanner = c(
      "**All samples**",
      "**Females**",
      "**Males**")) %>%
  bold_labels() %>%
  italicize_levels() %>%
  as_flex_table()

save_as_docx(
  "Supplementary Table S22 - Main Cox proportional hazards model by gender (90th only)" = final,
  path = "sup_tab22.docx")

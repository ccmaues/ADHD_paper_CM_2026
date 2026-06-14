# Main Cox proportional hazards model (stratified by history fig6F)
pacman::p_load(flextable, gtsummary, survival, broom)
source("C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

hist <-
  readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0))

wd <-
  inner_join(survival_data, hist, by = "IID") %>%
  select(-IID)

all_hist <-
  coxph(
    Surv(time, status) ~ percentile + site + gender,
    data = filter(wd, any_hist == 1))

all_no_hist <-
  coxph(
    Surv(time, status) ~ percentile + site + gender,
    data = filter(wd, any_hist == 0))

females_hist <- 
  coxph(
    Surv(time, status) ~ percentile + site,
    data = filter(wd, gender == "Female" & any_hist == 1))

females_no_hist <-
  coxph(
    Surv(time, status) ~ percentile + site,
    data = filter(wd, gender == "Female" & any_hist == 0))

males_hist <-
  coxph(
    Surv(time, status) ~ percentile + site,
    data = filter(wd, gender == "Male" & any_hist == 1))

males_no_hist <-
  coxph(
    Surv(time, status) ~ percentile + site,
    data = filter(wd, gender == "Male" & any_hist == 0))

final <- tbl_merge(
    list(
      tbl_regression(all_hist, exponentiate = TRUE, conf.int = TRUE), 
      tbl_regression(all_no_hist, exponentiate = TRUE, conf.int = TRUE), 
      tbl_regression(females_hist, exponentiate = TRUE, conf.int = TRUE), 
      tbl_regression(females_no_hist, exponentiate = TRUE, conf.int = TRUE), 
      tbl_regression(males_hist, exponentiate = TRUE, conf.int = TRUE), 
      tbl_regression(males_no_hist, exponentiate = TRUE, conf.int = TRUE)), 
    tab_spanner = c(
      "**All samples (+hist)**",
      "**All samples (-hist)**",
      "**Females (+hist)**",
      "**Females (-hist)**",
      "**Males (+hist)**",
      "**Males (-hist)**")) %>%
  bold_labels() %>%
  italicize_levels() %>%
  as_flex_table()

save_as_docx(
  "Supplementary Table S21 - Main Cox proportional hazards model by gender (family history only)" = final,
  path = "sup_tab21.docx")

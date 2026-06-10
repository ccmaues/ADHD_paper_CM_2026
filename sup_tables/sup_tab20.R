# Main Cox proportional hazards model (fig6 B and D)
pacman::p_load(flextable, gtsummary)
source("C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

hist <-
  readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0))

wd <-
  inner_join(survival_data, hist, by = "IID") %>%
  select(-IID)

all <- coxph(Surv(time, status) ~ percentile + gender + site + any_hist, data = wd)
females <- coxph(Surv(time, status) ~ percentile + any_hist, data = filter(wd, gender == "Female"))
males <- coxph(Surv(time, status) ~ percentile + site + any_hist, data = filter(wd, gender == "Male"))

final <- tbl_merge(
    list(
      tbl_regression(all, exponentiate = TRUE, conf.int = TRUE), 
      tbl_regression(females, exponentiate = TRUE, conf.int = TRUE), 
      tbl_regression(males, exponentiate = TRUE, conf.int = TRUE)),
    tab_spanner = c("**All samples**", "**Females**", "**Males**")) %>%
  bold_labels() %>%
  italicize_levels() %>%
  as_flex_table()

save_as_docx(
  "Supplementary Table S20 - Main Cox proportional hazards model" = final,
  path = "sup_tab20.docx")

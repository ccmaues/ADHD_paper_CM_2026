# Main Cox proportional hazards model
pacman::p_load(dplyr, tidyr, ggsurvfit, survminer, survival, broom, flextable, gtsummary)

source("C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

# Family history
hist <-
  readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0))

wd <-
  inner_join(survival_data, hist, by = "IID") %>%
  select(-IID)

# ALL samples (UPDATE THIS ONE ON THE TABLE)
cox1 <- coxph(
    Surv(time, status) ~ strata(percentile) + any_hist + gender + site,
    data = wd)

final <- 
    tbl_regression(cox1) %>%
    as_flex_table()

save_as_docx(
  "Supplementary Table S20" = final,
  path = "sup_tab20.docx")

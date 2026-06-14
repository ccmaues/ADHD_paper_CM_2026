# final time-to-event model report
pacman::p_load(broom, survival, flextable)
source("C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

hist <-
  readRDS("D:/cass_HD/DD_CM_backup/cass_BHRC_28042025_ARTICLE.RDS")$family_history %>%
  mutate(any_hist = if_else(if_any(starts_with("parent_"), ~ . == 1), 1, 0))

wd <-
  inner_join(survival_data, hist, by = "IID") %>%
  select(-IID)

cox1 <- coxph(Surv(time, status) ~ strata(percentile) + any_hist + gender + site, data = wd)
cox2 <- coxph(Surv(time, status) ~ strata(percentile) + any_hist, data = filter(wd, gender == "Female"))
cox3 <- coxph(Surv(time, status) ~ strata(percentile) + any_hist + site, data = filter(wd, gender == "Male"))

final <-
	tbl_merge(
		list(
			tbl_regression(cox1, exponentiate = TRUE, conf.int = TRUE),
			tbl_regression(cox2, exponentiate = TRUE, conf.int = TRUE),
			tbl_regression(cox3, exponentiate = TRUE, conf.int = TRUE)),
		tab_spanner = c(
      "**All samples**",
      "**Females**",
      "**Males**")) %>%
  bold_labels() %>%
  italicize_levels() %>%
  as_flex_table()

save_as_docx(
  "Supplementary Table S23 - Main Cox proportional hazards report" = final,
  path = "sup_tab23.docx")


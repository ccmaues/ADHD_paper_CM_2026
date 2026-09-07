# Final time-to-event model report
# Supplementary Table S23

pacman::p_load(
	dplyr,
	survival,
	gtsummary,
	flextable)

source("C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

# --------------------------------------------------
# Working data
# --------------------------------------------------

wd <-
	survival_data %>%
	inner_join(
		select(
			hist,
			IID,
			any_hist),
		by = "IID") %>%
	select(-IID)

# --------------------------------------------------
# Cox proportional hazards models
# --------------------------------------------------

cox1 <-
	coxph(
		Surv(time, status) ~
			strata(percentile) +
			any_hist +
			gender +
			site,
		data = wd)

cox2 <-
	coxph(
		Surv(time, status) ~
			strata(percentile) +
			any_hist,
		data = filter(
			wd,
			gender == "Female"))

cox3 <-
	coxph(
		Surv(time, status) ~
			strata(percentile) +
			any_hist +
			site,
		data = filter(
			wd,
			gender == "Male"))

# --------------------------------------------------
# Final table
# --------------------------------------------------

final <-
	tbl_merge(
		tbls = list(
			tbl_regression(
				cox1,
				exponentiate = TRUE,
				conf.int = TRUE),

			tbl_regression(
				cox2,
				exponentiate = TRUE,
				conf.int = TRUE),

			tbl_regression(
				cox3,
				exponentiate = TRUE,
				conf.int = TRUE)),

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

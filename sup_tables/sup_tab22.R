# Main Cox proportional hazards model
# Top 10% PRS only
# Supplementary Table S22

pacman::p_load(
	dplyr,
	survival,
	flextable,
	gtsummary)

source(
	"C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

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
	filter(
		percentile == "90th")

# --------------------------------------------------
# Cox proportional hazards models
# Top 10% PRS only
# --------------------------------------------------

all_90 <-
	coxph(
		Surv(time, status) ~
			site +
			gender +
			any_hist,
		data = wd)

females_90 <-
	coxph(
		Surv(time, status) ~
			site +
			any_hist,
		data = filter(
			wd,
			gender == "Female"))

males_90 <-
	coxph(
		Surv(time, status) ~
			site +
			any_hist,
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
				all_90,
				exponentiate = TRUE,
				conf.int = TRUE),

			tbl_regression(
				females_90,
				exponentiate = TRUE,
				conf.int = TRUE),

			tbl_regression(
				males_90,
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
	"Supplementary Table S22 - Hazard ratio by gender (90th only)" = final,
	path = "sup_tab22.docx")
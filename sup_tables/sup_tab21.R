# Main Cox proportional hazards model
# Stratified by family history
# Supplementary Table S21

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
	select(-IID)

# --------------------------------------------------
# All samples
# --------------------------------------------------

all_hist <-
	coxph(
		Surv(time, status) ~
			percentile +
			site +
			gender,
		data = filter(
			wd,
			any_hist == 1))

all_no_hist <-
	coxph(
		Surv(time, status) ~
			percentile +
			site +
			gender,
		data = filter(
			wd,
			any_hist == 0))

# --------------------------------------------------
# Females
# --------------------------------------------------

females_hist <-
	coxph(
		Surv(time, status) ~
			percentile +
			site,
		data = filter(
			wd,
			gender == "Female",
			any_hist == 1))

females_no_hist <-
	coxph(
		Surv(time, status) ~
			percentile +
			site,
		data = filter(
			wd,
			gender == "Female",
			any_hist == 0))

# --------------------------------------------------
# Males
# --------------------------------------------------

males_hist <-
	coxph(
		Surv(time, status) ~
			percentile +
			site,
		data = filter(
			wd,
			gender == "Male",
			any_hist == 1))

males_no_hist <-
	coxph(
		Surv(time, status) ~
			percentile +
			site,
		data = filter(
			wd,
			gender == "Male",
			any_hist == 0))

# --------------------------------------------------
# Final table
# --------------------------------------------------

final <-
	tbl_merge(
		tbls = list(
			tbl_regression(
				all_hist,
				exponentiate = TRUE,
				conf.int = TRUE),

			tbl_regression(
				all_no_hist,
				exponentiate = TRUE,
				conf.int = TRUE),

			tbl_regression(
				females_hist,
				exponentiate = TRUE,
				conf.int = TRUE),

			tbl_regression(
				females_no_hist,
				exponentiate = TRUE,
				conf.int = TRUE),

			tbl_regression(
				males_hist,
				exponentiate = TRUE,
				conf.int = TRUE),

			tbl_regression(
				males_no_hist,
				exponentiate = TRUE,
				conf.int = TRUE)),

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
	"Supplementary Table S21 - Hazard ratio by gender (family history only)" = final,
	path = "sup_tab21.docx")
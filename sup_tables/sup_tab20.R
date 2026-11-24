# Main Cox proportional hazards model
# Supplementary Table S20

pacman::p_load(
	dplyr,
	survival,
	flextable,
	gtsummary)

source(
	"C:/Users/Cass/Documents/work/ADHD_paper/other_files/survival_object.R")

# --------------------------------------------------
# Working data
# --------------------------------------------------
wd <-
	survival_data %>%
	inner_join(
		select(hist, IID, any_hist),
		by = "IID") %>%
	select(-IID)

wd_females <-
	survival_females %>%
	inner_join(
		select(hist, IID, any_hist),
		by = "IID") %>%
	select(-IID)

wd_males <-
	survival_males %>%
	inner_join(
		select(hist, IID, any_hist),
		by = "IID") %>%
	select(-IID)

# --------------------------------------------------
# Cox proportional hazards models
# --------------------------------------------------

all <-
	coxph(
		Surv(time, status) ~
			percentile +
			gender +
			site +
			any_hist,
		data = wd)

females <-
	coxph(
		Surv(time, status) ~
			percentile +
			any_hist,
		data = wd_females)

males <-
	coxph(
		Surv(time, status) ~
			percentile +
			site +
			any_hist,
		data = wd_males)

# --------------------------------------------------
# Proportional hazards assumption
# --------------------------------------------------

ph_all <- cox.zph(all)
ph_fem <- cox.zph(females)
ph_man <- cox.zph(males)

print(ph_all)
print(ph_fem)
print(ph_man)

# --------------------------------------------------
# Final table
# --------------------------------------------------

final <-
	tbl_merge(
		tbls = list(
			tbl_regression(
				all,
				exponentiate = TRUE,
				conf.int = TRUE),

			tbl_regression(
				females,
				exponentiate = TRUE,
				conf.int = TRUE),

			tbl_regression(
				males,
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
	"Supplementary Table S20 - Main Cox proportional hazards model" = final,
	path = "sup_tab20.docx")
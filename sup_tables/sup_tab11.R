# Proportional hazards assumption for Cox model by sex
# Schoenfeld residual test

pacman::p_load(
	dplyr,
	survival,
	flextable)

source(
	"C:/Users/cassi/Documents/work/ADHD_paper/other_files/survival_object.R")

# Working data
wd <-
	survival_data %>%
	inner_join(
		select(hist, IID, any_hist),
		by = "IID")

# --------------------------------------------------
# Females
# --------------------------------------------------

cox_fem <-
	coxph(
		Surv(time, status) ~
			strata(percentile) +
			any_hist +
			site,
		data = filter(
			wd,
			gender == "Female"))

ph_fem <- cox.zph(cox_fem)

# --------------------------------------------------
# Males
# --------------------------------------------------

cox_man <-
	coxph(
		Surv(time, status) ~
			strata(percentile) +
			any_hist +
			site,
		data = filter(
			wd,
			gender == "Male"))

ph_man <- cox.zph(cox_man)

# --------------------------------------------------
# Tables
# --------------------------------------------------

tab_fem <-
	data.frame(
		Variable = rownames(ph_fem$table),
		Chisq_F = ph_fem$table[, "chisq"],
		P_F = ph_fem$table[, "p"])

tab_man <-
	data.frame(
		Variable = rownames(ph_man$table),
		Chisq_M = ph_man$table[, "chisq"],
		P_M = ph_man$table[, "p"])

final <-
	left_join(
		tab_fem,
		tab_man,
		by = "Variable") %>%
	mutate(
		Chisq_F = sprintf(
			"%.2f",
			Chisq_F),

		P_F = ifelse(
			P_F < 0.001,
			"<0.001",
			sprintf("%.3f", P_F)),

		Chisq_M = sprintf(
			"%.2f",
			Chisq_M),

		P_M = ifelse(
			P_M < 0.001,
			"<0.001",
			sprintf("%.3f", P_M))) %>%
	flextable() %>%
	autofit()

save_as_docx(
	"Supplementary Table S11" = final,
	path = "sup_tab11.docx")
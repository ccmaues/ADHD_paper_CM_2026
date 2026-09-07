# Proportional hazards assumption for Cox model
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
		by = "IID") %>%
	select(-IID)

# Cox model
cox <-
	coxph(
		Surv(time, status) ~
			strata(percentile) +
			any_hist +
			gender +
			site,
		data = wd)

# Proportional hazards assumption
ph_test <- cox.zph(cox)

# Final table
final <-
	data.frame(
		Variable = rownames(ph_test$table),
		ChiSq = ph_test$table[, "chisq"],
		DF = ph_test$table[, "df"],
		P = ph_test$table[, "p"]) %>%
	mutate(
		ChiSq = sprintf("%.2f", ChiSq),
		P = ifelse(
			P < 0.001,
			"<0.001",
			sprintf("%.3f", P))) %>%
	flextable() %>%
	autofit()

save_as_docx(
	"Supplementary Table S10" = final,
	path = "sup_tab10.docx")
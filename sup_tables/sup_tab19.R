# Survival model selection by sex
# Supplementary Table S19

pacman::p_load(
	dplyr,
	purrr,
	tibble,
	survival,
	flextable,
	scales)

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
# Model comparison function
# --------------------------------------------------

compare_cox_models <- function(data, sex) {

	temp <-
		data %>%
		filter(
			gender == sex)

	# Models
	cox1 <-
		coxph(
			Surv(time, status) ~
				strata(percentile),
			data = temp)

	cox2 <-
		coxph(
			Surv(time, status) ~
				strata(percentile) +
				site,
			data = temp)

	cox3 <-
		coxph(
			Surv(time, status) ~
				strata(percentile) +
				any_hist,
			data = temp)

	cox4 <-
		coxph(
			Surv(time, status) ~
				strata(percentile) +
				site +
				any_hist,
			data = temp)

	models <- list(
		"Stratified baseline" = cox1,
		"+ Site" = cox2,
		"+ Family history" = cox3,
		"+ Site + Family history" = cox4)

	# ----------------------------------------------
	# Model-fit statistics
	# ----------------------------------------------

	comparison_table <-
		imap_dfr(
			models,
			~ {
				s <- summary(.x)

				tibble(
					Model = .y,
					Parameters = length(coef(.x)),
					logLik = as.numeric(logLik(.x)),
					AIC = AIC(.x),
					Concordance = s$concordance[1])
			}) %>%
		mutate(
			delta_AIC =
				AIC - min(AIC))

	# ----------------------------------------------
	# Likelihood-ratio tests vs baseline
	# ----------------------------------------------

	lrt_site <-
		anova(
			cox1,
			cox2,
			test = "LRT")

	lrt_hist <-
		anova(
			cox1,
			cox3,
			test = "LRT")

	lrt_site_hist <-
		anova(
			cox1,
			cox4,
			test = "LRT")

	lrt_table <-
		tibble(
			Model = c(
				"Stratified baseline",
				"+ Site",
				"+ Family history",
				"+ Site + Family history"),

			LRT_ChiSq = c(
				NA,
				lrt_site$Chisq[2],
				lrt_hist$Chisq[2],
				lrt_site_hist$Chisq[2]),

			LRT_p = c(
				NA,
				lrt_site$`Pr(>|Chi|)`[2],
				lrt_hist$`Pr(>|Chi|)`[2],
				lrt_site_hist$`Pr(>|Chi|)`[2]))

	# ----------------------------------------------
	# Final table
	# ----------------------------------------------

	comparison_table %>%
		left_join(
			lrt_table,
			by = "Model") %>%
		relocate(Model) %>%
		mutate(
			across(
				c(
					logLik,
					AIC,
					delta_AIC,
					Concordance,
					LRT_ChiSq),
				~ round(.x, 3)),

			LRT_p =
				scales::pvalue(
					LRT_p)) %>%
		flextable() %>%
		bold(
			part = "header") %>%
		align(
			part = "all",
			align = "center") %>%
		theme_booktabs() %>%
		autofit()
}

# --------------------------------------------------
# Males
# --------------------------------------------------

final_male <-
	compare_cox_models(
		wd,
		"Male")

save_as_docx(
	"Supplementary Table S19. Model selection - Males" =
		final_male,
	path = "sup_tab19_part1.docx")

# --------------------------------------------------
# Females
# --------------------------------------------------

final_female <-
	compare_cox_models(
		wd,
		"Female")

save_as_docx(
	"Supplementary Table S19. Model selection - Females" =
		final_female,
	path = "sup_tab19_part2.docx")
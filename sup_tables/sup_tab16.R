# Figure 5 - part 2 - table
# Nagelkerke pseudo-R² by PRS quintile and wave

pacman::p_load(
	dplyr,
	tidyr,
	purrr,
	DescTools,
	flextable)

data <- readRDS(
	"C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

# --------------------------------------------------
# Proband data
# --------------------------------------------------

database <-
	data$proband_data %>%
	mutate(
		across(
			c(W0, W1, W2, W3),
			~ ifelse(.x == 2, 1, .x)))

# --------------------------------------------------
# Family history
# --------------------------------------------------

hist <-
	data$family_history %>%
	mutate(
		any_hist = if_else(
			if_any(
				starts_with("parent_"),
				~ . == 1),
			1,
			0)) %>%
	select(
		IID,
		any_hist)

# --------------------------------------------------
# PCA-adjusted PRS
# Used only to define risk quintiles
# --------------------------------------------------

all_pcs <-
	data$PCA_all_samples %>%
	inner_join(
		select(
			database,
			IID,
			PRS),
		by = "IID") %>%
	select(
		IID,
		PC1,
		PC2,
		PC3,
		PC4,
		PRS)

prs_for_risk <-
	all_pcs %>%
	mutate(
		PRS_for_risk = residuals(
			glm(
				PRS ~ PC1 + PC2 + PC3 + PC4,
				family = "gaussian",
				data = all_pcs))) %>%
	select(
		IID,
		PRS_for_risk)

database <-
	database %>%
	inner_join(
		prs_for_risk,
		by = "IID") %>%
	mutate(
		risk = ntile(
			PRS_for_risk,
			5)) %>%
	inner_join(
		hist,
		by = "IID")

# --------------------------------------------------
# Pseudo-R² function
# --------------------------------------------------

calc_r2 <- function(
	df,
	wave,
	risk_group,
	sex = NULL) {

	age_var <- paste0(
		"age_",
		wave)

	temp <-
		df %>%
		filter(
			risk == risk_group)

	if (!is.null(sex)) {

		temp <-
			temp %>%
			filter(
				gender == sex)

		formula <-
			as.formula(
				paste0(
					wave,
					" ~ PRS + any_hist + ",
					age_var))

	} else {

		formula <-
			as.formula(
				paste0(
					wave,
					" ~ PRS + gender + any_hist + ",
					age_var))
	}

	model <-
		glm(
			formula,
			family = binomial,
			data = temp)

	tibble(
		risk = risk_group,
		wave = wave,
		R2 = as.numeric(
			PseudoR2(
				model,
				which = "Nagelkerke")))
}

# --------------------------------------------------
# Analysis
# --------------------------------------------------

waves <-
	c(
		"W0",
		"W1",
		"W2",
		"W3")

risk_groups <- 1:5

# All samples
tab_all <-
	map_dfr(
		waves,
		~ map_dfr(
			risk_groups,
			\(r)
			calc_r2(
				database,
				wave = .x,
				risk_group = r))) %>%
	mutate(
		cohort = "All")

# Males
tab_male <-
	map_dfr(
		waves,
		~ map_dfr(
			risk_groups,
			\(r)
			calc_r2(
				database,
				wave = .x,
				risk_group = r,
				sex = "Male"))) %>%
	mutate(
		cohort = "Male")

# Females
tab_female <-
	map_dfr(
		waves,
		~ map_dfr(
			risk_groups,
			\(r)
			calc_r2(
				database,
				wave = .x,
				risk_group = r,
				sex = "Female"))) %>%
	mutate(
		cohort = "Female")

# --------------------------------------------------
# Final table
# --------------------------------------------------

final <-
	bind_rows(
		tab_all,
		tab_male,
		tab_female) %>%
	mutate(
		wave = factor(
			wave,
			levels = waves),

		risk = factor(
			risk,
			levels = 1:5),

		R2 = round(
			R2 * 100,
			2)) %>%
	pivot_wider(
		names_from = risk,
		values_from = R2,
		names_prefix = "Q") %>%
	relocate(
		cohort,
		wave) %>%
	flextable() %>%
	bold(
		part = "header") %>%
	align(
		part = "all",
		align = "center") %>%
	autofit()

save_as_docx(
	"Supplementary Table S16" = final,
	path = "sup_tab16.docx")
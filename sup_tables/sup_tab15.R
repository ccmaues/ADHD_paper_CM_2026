# Odds ratio for low, intermediate, and high PRS groups
# all samples, females, and males

pacman::p_load(
	dplyr,
	tidyr,
	purrr,
	broom,
	flextable)

options(scipen = 999)

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

females <-
	database %>%
	filter(gender == "Female")

males <-
	database %>%
	filter(gender == "Male")

# --------------------------------------------------
# PCA
# --------------------------------------------------

all_pcs <-
	data$PCA_all_samples %>%
	inner_join(
		select(
			database,
			IID,
			PRS,
			gender),
		by = "IID")

fem_pcs <-
	data$PCA_by_sex %>%
	inner_join(
		select(
			females,
			IID,
			PRS),
		by = "IID")

man_pcs <-
	data$PCA_by_sex %>%
	inner_join(
		select(
			males,
			IID,
			PRS),
		by = "IID")

# --------------------------------------------------
# PRS correction
# --------------------------------------------------

prs_all <-
	all_pcs %>%
	mutate(
		PRS_adj = residuals(
			glm(
				PRS ~ PC1 + PC2 + PC3 + PC4 + gender,
				family = "gaussian",
				data = all_pcs))) %>%
	select(
		IID,
		PRS = PRS_adj)

prs_fem <-
	fem_pcs %>%
	mutate(
		PRS_adj = residuals(
			glm(
				PRS ~ PC1 + PC2 + PC3 + PC4,
				family = "gaussian",
				data = fem_pcs))) %>%
	select(
		IID,
		PRS = PRS_adj)

prs_man <-
	man_pcs %>%
	mutate(
		PRS_adj = residuals(
			glm(
				PRS ~ PC1 + PC2 + PC3 + PC4,
				family = "gaussian",
				data = man_pcs))) %>%
	select(
		IID,
		PRS = PRS_adj)

database <-
	database %>%
	select(-PRS) %>%
	inner_join(
		prs_all,
		by = "IID")

females <-
	females %>%
	select(-PRS) %>%
	inner_join(
		prs_fem,
		by = "IID")

males <-
	males %>%
	select(-PRS) %>%
	inner_join(
		prs_man,
		by = "IID")

# --------------------------------------------------
# Long-format data
# --------------------------------------------------

make_long <- function(df) {

	df %>%
		select(
			W0, W1, W2, W3,
			PRS) %>%
		mutate(
			risk = ntile(PRS, 100)) %>%
		pivot_longer(
			cols = c(W0, W1, W2, W3),
			names_to = "wave",
			values_to = "diagnosis") %>%
		mutate(
			wave = factor(
				wave,
				levels = c(
					"W0",
					"W1",
					"W2",
					"W3")),

			type1 = case_when(
				risk >= 90 ~ "high",
				risk <= 10 ~ "low",
				TRUE ~ "else"),

			type2 = ifelse(
				risk <= 10,
				"low",
				"else"),

			type3 = ifelse(
				risk >= 90,
				"high",
				"else"),

			type1 = factor(
				type1,
				levels = c(
					"low",
					"else",
					"high")),

			type2 = factor(
				type2,
				levels = c(
					"else",
					"low")),

			type3 = factor(
				type3,
				levels = c(
					"else",
					"high")))
}

database_long <- make_long(database)
females_long <- make_long(females)
males_long <- make_long(males)

# --------------------------------------------------
# Interaction tests
# all samples only, as in original script
# --------------------------------------------------

# High vs Low
full_model <-
	glm(
		diagnosis ~ type1 * wave,
		family = binomial,
		data = database_long)

reduced_model <-
	glm(
		diagnosis ~ type1,
		family = binomial,
		data = database_long)

anova(
	reduced_model,
	full_model,
	test = "LRT")

# Low vs Intermediate
full_model <-
	glm(
		diagnosis ~ type2 * wave,
		family = binomial,
		data = database_long)

reduced_model <-
	glm(
		diagnosis ~ type2,
		family = binomial,
		data = database_long)

anova(
	reduced_model,
	full_model,
	test = "LRT")

# High vs Intermediate
full_model <-
	glm(
		diagnosis ~ type3 * wave,
		family = binomial,
		data = database_long)

reduced_model <-
	glm(
		diagnosis ~ type3,
		family = binomial,
		data = database_long)

anova(
	reduced_model,
	full_model,
	test = "LRT")

# --------------------------------------------------
# Wave-specific ORs
# --------------------------------------------------

get_or_results <- function(df, strata_name) {

	# High vs Low
	type1_results <-
		df %>%
		filter(
			type1 %in% c(
				"high",
				"low")) %>%
		group_by(wave) %>%
		group_modify(
			~ tidy(
				glm(
					diagnosis ~ type1,
					family = binomial,
					data = .x),
				exponentiate = TRUE,
				conf.int = TRUE)) %>%
		filter(term == "type1high") %>%
		transmute(
			wave,
			OR = estimate,
			CI_low = conf.low,
			CI_high = conf.high,
			p.value,
			comparison = "90th vs. 10th",
			strata = strata_name)

	# Low vs Intermediate
	type2_results <-
		df %>%
		group_by(wave) %>%
		group_modify(
			~ tidy(
				glm(
					diagnosis ~ type2,
					family = binomial,
					data = .x),
				exponentiate = TRUE,
				conf.int = TRUE)) %>%
		filter(term == "type2low") %>%
		transmute(
			wave,
			OR = estimate,
			CI_low = conf.low,
			CI_high = conf.high,
			p.value,
			comparison = "10th vs. Else",
			strata = strata_name)

	# High vs Intermediate
	type3_results <-
		df %>%
		group_by(wave) %>%
		group_modify(
			~ tidy(
				glm(
					diagnosis ~ type3,
					family = binomial,
					data = .x),
				exponentiate = TRUE,
				conf.int = TRUE)) %>%
		filter(term == "type3high") %>%
		transmute(
			wave,
			OR = estimate,
			CI_low = conf.low,
			CI_high = conf.high,
			p.value,
			comparison = "90th vs. Else",
			strata = strata_name)

	bind_rows(
		type1_results,
		type2_results,
		type3_results)
}

fp1 <-
	get_or_results(
		database_long,
		"all")

fp2 <-
	get_or_results(
		females_long,
		"females")

fp3 <-
	get_or_results(
		males_long,
		"males")

# --------------------------------------------------
# Final table
# --------------------------------------------------

wd <-
	bind_rows(
		fp1,
		fp2,
		fp3) %>%
	mutate(
		strata = factor(
			strata,
			levels = c(
				"all",
				"males",
				"females")),

		comparison = factor(
			comparison,
			levels = c(
				"90th vs. 10th",
				"10th vs. Else",
				"90th vs. Else")))

final <-
	wd %>%
	mutate(
		estimate = sprintf(
			"%.2f (%.2f–%.2f)",
			OR,
			CI_low,
			CI_high),

		p.value = format.pval(
			p.value,
			digits = 3,
			eps = .001)) %>%

	select(
		strata,
		comparison,
		wave,
		estimate,
		p.value) %>%

	rename(
		Sample = strata,
		Comparison = comparison,
		Wave = wave,
		`OR (95% CI)` = estimate,
		`P-value` = p.value) %>%

	flextable() %>%
	bold(part = "header") %>%
	align(
		part = "all",
		align = "center") %>%
	theme_booktabs() %>%
	autofit()

save_as_docx(
	"Supplementary Table S15" = final,
	path = "sup_tab15.docx")
# Figure 5 - part 1 - table
# ADHD prevalence by PRS quintile

pacman::p_load(
	dplyr,
	tidyr,
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
# PCA-adjusted PRS
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

prs_adj <-
	all_pcs %>%
	mutate(
		PRS_adj = residuals(
			glm(
				PRS ~ PC1 + PC2 + PC3 + PC4,
				family = "gaussian",
				data = all_pcs))) %>%
	select(
		IID,
		PRS = PRS_adj)

database <-
	database %>%
	select(-PRS) %>%
	inner_join(
		prs_adj,
		by = "IID")

# --------------------------------------------------
# Stratify by sex
# --------------------------------------------------

females <-
	database %>%
	filter(
		gender == "Female")

males <-
	database %>%
	filter(
		gender == "Male")

# --------------------------------------------------
# Prevalence function
# --------------------------------------------------

calc_prev <- function(
	data,
	n,
	column_name,
	wave) {

	data %>%
	mutate(
		quantile_n = ntile(
			.data[[column_name]],
			n)) %>%
	group_by(
		quantile_n) %>%
	summarise(
		prev = mean(
			.data[[wave]] == 1,
			na.rm = TRUE),
		.groups = "drop") %>%
	rename(
		ntile = quantile_n)
}

# --------------------------------------------------
# All samples
# --------------------------------------------------

p1 <-
	calc_prev(
		database,
		5,
		"PRS",
		"W0") %>%
	rename(
		W0 = prev)

p2 <-
	calc_prev(
		database,
		5,
		"PRS",
		"W1") %>%
	rename(
		W1 = prev)

p3 <-
	calc_prev(
		database,
		5,
		"PRS",
		"W2") %>%
	rename(
		W2 = prev)

p4 <-
	calc_prev(
		database,
		5,
		"PRS",
		"W3") %>%
	rename(
		W3 = prev)

# --------------------------------------------------
# Females
# --------------------------------------------------

p1_fem <-
	calc_prev(
		females,
		5,
		"PRS",
		"W0") %>%
	rename(
		W0 = prev)

p2_fem <-
	calc_prev(
		females,
		5,
		"PRS",
		"W1") %>%
	rename(
		W1 = prev)

p3_fem <-
	calc_prev(
		females,
		5,
		"PRS",
		"W2") %>%
	rename(
		W2 = prev)

p4_fem <-
	calc_prev(
		females,
		5,
		"PRS",
		"W3") %>%
	rename(
		W3 = prev)

# --------------------------------------------------
# Males
# --------------------------------------------------

p1_male <-
	calc_prev(
		males,
		5,
		"PRS",
		"W0") %>%
	rename(
		W0 = prev)

p2_male <-
	calc_prev(
		males,
		5,
		"PRS",
		"W1") %>%
	rename(
		W1 = prev)

p3_male <-
	calc_prev(
		males,
		5,
		"PRS",
		"W2") %>%
	rename(
		W2 = prev)

p4_male <-
	calc_prev(
		males,
		5,
		"PRS",
		"W3") %>%
	rename(
		W3 = prev)

# --------------------------------------------------
# Plot data preparation
# --------------------------------------------------

for_plot_overall <-
	plyr::join_all(
		list(
			p1,
			p2,
			p3,
			p4),
		by = "ntile",
		type = "inner") %>%
	pivot_longer(
		cols = starts_with("W"),
		values_to = "prevalence",
		names_to = "wave")

for_plot_female <-
	plyr::join_all(
		list(
			p1_fem,
			p2_fem,
			p3_fem,
			p4_fem),
		by = "ntile",
		type = "inner") %>%
	pivot_longer(
		cols = starts_with("W"),
		values_to = "prevalence",
		names_to = "wave")

for_plot_male <-
	plyr::join_all(
		list(
			p1_male,
			p2_male,
			p3_male,
			p4_male),
		by = "ntile",
		type = "inner") %>%
	pivot_longer(
		cols = starts_with("W"),
		values_to = "prevalence",
		names_to = "wave")

# --------------------------------------------------
# Long-format data
# --------------------------------------------------

database_long <-
	database %>%
	select(
		W0,
		W1,
		W2,
		W3,
		PRS) %>%
	mutate(
		risk = ntile(
			PRS,
			5)) %>%
	pivot_longer(
		cols = c(
			W0,
			W1,
			W2,
			W3),
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

		diagnosis = factor(
			diagnosis,
			levels = c(
				0,
				1)),

		risk = factor(
			risk,
			levels = 1:5))

females_long <-
	database %>%
	filter(
		gender == "Female") %>%
	select(
		W0,
		W1,
		W2,
		W3,
		PRS) %>%
	mutate(
		risk = ntile(
			PRS,
			5)) %>%
	pivot_longer(
		cols = c(
			W0,
			W1,
			W2,
			W3),
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

		diagnosis = factor(
			diagnosis,
			levels = c(
				0,
				1)),

		risk = factor(
			risk,
			levels = 1:5))

males_long <-
	database %>%
	filter(
		gender == "Male") %>%
	select(
		W0,
		W1,
		W2,
		W3,
		PRS) %>%
	mutate(
		risk = ntile(
			PRS,
			5)) %>%
	pivot_longer(
		cols = c(
			W0,
			W1,
			W2,
			W3),
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

		diagnosis = factor(
			diagnosis,
			levels = c(
				0,
				1)),

		risk = factor(
			risk,
			levels = 1:5))

# --------------------------------------------------
# Final table
# --------------------------------------------------

final <-
	bind_rows(

		mutate(
			plyr::join_all(
				list(
					p1,
					p2,
					p3,
					p4),
				by = "ntile"),
			cohort = "All"),

		mutate(
			plyr::join_all(
				list(
					p1_male,
					p2_male,
					p3_male,
					p4_male),
				by = "ntile"),
			cohort = "Male"),

		mutate(
			plyr::join_all(
				list(
					p1_fem,
					p2_fem,
					p3_fem,
					p4_fem),
				by = "ntile"),
			cohort = "Female")) %>%

	mutate(
		across(
			starts_with("W"),
			~ round(
				.x * 100,
				2)),

		ntile = factor(
			ntile,
			levels = 1:5,
			labels = c(
				"1st",
				"2nd",
				"3rd",
				"4th",
				"5th"))) %>%

	relocate(
		cohort,
		ntile) %>%

	mutate(
		delta = round(
			W3 - W0,
			2)) %>%

	rename(
		`Δ W0–W3` = delta) %>%

	flextable() %>%
	bold(
		part = "header") %>%
	align(
		part = "all",
		align = "center") %>%
	add_header_row(
		values = c(
			"",
			"",
			"Prevalence (%)",
			""),
		colwidths = c(
			1,
			1,
			4,
			1)) %>%
	autofit()

save_as_docx(
	"Supplementary Table S17" = final,
	path = "sup_tab17.docx")
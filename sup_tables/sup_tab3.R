# PRS table characteristics (all samples)
pacman::p_load(dplyr, flextable, officer)

data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

# Family history
hist <-
	data$family_history %>%
	mutate(
		any_hist = if_else(
			if_any(starts_with("parent_"), ~ . == 1),
			1,
			0)) %>%
	select(IID, any_hist)

# Proband data
database <-
	data$proband_data %>%
	mutate(
		across(
			c(W0, W1, W2, W3),
			~ ifelse(.x == 2, 1, .x))) %>%
	inner_join(hist, by = "IID")

females <- filter(database, gender == "Female")
males <- filter(database, gender == "Male")

# --------------------------------------------------
# PRS correction
# --------------------------------------------------

correct_prs <- function(db, pca) {

	pcs <-
		pca %>%
		inner_join(
			select(db, IID, PRS),
			by = "IID") %>%
		select(IID, PC1, PC2, PC3, PC4, PRS)

	prs_adj <-
		pcs %>%
		mutate(
			PRS_adj = residuals(
				glm(
					PRS ~ PC1 + PC2 + PC3 + PC4,
					family = "gaussian",
					data = pcs))) %>%
		select(IID, PRS = PRS_adj)

	db %>%
		select(-PRS) %>%
		inner_join(prs_adj, by = "IID")
}

wd <- list(
	all = correct_prs(database, data$PCA_all_samples),
	fem = correct_prs(females, data$PCA_by_sex),
	man = correct_prs(males, data$PCA_by_sex)) %>%

	lapply(
		.,
		mutate,
		percentile = ntile(PRS, 100),
		risk_strata = case_when(
			percentile <= 10 ~ "Low",
			percentile >= 90 ~ "High",
			TRUE ~ "Intermediate"),
		risk_strata = factor(
			risk_strata,
			levels = c(
				"Low",
				"Intermediate",
				"High")))

# --------------------------------------------------
# Table
# --------------------------------------------------

final <-
	wd$all %>%
	group_by(risk_strata) %>%
	summarise(
		N = n(),

		Female =
			sprintf(
				"%d (%.1f%%)",
				sum(gender == "Female"),
				100 * mean(gender == "Female")),

		SP =
			sprintf(
				"%d (%.1f%%)",
				sum(site == "SP"),
				100 * mean(site == "SP")),

		Family_history =
			sprintf(
				"%d (%.1f%%)",
				sum(any_hist == 1),
				100 * mean(any_hist == 1)),

		ADHD_W0 =
			sprintf(
				"%d (%.1f%%)",
				sum(W0 == 1),
				100 * mean(W0 == 1)),

		ADHD_W1 =
			sprintf(
				"%d (%.1f%%)",
				sum(W1 == 1),
				100 * mean(W1 == 1)),

		ADHD_W2 =
			sprintf(
				"%d (%.1f%%)",
				sum(W2 == 1),
				100 * mean(W2 == 1)),

		ADHD_W3 =
			sprintf(
				"%d (%.1f%%)",
				sum(W3 == 1),
				100 * mean(W3 == 1)),

		Age_W0 =
			sprintf(
				"%.1f [%.1f, %.1f]",
				median(age_W0),
				quantile(age_W0, .25),
				quantile(age_W0, .75)),

		Age_W1 =
			sprintf(
				"%.1f [%.1f, %.1f]",
				median(age_W1),
				quantile(age_W1, .25),
				quantile(age_W1, .75)),

		Age_W2 =
			sprintf(
				"%.1f [%.1f, %.1f]",
				median(age_W2),
				quantile(age_W2, .25),
				quantile(age_W2, .75)),

		Age_W3 =
			sprintf(
				"%.1f [%.1f, %.1f]",
				median(age_W3),
				quantile(age_W3, .25),
				quantile(age_W3, .75)),

		.groups = "drop") %>%

	flextable() %>%
	bold(part = "header") %>%
	autofit() %>%

	add_footer_lines(
		values = c(
			"Note. Continuous variables are reported as median [Q1, Q3].",
			"Categorical variables are reported as n (%).",
			"PRS strata were defined as Low (1st–10th percentile), Intermediate (11th–89th percentile), and High (90th–100th percentile)."))

save_as_docx(
	"Supplementary Table S3. Cohort characteristics by PRS strata." = final,
	path = "sup_tab3.docx")
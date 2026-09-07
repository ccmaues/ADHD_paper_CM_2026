# Evaluation of generated PRS in this study
# per wave (all, males, females)
pacman::p_load(
	dplyr, tidyr, nsROC, PRROC,
	DescTools, purrr, gt)

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
	inner_join(hist, by = "IID") %>%
	mutate(
		across(
			c(W0, W1, W2, W3),
			~ ifelse(.x == 2, 1, .x)))

# Datasets
datasets <- list(
	all =
		database %>%
		inner_join(
			data$PCA_all_samples,
			by = "IID"),

	females =
		database %>%
		inner_join(
			data$PCA_by_sex,
			by = "IID") %>%
		filter(gender == "Female"),

	males =
		database %>%
		inner_join(
			data$PCA_by_sex,
			by = "IID") %>%
		filter(gender == "Male"))

# --------------------------------------------------
# Functions
# --------------------------------------------------

# PRS correction
correct_prs <- function(df) {

	new_prs <-
		residuals(
			glm(
				PRS ~ PC1 + PC2 + PC3 + PC4,
				family = "gaussian",
				data = df))

	cbind(
		select(df, -PRS),
		PRS = new_prs)
}

# Evaluation
get_metrics <- function(wave, df, sex_adjust = TRUE) {

	age_var <- paste0("age_", wave)

	formula <-
		if (sex_adjust) {
			as.formula(
				paste0(
					wave,
					" ~ PRS + any_hist + gender + ",
					age_var))
		} else {
			as.formula(
				paste0(
					wave,
					" ~ PRS + any_hist + ",
					age_var))
		}

	model <-
		glm(
			formula,
			family = "binomial",
			data = df)

	r2 <-
		as.numeric(
			PseudoR2(
				model,
				which = "Nagelkerke"))

	auroc <-
		as.numeric(
			gROC(
				X = df$PRS,
				D = df[[wave]],
				pvac.auc = TRUE,
				side = "auto")$auc)

	aucpr <-
		as.numeric(
			pr.curve(
				scores.class0 = df$PRS,
				weights.class0 = df[[wave]],
				curve = TRUE,
				sorted = FALSE,
				max.compute = TRUE,
				min.compute = TRUE,
				rand.compute = TRUE)$auc.integral)

	tibble(
		wave = wave,
		predictor = c("R2", "AUROC", "AUCPR"),
		value = c(r2, auroc, aucpr))
}

# --------------------------------------------------
# Analysis
# --------------------------------------------------

waves <- c("W0", "W1", "W2", "W3")

tabs <- map(datasets, correct_prs)

results <-
	imap_dfr(
		tabs,
		~ map_dfr(
			waves,
			get_metrics,
			df = .x,
			sex_adjust = (.y == "all")) %>%
			mutate(subset = .y)) %>%
	mutate(
		predictor = recode(
			predictor,
			"R2" = "R²")) %>%
	pivot_wider(
		names_from = c(subset, predictor),
		values_from = value)

# --------------------------------------------------
# Table
# --------------------------------------------------

final <-
	results %>%
	gt(rowname_col = "wave") %>%

	tab_spanner(
		label = "All",
		columns = starts_with("all")) %>%

	tab_spanner(
		label = "Females",
		columns = starts_with("females")) %>%

	tab_spanner(
		label = "Males",
		columns = starts_with("males")) %>%

	cols_label(
		`all_R²` = "R²",
		all_AUROC = "AUROC",
		all_AUCPR = "AUCPR",

		`females_R²` = "R²",
		females_AUROC = "AUROC",
		females_AUCPR = "AUCPR",

		`males_R²` = "R²",
		males_AUROC = "AUROC",
		males_AUCPR = "AUCPR") %>%

	fmt_number(
		columns = where(is.numeric),
		decimals = 3)

gtsave(
	data = final,
	filename = "sup_tab4.docx")
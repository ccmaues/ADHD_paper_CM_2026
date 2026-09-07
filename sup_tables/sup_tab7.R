# Odds ratio for each PRS decile per wave
pacman::p_load(dplyr, broom, flextable, gtsummary)

data <- readRDS("C:/Users/cassi/Documents/work/cass_07092026_ARTICLE.rds")

# Proband data
database <-
	data$proband_data %>%
	mutate(
		across(
			c(W0, W1, W2, W3),
			~ ifelse(.x == 2, 1, .x))) %>%
	select(
		IID,
		gender,
		W0, W1, W2, W3,
		age_W0, age_W1, age_W2, age_W3)

# PCA + raw PRS
all_pcs <-
	data$PCA_all_samples %>%
	select(
		IID,
		PC1, PC2, PC3, PC4) %>%
	inner_join(
		select(data$proband_data, IID, PRS),
		by = "IID")

# PRS correction
new_PRS <-
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

# Family history
hist <-
	data$family_history %>%
	mutate(
		any_hist = if_else(
			if_any(starts_with("parent_"), ~ . == 1),
			1,
			0)) %>%
	select(IID, any_hist)

# Working data
wd <-
	plyr::join_all(
		list(
			database,
			new_PRS,
			hist),
		by = "IID",
		type = "inner") %>%
	mutate(
		decile = ntile(PRS, 10)) %>%
	select(-PRS, -IID)

# Models
mod_w0 <-
	glm(
		W0 ~ factor(decile) + gender + any_hist + age_W0,
		family = binomial,
		data = wd) %>%
	tbl_regression(exponentiate = TRUE) %>%
	bold_p()

mod_w1 <-
	glm(
		W1 ~ factor(decile) + gender + any_hist + age_W1,
		family = binomial,
		data = wd) %>%
	tbl_regression(exponentiate = TRUE) %>%
	bold_p()

mod_w2 <-
	glm(
		W2 ~ factor(decile) + gender + any_hist + age_W2,
		family = binomial,
		data = wd) %>%
	tbl_regression(exponentiate = TRUE) %>%
	bold_p()

mod_w3 <-
	glm(
		W3 ~ factor(decile) + gender + any_hist + age_W3,
		family = binomial,
		data = wd) %>%
	tbl_regression(exponentiate = TRUE) %>%
	bold_p()

# Final table
final <-
	tbl_merge(
		tbls = list(
			mod_w0,
			mod_w1,
			mod_w2,
			mod_w3),
		tab_spanner = c(
			"**W0**",
			"**W1**",
			"**W2**",
			"**W3**")) %>%
	bold_labels() %>%
	italicize_levels() %>%
	as_flex_table()

save_as_docx(
	"Supplementary Table S7" = final,
	path = "sup_tab7.docx")